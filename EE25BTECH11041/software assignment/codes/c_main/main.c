#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <jpeglib.h>
#include <time.h>
#include <png.h>

#define eigens 25
#define error 1e-8

int matrix_multiply(const unsigned char *mat1, int w1, int h1,
                    const unsigned char *mat2, int w2, int h2,
                    int *result)
{
    if (w1 != h2) return 1;

    for (int i = 0; i < h1; i++) {
        for (int j = 0; j < w2; j++) {
            int sum = 0;
            for (int k = 0; k < w1; k++) {
                sum += mat1[i * w1 + k] * mat2[k * w2 + j];
            }
            result[i * w2 + j] = sum;
        }
    }
    return 0;
}


void transpose(const unsigned char *mat, int w, int h, int *result) {
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            result[j * h + i] = mat[i * w + j];
        }
    }
}



//jpeg read
unsigned char *read_image(const char *filename, int *width, int *height) {
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *f = fopen(filename, "rb");
    if (!f) { perror("open"); return NULL; }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, f);
    jpeg_read_header(&cinfo, TRUE);
    cinfo.out_color_space = JCS_GRAYSCALE;
    jpeg_start_decompress(&cinfo);

    *width = cinfo.output_width;
    *height = cinfo.output_height;

    unsigned char *pixels = malloc((*width) * (*height));
    unsigned char *row = malloc(*width);

    while (cinfo.output_scanline < cinfo.output_height) {
        jpeg_read_scanlines(&cinfo, &row, 1);
        int y = cinfo.output_scanline - 1;
        for (int x = 0; x < *width; x++)
            pixels[y * (*width) + x] = row[x];
    }

    free(row);
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(f);

    return pixels;
}


//write jpeg
int write_jpeg(const char *filename, unsigned char *gray_data, int width, int height) {//, int quality
    FILE *outfile = fopen(filename, "wb");
    if (!outfile) {
        perror("Cannot open output JPEG file");
        return 1;
    }

    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);

    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = 1;
    cinfo.in_color_space = JCS_GRAYSCALE;

    jpeg_set_defaults(&cinfo);
    //jpeg_set_quality(&cinfo, quality, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    int row_stride = width;
    while (cinfo.next_scanline < cinfo.image_height) {
        unsigned char *row_pointer = &gray_data[cinfo.next_scanline * row_stride];
        jpeg_write_scanlines(&cinfo, &row_pointer, 1);
    }

    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);
    fclose(outfile);

    return 0;
}



//read pgm
unsigned char *read_pgm(const char *filename, int *width, int *height)
{
    FILE *f = fopen(filename, "rb");
    if (!f) {
        perror("Cannot open PGM file");
        return NULL;
    }

    char format[3];
    if (fscanf(f, "%2s", format) != 1) {
        fprintf(stderr, "Invalid PGM header\n");
        fclose(f);
        return NULL;
    }

    if (format[0] != 'P' || (format[1] != '5' && format[1] != '2')) {
        fprintf(stderr, "Unsupported PGM format (use P2 or P5)\n");
        fclose(f);
        return NULL;
    }
    int c = fgetc(f);
    while (c == '#') {
        while (fgetc(f) != '\n');
        c = fgetc(f);
    }
    ungetc(c, f);

    int maxval;
    if (fscanf(f, "%d %d %d", width, height, &maxval) != 3) {
        fprintf(stderr, "Invalid PGM metadata\n");
        fclose(f);
        return NULL;
    }
    fgetc(f); 

    unsigned char *data = malloc((*width) * (*height));
    if (!data) {
        perror("malloc");
        fclose(f);
        return NULL;
    }

    if (format[1] == '5') {
        fread(data, 1, (*width) * (*height), f);
    } else { 
        for (int i = 0; i < (*width) * (*height); i++) {
            int val;
            fscanf(f, "%d", &val);
            data[i] = (unsigned char)val;
        }
    }

    fclose(f);
    return data;
}
//write pgm
int write_pgm(const char *filename, const unsigned char *data, int width, int height)
{
    FILE *f = fopen(filename, "wb");
    if (!f) {
        perror("Cannot open output PGM file");
        return 1;
    }

    fprintf(f, "P5\n%d %d\n255\n", width, height);
    fwrite(data, 1, width * height, f);
    fclose(f);
    return 0;
}

// PNG READ
unsigned char *read_png(const char *filename, int *width, int *height) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("Cannot open PNG file");
        return NULL;
    }

    unsigned char header[8];
    fread(header, 1, 8, fp);
    if (png_sig_cmp(header, 0, 8)) {
        fprintf(stderr, "Error: %s is not a valid PNG file.\n", filename);
        fclose(fp);
        return NULL;
    }

    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) { fclose(fp); return NULL; }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) { png_destroy_read_struct(&png_ptr, NULL, NULL); fclose(fp); return NULL; }

    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        fclose(fp);
        return NULL;
    }

    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);
    png_read_info(png_ptr, info_ptr);

    *width = png_get_image_width(png_ptr, info_ptr);
    *height = png_get_image_height(png_ptr, info_ptr);
    int color_type = png_get_color_type(png_ptr, info_ptr);
    int bit_depth = png_get_bit_depth(png_ptr, info_ptr);

    if (bit_depth == 16)
        png_set_strip_16(png_ptr);

    if (color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png_ptr);

    if (color_type == PNG_COLOR_TYPE_RGB ||
        color_type == PNG_COLOR_TYPE_RGB_ALPHA)
        png_set_rgb_to_gray_fixed(png_ptr, 1, -1, -1); // Convert RGB to Gray

    if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_strip_alpha(png_ptr);

    png_read_update_info(png_ptr, info_ptr);

    unsigned char *pixels = malloc((*width) * (*height));
    png_bytep *row_pointers = malloc((*height) * sizeof(png_bytep));
    for (int y = 0; y < *height; y++)
        row_pointers[y] = pixels + y * (*width);

    png_read_image(png_ptr, row_pointers);

    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    free(row_pointers);
    fclose(fp);

    return pixels;
}

// A mxn
void jacobi_svd(double *A, double *V ,double *sig, int m, int n){
    // V= I
    for(int i=0;i<n*n;i++) V[i]=0;
    for(int i=0; i<n; i++) V[i*n+i]=1;


    for(int i=0; i<eigens ;i++){
        double max_off_diag=0;
        for(int j=0; j<n-1 ; j++){
            for(int k=j+1;k<n;k++){
                double alpha =0,beta=0,gamma=0;
                for (int l = 0; l < m; l++)
                {
                    double ap = A[n*l+j];
                    double aq = A[n*l+k];

                    alpha +=ap*ap;
                    beta += aq*aq;
                    gamma+= aq*ap;
                }

                max_off_diag = fmax(max_off_diag,fabs(gamma));
                if ((fabs(gamma) < error * sqrt(alpha * beta))|| gamma==0.0) continue;

                double tao = (beta - alpha) / (2.0 * gamma);
                double t = (tao>=0?1:-1)/(fabs(tao)+sqrt(1+tao*tao));
                double c = 1/sqrt(1+t*t);
                double s = c*t;
                
                for (int l = 0; l < m; l++)
                {
                    double ap = A[l*n+j], aq = A[l*n+k];

                    A[l*n+j]=c*ap-s*aq;
                    A[l*n+k]=s*ap+c*aq;
                }

                for (int l = 0; l < n; l++)
                {
                    double vp = V[l*n+j], vq = V[l*n+k];

                    V[l*n+j]=c*vp-s*vq;
                    V[l*n+k]=s*vp+c*vq;
                }
            }
        }


        if(max_off_diag<error) break;
    }
    for (int j = 0; j < n; j++) {
        double norm = 0;
        for (int i = 0; i < m; i++) {
            norm += A[i*n + j] * A[i*n + j];
        }
        sig[j] = sqrt(norm);
    }   
    for (int j = 0; j < n; j++) {
        if (sig[j] > 1e-12)
            for (int i = 0; i < m; i++)
                A[i*n + j] /= sig[j];
    }

    for (int i = 0; i < n - 1; i++) {
        int max_idx = i;
        for (int j = i + 1; j < n; j++) {
            if (sig[j] > sig[max_idx]) max_idx = j;
        }

        if (max_idx != i) {
            // swap σ
            double temp = sig[i];
            sig[i] = sig[max_idx];
            sig[max_idx] = temp;

            // swap U columns
            for (int r = 0; r < m; r++) {
                double tmp = A[r * n + i];
                A[r * n + i] = A[r * n + max_idx];
                A[r * n + max_idx] = tmp;
            }

            // swap V columns
            for (int r = 0; r < n; r++) {
                double tmp = V[r * n + i];
                V[r * n + i] = V[r * n + max_idx];
                V[r * n + max_idx] = tmp;
            }
        }
    }

}


//forbenius norm
double forbenius(double *A, int m, int n){
    double sum=0;
    for(int i=0; i<m;i++){
        for (int j = 0; j < n; j++)
        {
            sum+=A[i*n+j]*A[i*n+j];
        }
    }

    return sqrt(sum);
}

int main() {
    int width, height;
    char filename[256];

    clock_t start, end;
    double cpu_time_used;
    start = clock();

    int response;
    printf("Enter image is jpeg(1) or png(0): ");
    scanf("%d", &response);

    if(response){
        printf("Enter JPG path (e.g., ../../figs/name.jpg): ");
        scanf("%255s", filename); 
    
        unsigned char *jpg = read_image(filename, &width, &height);
        if (!jpg) return 1;
    
        if (write_pgm("../../figs/input_pgm/input.pgm", jpg, width, height) == 0)
            printf("Successfully converted JPEG to PGM!\n");
        else
            printf("Failed to write PGM.\n");
            free(jpg);
    }
    else{
        printf("Enter PNG path (e.g., ../../figs/name.png): ");
        scanf("%255s", filename); 
    
        unsigned char *png =  read_png(filename, &width, &height);
        if (!png) return 1;
    
        if (write_pgm("../../figs/input_pgm/input.pgm", png, width, height) == 0)
            printf("Successfully converted PNG to PGM!\n");
        else
            printf("Failed to write PGM.\n");
            free(png);
    }

        
    unsigned char *img = read_pgm("../../figs/input_pgm/input.pgm", &width, &height);
    if (!img) return 1;
    
    // making char matrix into double matrix for jacobi
    double *A = malloc(width * height * sizeof(double));
    double *error_matrix = malloc(width * height * sizeof(double));
    for (int i=0;i<height*width;i++) A[i]=img[i];


    
    // Build identity matrix
    //unsigned char *ide = malloc(width * width * sizeof(unsigned char));
    // for (int i = 0; i < width; i++) {
    //     for (int j = 0; j < width; j++) {
    //         ide[i * width + j] = (i == j) ? 1 : 0;
    //     }
    // }
    
    int minm = (height < width) ? height : width;
    double *sig = calloc(minm, sizeof(double));

    double *V = malloc(width*width*sizeof(double));
    jacobi_svd(A,V, sig, height, width );



    // now im=U , V=V, sig=sigma
    int k;
    printf("Enter values of k ");
    scanf("%d", &k);
    double *Ak = calloc(height*width, sizeof(double));
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            double sum=0;
            for (int r = 0; r < k; r++)
            {
                sum+=A[i*width+r]*sig[r]*V[j*width+r];
            }
            if(sum<0) sum=0;
            if(sum>255) sum=255;
            Ak[i*width+j]=sum;
        }
        
    }
    // round of every entery
    unsigned char *recon = malloc(height*width);
    for (int i=0;i<height*width;i++) recon[i]=(unsigned char)(Ak[i]+0.5);

    for (int i=0;i<height*width;i++) error_matrix[i]=fabs(recon[i]-img[i]);


    write_pgm("../../figs/output_pgm/compressed_output.pgm", recon, width, height);
    printf("Compressed image written.\n");
    write_jpeg("../../figs/output_jpg_png/compressed_output.jpg",recon, width,height);
    
    printf("%lf\n", forbenius(error_matrix, height,width));


    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("Runtime: %.6f seconds\n", cpu_time_used);

    free(img);
    free(A);
    free(Ak);free(V);free(recon);free(sig);

    return 0;
}

// gcc main.c -lpng -ljpeg -lm -o image
// ../../figs/globe.jpg