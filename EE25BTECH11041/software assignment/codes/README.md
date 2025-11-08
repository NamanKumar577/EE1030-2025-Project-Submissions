# SVD-Based Image Compression (C and Python)

## Overview
This project demonstrates image compression using Singular Value Decomposition (SVD) in both C and Python.
It converts grayscale images into a matrix form, applies SVD to decompose the image, and reconstructs it using only the top k singular values.
The value of k controls the trade-off between image quality and compression efficiency.

## Files Included

### 1. main.c
A C implementation of SVD-based image compression using the Jacobi method.

**Key features:**
- Reads both JPEG and PNG files (converted internally to grayscale PGM).
- Computes SVD using a one-sided Jacobi rotation algorithm.
- Reconstructs the compressed image using the top k singular values.
- Writes output in PGM and JPEG format.
- Prints runtime and Frobenius norm (error) for quality comparison.

**Compile and Run:**
```bash
gcc main.c -lpng -ljpeg -lm -o image_compress
./image_compress
```

**Input Example:**
```
Enter image is jpeg(1) or png(0): 1
Enter JPG path (e.g., ../figs/image.jpg): ../figs/globe.jpg
Enter values of k: 50
```

**Output Files:**
- ../input.pgm – Intermediate grayscale image.
- ../figs/compressed_output.pgm – Compressed grayscale image.
- ../figs/compressed_output.jpg – Final compressed JPEG output.

### 2. main.py
A Python implementation of the same concept using NumPy’s built-in SVD function.

**Key features:**
- Reads grayscale .pgm images.
- Performs SVD decomposition and reconstructs the image using top k singular values.
- Displays both original and compressed images using Matplotlib.
- Measures runtime and reconstruction error using the Frobenius norm.

**Run the Script:**
```bash
python3 main.py
```

**Input Example:**
```
Enter PGM file path (e.g., ../../figs/input_pgm/input.pgm): ../../figs/input_pgm/globe.pgm
Enter k: 50
```

**Output Files:**
- compressed_output_python.pgm – Compressed grayscale image.

## How It Works

1. Image to Matrix Conversion:
   The image is represented as a 2D matrix where each element corresponds to a grayscale intensity value (0–255).

2. SVD Decomposition:
   The image matrix A is decomposed as:
   A = U Σ Vᵀ
   where U and V are orthogonal matrices and Σ contains singular values.

3. Truncated Reconstruction:
   Only the top k singular values are used to reconstruct the image:
   Aₖ = Uₖ Σₖ Vₖᵀ

4. Compression-Quality Trade-off:
   Smaller k gives higher compression but lower visual quality, while larger k preserves more details but increases file size.

## Example Results

| k | C Runtime (s) | Python Runtime (s) | C Error | Python Error |
|---|----------------|--------------------|----------|---------------|
| 5  | 115.26 | 2.14 | 20501.10 | 20410.96 |
| 25 | 112.26 | 2.14 | 9355.04  | 9162.20  |
| 50 | 106.07 | 2.19 | 6124.03  | 5837.31  |

As k increases, the Frobenius error decreases and the reconstructed image quality improves.

## Visualization
The Python script automatically displays:
- The original grayscale image.
- The compressed image for the selected k.

You can also compare the results from both implementations using your saved output images.

## Dependencies

### For C
- libjpeg-dev
- libpng-dev
- math.h

Install with:
```bash
sudo apt install libjpeg-dev libpng-dev
```

### For Python
```bash
pip install numpy matplotlib
```

## Notes
- The C code uses the Jacobi SVD method for educational clarity but runs slower.
- The Python code uses optimized NumPy routines for faster performance.
- Both versions use the Frobenius norm to measure reconstruction accuracy.

## Final
Developed for demonstrating matrix decomposition and image compression concepts using SVD in both low-level (C) and high-level (Python) programming environments.
