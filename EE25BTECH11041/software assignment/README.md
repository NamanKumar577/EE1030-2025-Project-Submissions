# Image Compression Using Singular Value Decomposition (SVD)

## Overview
This project demonstrates image compression using Singular Value Decomposition (SVD) implemented in both C and Python.
It explores how SVD can be applied to reduce image size by reconstructing it from a limited number of singular values while maintaining acceptable visual quality.
The C version implements the Jacobi one-sided SVD algorithm, while the Python version uses NumPy’s optimized SVD function for performance comparison.

## Concept Summary

### Singular Value Decomposition
SVD decomposes any m × n matrix A into three matrices:

A = U Σ Vᵀ

- U: Orthogonal matrix of left singular vectors (column space basis)
- Σ: Diagonal matrix of singular values σ₁ ≥ σ₂ ≥ ... ≥ 0
- V: Orthogonal matrix of right singular vectors (row space basis)

In image compression, only the top k singular values and corresponding vectors are used to reconstruct an approximate version of the image:

Aₖ = Uₖ Σₖ Vₖᵀ

This reduces the amount of stored data while keeping the main visual information.

## Algorithms Used

### 1. Jacobi One-Sided SVD (C Implementation)
- Iteratively applies Givens rotations to orthogonalize matrix columns.
- Avoids forming AᵀA, improving numerical stability.
- Produces accurate singular values and vectors for medium-sized matrices.
- Computationally slower but ideal for low-level implementation and learning.

### 2. NumPy-Based SVD (Python Implementation)
- Uses numpy.linalg.svd, an optimized LAPACK-based routine.
- Much faster and scalable for large images.
- Ideal for validation and comparison with the C implementation.

## Project Structure

📦 SVD_Image_Compression
├── main.c                # Jacobi SVD in C (JPEG & PNG supported)
├── main.py               # Python implementation using NumPy SVD
├── report.pdf            # Full technical report and analysis
├── figs/                 # Contains input and output images
│   ├── input_pgm/
│   ├── output_pgm/
│   └── output_jpg_png/
└── README.md             # This documentation file

## How It Works

1. Convert Image to Matrix:
   The input image (JPEG/PNG) is converted to a grayscale .pgm file and read as a 2D matrix of pixel intensities (0–255).

2. Compute SVD:
   The image matrix A is decomposed using SVD into U, Σ, Vᵀ.

3. Reconstruct Using Top k Singular Values:
   The compressed image is reconstructed using:
   Aₖ = Uₖ Σₖ Vₖᵀ

4. Measure Quality and Runtime:
   The Frobenius norm quantifies reconstruction error:
   ||A - Aₖ||_F
   Runtime is measured to compare efficiency between the C and Python implementations.

## How to Run

### C Version
gcc main.c -lpng -ljpeg -lm -o image_compress
./image_compress

Example Input:
Enter image is jpeg(1) or png(0): 1
Enter JPG path (e.g., ../figs/image.jpg): ../figs/globe.jpg
Enter values of k: 50

Outputs:
- ../figs/output_pgm/compressed_output.pgm
- ../figs/output_jpg_png/compressed_output.jpg

### Python Version
python3 main.py

Example Input:
Enter PGM file path (e.g., ../../figs/input_pgm/input.pgm)
Enter k: 50

Outputs:
- compressed_output_python.pgm
- Displays original and reconstructed images for comparison.

## Results

| k  | C Error | Python Error | C Runtime (s) | Python Runtime (s) |
|----|----------|---------------|----------------|--------------------|
| 5  | 20501.10 | 20410.96 | 115.26 | 2.14 |
| 25 | 9355.04  | 9162.20  | 112.26 | 2.14 |
| 50 | 6124.03  | 5837.31  | 106.07 | 2.19 |
| 100| 3630.36  | 3149.08  | 118.62 | 2.12 |

### Observations
- Error decreases as k increases — better reconstruction.
- C implementation gives accurate results but slower runtime.
- Python is significantly faster due to optimized libraries.
- Both methods preserve orthogonality and stability.

## Image Quality vs k
As k increases:
- For k = 5: The image is heavily blurred and loses details.
- For k = 25: The image becomes recognizable but not sharp.
- For k = 50: The image is nearly identical to the original, with good clarity.

Thus, k determines the balance between compression and image fidelity.

## Dependencies

### For C
sudo apt install libjpeg-dev libpng-dev

### For Python
pip install numpy matplotlib

## Key Takeaways
- The Jacobi SVD algorithm, though slower, is simple and reliable for educational or experimental use.
- SVD effectively compresses images by eliminating redundant information in pixel space.
- Increasing k improves image quality at the cost of larger storage.
- Python’s NumPy-based approach offers faster computation with comparable accuracy.

## References
- G. H. Golub and C. F. Van Loan, Matrix Computations, 4th Edition.
- Gilbert Strang, Linear Algebra and Its Applications.
- Wikipedia: Singular Value Decomposition and Jacobi Method.
