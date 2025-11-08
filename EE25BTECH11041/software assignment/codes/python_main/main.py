import numpy as np
import matplotlib.pyplot as plt
import time


def read_pgm(filename):
    with open(filename, 'rb') as f:
        header = f.readline().decode().strip()
        if header not in ['P2', 'P5']:
            raise ValueError("Unsupported PGM format (only P2 or P5 supported).")

        line = f.readline().decode()
        while line.startswith('#'):
            line = f.readline().decode()

        width, height = map(int, line.split())
        f.readline()  # skip maxval line
        if header == 'P5':
            data = np.frombuffer(f.read(width * height), dtype=np.uint8)
        else:
            data = np.loadtxt(f, dtype=np.uint8)

    return data.reshape((height, width))


def write_pgm(filename, img):
    h, w = img.shape
    img = np.clip(img, 0, 255).astype(np.uint8)
    with open(filename, 'wb') as f:
        f.write(f'P5\n{w} {h}\n255\n'.encode())
        f.write(img.tobytes())


def compress_image_svd(A, k):
    U, S, VT = np.linalg.svd(A, full_matrices=False)
    A_k = U[:, :k] @ np.diag(S[:k]) @ VT[:k, :]
    return A_k, S


def frobenius_norm(A, B):
    return np.sqrt(np.sum((A - B) ** 2))

filename = input("Enter PGM file path (e.g., ../../figs/input_pgm/input.pgm): ").strip()
img = read_pgm(filename).astype(float)
height, width = img.shape
print(f"Loaded image: {width}x{height}")

k = int(input("Enter k: "))
start = time.time()

print("Performing SVD compression...")
compressed, sigma = compress_image_svd(img, k)

compressed = np.clip(compressed, 0, 255)
write_pgm("../../figs/output_pgm/compressed_output_python.pgm", compressed)

err = frobenius_norm(img, compressed)
print(f"Frobenius norm of error matrix = {err:.6f}")
print("Compressed image saved as 'compressed_output_python.pgm'")

end = time.time()
print(f"Runtime: {end - start:.6f} seconds")

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.title("Original")
plt.imshow(img, cmap='gray')
plt.axis('off')

plt.subplot(1, 2, 2)
plt.title(f"Compressed (k={k})")
plt.imshow(compressed, cmap='gray')
plt.axis('off')

plt.tight_layout()
plt.show()
