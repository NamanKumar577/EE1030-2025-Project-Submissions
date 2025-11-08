import numpy as np
import matplotlib.pyplot as plt

# Data
k = np.array([1,5,10,25,50,100,125,150,200])
C = np.array([39045.20647,20501.102,14918.67,9355.04,6124.03,3630.36,2969.53,2469.45,1775.52])
P = np.array([38989.92,20410.96,14791.37,9162.2,5837.31,3149.08,2371.01,1716.63,198.63])

# Take log of runtime (natural log)
C_log = np.log(C)
P_log = np.log(P)

# Plot log values vs k (semi-log linearization)
plt.figure(figsize=(8,6))
plt.plot(k, C_log, marker='o', label="C code (ln runtime)")
plt.plot(k, P_log, marker='s', label="Python code (ln runtime)")

plt.grid(True, linestyle='--', linewidth=0.6)
plt.legend()
plt.xlabel("Values of k")
plt.ylabel("ln(Runtime)")
plt.title("Linearized Plot of Runtime vs k (ln scale)")
plt.savefig("plot_error_jpg.png", dpi=300)
plt.show()
