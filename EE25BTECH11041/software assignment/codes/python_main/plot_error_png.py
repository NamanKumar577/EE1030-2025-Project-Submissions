import numpy as np
import matplotlib.pyplot as plt

k = np.array([1,5,10,25,50,100,125,150,200])
C = np.array([39078.86,10995.77,6839.97,2706.81,1079.06,577,548.86,522.92,471.43])
P = np.array([39079.59,10992.788,6834.69,2692.668,1044.17,516.06,485,459.41,411.59])

# Take natural log of errors
C_log = np.log(C)
P_log = np.log(P)

plt.figure(figsize=(8,6))
plt.plot(k, C_log, marker='o', label="C code (ln error)")
plt.plot(k, P_log, marker='s', label="Python code (ln error)")
plt.grid(True, linestyle='--', linewidth=0.6)
plt.legend()
plt.xlabel("Values of k")
plt.ylabel("ln(Error)")
plt.title("Linearized Plot of Error vs k (ln scale)")
plt.savefig("plot_error_png.png", dpi=300)
plt.show()
