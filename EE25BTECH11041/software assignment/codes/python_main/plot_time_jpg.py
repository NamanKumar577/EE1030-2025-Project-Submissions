import numpy as np

import matplotlib.pyplot as plt

k = np.array([1,5,10,25,50,100,125,150,200])
C = np.array([118.922363,115.26,111.41,112.26,106.07,118.62,117.52,117.56,126.46])
P = np.array([2.13,2.14,2.12,2.14,2.19,2.2,2.2,2.3,2.3])

plt.figure(figsize=(8,6))
plt.plot(k,C,label ="c code" )
plt.plot(k,P,label="python code")
plt.grid()
plt.legend()
plt.xlabel("Values of k")
plt.ylabel("Time taken to compress")
plt.savefig("plot_time_jpg.png",dpi=300)
plt.show()
