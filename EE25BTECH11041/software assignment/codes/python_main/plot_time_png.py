import numpy as np

import matplotlib.pyplot as plt

k = np.array([1,5,10,25,50,100,125,150,200])
C = np.array([274.25,269.09,267.75,271.69,280.33,294.7,276.97,290.01,282.84])
P = np.array([4.17,4.244,4.14,4.1,4.32,4.2,4.4,4.22,4.42])

plt.figure(figsize=(8,6))
plt.plot(k,C,label ="c code" )
plt.plot(k,P,label="python code")
plt.grid()
plt.legend()
plt.xlabel("Values of k")
plt.ylabel("Time taken to compress")
plt.savefig("plot_time_png.png",dpi=300)
plt.show()

