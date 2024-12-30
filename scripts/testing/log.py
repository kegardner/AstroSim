import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(1,100, 100)
y = x**3

plt.plot(np.log10(x),np.log10(y))
plt.yscale("log")
plt.xscale("log")
plt.show()