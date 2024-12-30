import matplotlib.pyplot as plt
import numpy as np 


x = np.linspace(0,100,1000)
y = np.random.random(x.shape)

fig,axs = plt.subplots(3,2,figsize=(15,10))
print(axs)
axs[1][1].plot(x,y)
axs[1][1].set_title("Small lil tittle")
fig.suptitle("Big ol' title")
plt.show()