import numpy as np


MY_SEED = 0x5EED

ex_count = 4
possibilities = [i for i in np.arange(-5,5,0.1)]
getem = lambda poss : np.random.Generator(np.random.PCG64(MY_SEED)).choice(poss, (4,ex_count), False)
xlist,ylist,xgal,ygal = getem(possibilities)