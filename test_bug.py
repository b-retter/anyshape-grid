import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt


Nyso = np.array([1,2,3])
yso_x = []
yso_y = []
for N in Nyso:
    y,x = rnd.rand(N),rnd.rand(N)

    yso[0].append(x.tolist())
    yso[1].append(y.tolist())

yso = np.array(yso)
