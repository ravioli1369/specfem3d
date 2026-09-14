import numpy as np

np.random.seed(10)

nsta = 10
x = -1500 + np.random.rand(nsta) * 3000
y = -1500 + np.random.rand(nsta) * 3000

fio = open("DATA/STATIONS","w")
for i in range(nsta):
    fio.write(f"A{i} XA %f %f 0. 0.\n" %(y[i],x[i]))
fio.close()