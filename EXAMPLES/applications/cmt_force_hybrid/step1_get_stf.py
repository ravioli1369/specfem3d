import numpy as np
import matplotlib.pyplot as plt

def ricker_wavelet(t, f_c):
    """
    Compute the Ricker wavelet.

    Parameters:
        t (numpy array): Time array.
        f_c (float): Central frequency.

    Returns:
        numpy array: Ricker wavelet values.
    """
    pi2_fc2_t2 = (np.pi * f_c * t) ** 2
    wavelet = (1 - 2 * pi2_fc2_t2) * np.exp(-pi2_fc2_t2)
    return wavelet

nt = 3000
dt = 2.0e-3
t = np.arange(nt) * dt

NS_CMT = 10
NS_FORCE = 5

# stf_force
stf = ricker_wavelet(t-0.3,5)

for i in range(NS_FORCE):
    with open(f"DATA/stf.force.{i}.txt","w") as fio:
        for it in range(nt):
            fio.write("%f\n" %(stf[it]))

# cmt
stf1 = np.cumsum(stf) * dt
for i in range(NS_CMT):
    with open(f"DATA/stf.cmt.{i}.txt","w") as fio:
        for it in range(nt):
            fio.write("%f\n" %(stf1[it]))

plt.plot(t,stf)
plt.plot(t,stf1)
plt.show()