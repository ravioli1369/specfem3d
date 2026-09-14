import numpy as np
import matplotlib.pyplot as plt

# load data
stnm = 'XA.A2.CXZ.semd'
d_txt = np.loadtxt(f"OUTPUT_FILES.txt/{stnm}")
d_bin = np.loadtxt(f"OUTPUT_FILES.bin/{stnm}")

# cmt + force separately
d_sep = np.loadtxt(f"OUTPUT_FILES.force/{stnm}")
d_sep1 = np.loadtxt(f"OUTPUT_FILES.moment/{stnm}")
d_sum = d_sep[:,1] + d_sep1[:,1]

# plot
plt.figure(1,figsize=(14,5))
plt.plot(d_bin[:,0],d_bin[:,1],label='force + cmt with binary')
plt.plot(d_txt[:,0],d_txt[:,1],label='force + cmt with txt',ls='--')
plt.plot(d_sep[:,0],d_sum,label='separate',ls=':',color='k')
plt.legend()
plt.savefig("seismo.jpg",dpi=300)
