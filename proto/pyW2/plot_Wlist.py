import numpy as np
import matplotlib.pyplot as plt

# read Wlist.csv
ulist, Wlist = np.genfromtxt('Wlist.csv', delimiter=',', unpack=True)

# plot Wlist vs ulist
plt.plot(ulist, Wlist)
plt.xlim(0, 1.5)
plt.ylim(-3, 0)

plt.xlabel('u')
plt.ylabel('W(u)')

plt.savefig('Wlist.png')