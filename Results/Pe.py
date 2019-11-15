import matplotlib.pyplot as plt
import numpy as np
import sys
E_dis_1 = []
E_dis_24 = []

import matplotlib
print(matplotlib.__version__)

with open("Pe_temp1.000000.dat", 'r') as infile:
    print("Reading from "+"Pe_temp1.000000.dat")
    for i in infile:
        data = i.split()
        E_dis_1.append(data)



with open("Pe_temp2.400000.dat", 'r') as infile:
    print("Reading from "+"Pe_temp2.400000.dat")
    for i in infile:
        data = i.split()
        E_dis_24.append(data)


E_dis_1 = np.array(E_dis_1).astype(float)
skip = 10**4
E_dis_1 = E_dis_1[skip:]
#np.histogram(E_distrub)
plt.hist(E_dis_1,bins = 'auto',normed = True)
plt.title("Sansynlighets distrubisjon T = 1\n",size=17)
plt.ylabel("Prosent ; 100 %",size=15)
plt.xlabel("Energi ; ",size=15)
plt.grid()
plt.show()

E_dis_24 = np.array(E_dis_24).astype(float)
skip = 10**4
E_dis_24 = E_dis_24[skip:]
#np.histogram(E_distrub)
plt.hist(E_dis_24,bins = 'auto',normed = True)
plt.title("Sansynlighets distrubisjon T = 2.4\n",size=17)
plt.ylabel("Prosent ; 100 %",size=15)
plt.xlabel("Energi ; ",size=15)
plt.grid()
plt.show()
