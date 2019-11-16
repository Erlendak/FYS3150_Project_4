import matplotlib.pyplot as plt
import numpy as np
import sys
E_dis_1_cold = []
E_dis_24_cold = []

#import matplotlib
#print(matplotlib.__version__)

with open("Pe_temp1.000000_cold.dat", 'r') as infile:
    print("Reading from "+"Pe_temp1.000000_cold.dat")
    for i in infile:
        data = i.split()
        E_dis_1_cold.append(data)



with open("Pe_temp2.400000_cold.dat", 'r') as infile:
    print("Reading from "+"Pe_temp2.400000_cold.dat")
    for i in infile:
        data = i.split()
        E_dis_24_cold.append(data)


E_dis_1_cold = np.array(E_dis_1_cold).astype(float)
skip = 10**4
E_dis_1_cold = E_dis_1_cold[skip:]
#np.histogram(E_distrub)
plt.subplot(2,2,1)
plt.hist(E_dis_1_cold,bins = 'auto',normed = True)
plt.title("Sannsynlighetsdistribusjon T = 1\nKald start",size=17)
plt.ylabel("Prosent ; 100 %",size=15)
plt.xlabel("Energi ; ",size=15)
plt.grid()

E_dis_24_cold = np.array(E_dis_24_cold).astype(float)
skip = 10**4
E_dis_24_cold = E_dis_24_cold[skip:]
#np.histogram(E_distrub)
plt.subplot(2,2,2)
plt.hist(E_dis_24_cold,bins = 'auto',normed = True)
plt.title("Sannsynlighetsdistribusjon T = 2.4\nKald start",size=17)
plt.ylabel("Prosent ; 100 %",size=15)
plt.xlabel("Energi ; ",size=15)
plt.grid()


E_dis_1_random = []
E_dis_24_random = []

#import matplotlib
#print(matplotlib.__version__)

with open("Pe_temp1.000000_random.dat", 'r') as infile:
    print("Reading from "+"Pe_temp1.000000_random.dat")
    for i in infile:
        data = i.split()
        E_dis_1_random.append(data)



with open("Pe_temp2.400000_random.dat", 'r') as infile:
    print("Reading from "+"Pe_temp2.400000_random.dat")
    for i in infile:
        data = i.split()
        E_dis_24_random.append(data)


E_dis_1_random = np.array(E_dis_1_random).astype(float)
skip = 10**4
E_dis_1_random = E_dis_1_random[skip:]
#np.histogram(E_distrub)
plt.subplot(2,2,3)
plt.hist(E_dis_1_random,bins = 'auto',normed = True)
plt.title("Sannsynlighetsdistribusjon T = 1\nTilfeldig start",size=17)
plt.ylabel("Prosent ; 100 %",size=15)
plt.xlabel("Energi ; ",size=15)
plt.grid()


E_dis_24_random = np.array(E_dis_24_random).astype(float)
skip = 10**4
E_dis_24_random = E_dis_24_random[skip:]
#np.histogram(E_distrub)
plt.subplot(2,2,4)
plt.hist(E_dis_24_random,bins = 'auto',normed = True)
plt.title("Sannsynlighetsdistribusjon T = 2.4\nTilfeldig start",size=17)
plt.ylabel("Prosent ; 100 %",size=15)
plt.xlabel("Energi ; ",size=15)
plt.grid()
plt.subplots_adjust(hspace = 0.4)
plt.show()
