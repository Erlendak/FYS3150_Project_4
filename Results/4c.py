import matplotlib.pyplot as plt
import numpy as np
import sys
result = []
#temp24_20x20 = []

with open("test.dat", 'r') as infile:
    print("Reading from "+"test.dat")
    for i in infile:
        data = i.split()
        result.append(data)

result = np.array(result).astype(float).transpose()
temp1_20x20 = np.zeros([int(result.shape[0]),int(result.shape[1]/2)])
temp24_20x20 = np.zeros(temp1_20x20.shape)




for i in range(int(result.shape[1]/2)):
    #print(result[0])
    temp1_20x20[:,i] = result[:,(i*2)]
    temp24_20x20[:,i] = result[:,(i*2)+1]
"""
maxtid = (max(gauss_legendre[:,2][-1], gauss_laguerre[:,2][-1],brute_force_monte_carlo
[:,2][-1], importance_sampling_monte_carlo [:,2][-1],
parallel_brute_force_monte_carlo[:,2][-1], parallel_importance_sampling_monte_carlo[:,2][-1]))
tidspan = np.linspace(0,float(maxtid),4)
_a = 5*(np.pi**2)/(16**2)
analytisk = np.zeros(4) +_a
"""
plt.plot(temp1_20x20[1], temp1_20x20[4] , label = "Temperatur ; 1\n<M>")#,marker ="x")
plt.plot(temp1_20x20[1], temp1_20x20[5] , label = "Temperatur ; 1\n<|M|>", marker ="x")
plt.plot(temp24_20x20[1], temp24_20x20[4] , label = "Temperatur ; 2.4\n<M>")#,marker ="x")
plt.plot(temp24_20x20[1], temp24_20x20[5] , label = "Temperatur ; 2.4\n<|M|>", marker ="x")
plt.title("Magnetisasjon\n <M> vs <|M|>",size=17)
plt.ylabel("Magnetisasjon ; ",size=15)
plt.xlabel("integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.ylim(-1,1.5)
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()

plt.plot(temp1_20x20[1], temp1_20x20[2] , label = "t = 1    <E>")#,marker ="o")
plt.plot(temp24_20x20[1], temp24_20x20[2] , label = "t = 2.4  <E>")#,marker ="o")
plt.title("Energien per atom",size=17)
plt.ylabel("Energi ; ",size=15)
plt.xlabel("integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()
