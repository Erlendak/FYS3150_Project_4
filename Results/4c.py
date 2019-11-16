import matplotlib.pyplot as plt
import numpy as np
import sys
#result = []
temp1_20x20_cold = []
temp24_20x20_cold = []
temp1_20x20_random = []
temp24_20x20_random = []

with open("20x20_cold_start_temp_2_4.dat", 'r') as infile:
    print("Reading from "+"20x20_cold_start_temp_2_4.dat")
    for i in infile:
        data = i.split()
        temp24_20x20_cold.append(data)

with open("20x20_cold_start_temp_1.dat", 'r') as infile:
    print("Reading from "+"20x20_cold_start_temp_1.dat")
    for i in infile:
        data = i.split()
        temp1_20x20_cold.append(data)

with open("20x20_random_start_temp_2_4.dat", 'r') as infile:
    print("Reading from "+"20x20_random_start_temp_2_4.dat")
    for i in infile:
        data = i.split()
        temp24_20x20_random.append(data)

with open("20x20_random_start_temp_1.dat", 'r') as infile:
    print("Reading from "+"20x20_random_start_temp_1.dat")
    for i in infile:
        data = i.split()
        temp1_20x20_random.append(data)

temp1_20x20_cold = np.array(temp1_20x20_cold).astype(float).transpose()
temp24_20x20_cold = np.array(temp24_20x20_cold).astype(float).transpose()
temp1_20x20_random = np.array(temp1_20x20_random).astype(float).transpose()
temp24_20x20_random = np.array(temp24_20x20_random).astype(float).transpose()

"""
maxtid = (max(gauss_legendre[:,2][-1], gauss_laguerre[:,2][-1],brute_force_monte_carlo
[:,2][-1], importance_sampling_monte_carlo [:,2][-1],
parallel_brute_force_monte_carlo[:,2][-1], parallel_importance_sampling_monte_carlo[:,2][-1]))
tidspan = np.linspace(0,float(maxtid),4)
_a = 5*(np.pi**2)/(16**2)
analytisk = np.zeros(4) +_a
"""

plt.plot(temp1_20x20_cold[0], temp1_20x20_cold[3] , label = "Temperatur ; 1\n<M>",marker ="o")
plt.plot(temp1_20x20_cold[0], temp1_20x20_cold[4] , label = "Temperatur ; 1\n<|M|>", marker ="x")
plt.plot(temp24_20x20_cold[0], temp24_20x20_cold[3] , label = "Temperatur ; 2.4\n<M>",marker ="o")
plt.plot(temp24_20x20_cold[0], temp24_20x20_cold[4] , label = "Temperatur ; 2.4\n<|M|>", marker ="x")
plt.title("Magnetisasjon ved kald start\n <M> vs <|M|>",size=17)
plt.ylabel("Magnetisasjon ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
#plt.ylim(-1,1.5)
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()


plt.plot(temp1_20x20_random[0], temp1_20x20_random[3] , label = "Temperatur ; 1\n<M>")#,marker ="x"_random)
plt.plot(temp1_20x20_random[0], temp1_20x20_random[4] , label = "Temperatur ; 1\n<|M|>", marker ="x")
plt.plot(temp24_20x20_random[0], temp24_20x20_random[3] , label = "Temperatur ; 2.4\n<M>")#,marker ="x")
plt.plot(temp24_20x20_random[0], temp24_20x20_random[4] , label = "Temperatur ; 2.4\n<|M|>", marker ="x")
plt.title("Magnetisasjon ved tilfeldig start\n <M> vs <|M|>",size=17)
plt.ylabel("Magnetisasjon ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
#plt.ylim(-1,1.5)
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()

plt.plot(temp1_20x20_cold[0], temp1_20x20_cold[1] , label = "Kald start,    t = 1    <E>",marker ="o")
plt.plot(temp24_20x20_cold[0], temp24_20x20_cold[1] , label = "Kald start,    t = 2.4  <E>",marker ="o")
plt.plot(temp1_20x20_random[0], temp1_20x20_random[1] , label = "Tilfeldig start,    t = 1    <E>",marker ="x")
plt.plot(temp24_20x20_random[0], temp24_20x20_random[1] , label = "Tilfeldig start,    t = 2.4  <E>",marker ="x")
plt.title("Gjennomsnittlig Energi\nKald start og tilfeldig start.",size=17)
plt.ylabel("Energi ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()


plt.subplot(1,2,1)
plt.plot(temp1_20x20_cold[0], (1/1)*temp1_20x20_cold[2] , label = "Kald start,    t = 1 "+ r'$C_v $',marker = 6)
plt.plot(temp24_20x20_cold[0],(1/(1*(2.4)**2))* temp24_20x20_cold[2] , label = "Kald start,    t = 2.4 "+r'$C_v$',marker =6)
plt.plot(temp1_20x20_random[0], (1/1)*temp1_20x20_random[2] , label = "Tilfeldig start,    t = 1 "+ r'$C_v $' , marker = 7)
plt.plot(temp24_20x20_random[0],(1/(1*(2.4)**2))* temp24_20x20_random[2] , label = "Tilfeldig start,    t = 2.4 "+r'$C_v$',marker =7)
plt.title("Varmen "r'$Cv$'+"\nKald start og tilfeldig start.",size=17)# Vi trenger kald start.
plt.ylabel("Energi utbytte ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])

plt.subplot(1,2,2)
plt.plot(temp1_20x20_cold[0], (1/1)*temp1_20x20_cold[5] , label = "Kald start,    t = 1 "+ r'$C_v $',marker = 6)
plt.plot(temp24_20x20_cold[0],(1/(1*(2.4)))* temp24_20x20_cold[5] , label = "Kald start,    t = 2.4 "+r'$C_v$',marker =6)
plt.plot(temp1_20x20_random[0], (1/1)*temp1_20x20_random[5] , label = "Tilfeldig start,    t = 1 "+ r'$C_v $',marker =7)
plt.plot(temp24_20x20_random[0],(1/(1*(2.4)))* temp24_20x20_random[5] , label = "Tilfeldig start,    t = 2.4 "+r'$C_v$',marker =7)
plt.title("Suseptibilitet "+r'$\chi$'+"\nKald start og tilfeldig start.",size=17)# Vi trenger kald start.
plt.ylabel("Magnetismasjons utbytte ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()

plt.plot(temp1_20x20_cold[0], temp1_20x20_cold[6]/temp1_20x20_cold[0] , label = "Kald start,    t = 1 " ,marker = 6)
plt.plot(temp24_20x20_cold[0], temp24_20x20_cold[6]/temp24_20x20_cold[0] , label = "Kald start,    t = 2.4 ", marker = 6)
plt.plot(temp1_20x20_random[0], temp1_20x20_random[6]/temp1_20x20_random[0] , label = "Tilfeldig start,    t = 1 ",marker = 7)
plt.plot(temp24_20x20_random[0], temp24_20x20_random[6]/temp24_20x20_random[0] , label = "Tilfeldig start,    t = 2.4 ", marker = 7)
plt.title("Godkjente konfigurasjons rate\nKald start og tilfeldig start.",size=17)
plt.ylabel("Antall godkjente konfigurasjoner ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')
plt.yscale('log')

#plt.ylim([0.15,0.20])
plt.show()

plt.plot(temp1_20x20_cold[0], temp1_20x20_cold[2] , label = "Kald start,    t = 1 "+ r'$\sigma$',marker = 6)
plt.plot(temp24_20x20_cold[0], temp24_20x20_cold[2] , label = "Kald start,    t = 2.4 "+r'$\sigma$',marker =6)
plt.plot(temp1_20x20_random[0], temp1_20x20_random[2] , label = "Tilfeldig start,    t = 1 "+ r'$\sigma $' , marker = 7)
plt.plot(temp24_20x20_random[0], temp24_20x20_random[2] , label = "Tilfeldig start,    t = 2.4 "+r'$\sigma$',marker =7)
plt.title("Varianse "r'$\sigma$'+"\nKald start og tilfeldig start.",size=17)# Vi trenger kald start.
plt.ylabel(r'$Energi^2$'" ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')
plt.show()





















plt.subplot(1,2,1)
plt.plot(temp1_20x20_cold[0], temp1_20x20_cold[7] , label = r'$\Delta E_i$'+" = -j8 ",marker =7)
plt.plot(temp1_20x20_cold[0], temp1_20x20_cold[8] , label = r'$\Delta E_i$'+" = -j4 ",marker =7)
plt.plot(temp1_20x20_cold[0], temp1_20x20_cold[9] , label = r'$\Delta E_i$'+" = j0 ",marker ="x")
plt.plot(temp1_20x20_cold[0], temp1_20x20_cold[10] , label = r'$\Delta E_i$'+" = j4 ",marker = 6 )
plt.plot(temp1_20x20_cold[0], temp1_20x20_cold[11] , label = r'$\Delta E_i$'+" = j8 ",marker =6)
plt.title(r'$\Delta E $'+" ved kald start \n T = 1",size=17) #Vi trenger kald start også
plt.ylabel("Prosent ;  %",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
#plt.show()
plt.subplot(1,2,2)
plt.plot(temp1_20x20_random[0], temp1_20x20_random[7] , label = r'$\Delta E_i$'+" = -j8 ",marker =7)
plt.plot(temp1_20x20_random[0], temp1_20x20_random[8] , label = r'$\Delta E_i$'+" = -j4 ",marker =7)
plt.plot(temp1_20x20_random[0], temp1_20x20_random[9] , label = r'$\Delta E_i$'+" = j0 ",marker ="x")
plt.plot(temp1_20x20_random[0], temp1_20x20_random[10] , label = r'$\Delta E_i$'+" = j4 ",marker = 6 )
plt.plot(temp1_20x20_random[0], temp1_20x20_random[11] , label = r'$\Delta E_i$'+" = j8 ",marker =6)
plt.title(r'$\Delta E $'+" ved tilfeldig start\nT = 1",size=17)
plt.ylabel("Prosent ;  %",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()


plt.subplot(1,2,1)
plt.plot(temp24_20x20_cold[0], temp24_20x20_cold[7] , label = r'$\Delta E_i$'+" = -j8 ",marker =7)
plt.plot(temp24_20x20_cold[0], temp24_20x20_cold[8] , label = r'$\Delta E_i$'+" = -j4 ",marker =7)
plt.plot(temp24_20x20_cold[0], temp24_20x20_cold[9] , label = r'$\Delta E_i$'+" = j0 ",marker ="x")
plt.plot(temp24_20x20_cold[0], temp24_20x20_cold[10] , label = r'$\Delta E_i$'+" = j4 ",marker = 6 )
plt.plot(temp24_20x20_cold[0], temp24_20x20_cold[11] , label = r'$\Delta E_i$'+" = j8 ",marker =6)
plt.title(r'$\Delta E $'+" ved kald start \n T = 2.4",size=17) #Vi trenger kald start også
plt.ylabel("Prosent ;  %",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
#plt.show()
plt.subplot(1,2,2)
plt.plot(temp24_20x20_random[0], temp24_20x20_random[7] , label = r'$\Delta E_i$'+" = -j8 ",marker =7)
plt.plot(temp24_20x20_random[0], temp24_20x20_random[8] , label = r'$\Delta E_i$'+" = -j4 ",marker =7)
plt.plot(temp24_20x20_random[0], temp24_20x20_random[9] , label = r'$\Delta E_i$'+" = j0 ",marker ="x")
plt.plot(temp24_20x20_random[0], temp24_20x20_random[10] , label = r'$\Delta E_i$'+" = j4 ",marker = 6 )
plt.plot(temp24_20x20_random[0], temp24_20x20_random[11] , label = r'$\Delta E_i$'+" = j8 ",marker = 6)
plt.title(r'$\Delta E $'+" ved tilfeldig start\nT = 2.4",size=17)
plt.ylabel("Prosent ;  %",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()
