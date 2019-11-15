import matplotlib.pyplot as plt
import numpy as np
import sys
#result = []
temp1_20x20 = []
temp24_20x20 = []

with open("20x20_temp_2_4.dat", 'r') as infile:
    print("Reading from "+"20x20_temp_2_4.dat")
    for i in infile:
        data = i.split()
        temp24_20x20.append(data)

with open("20x20_temp_1.dat", 'r') as infile:
    print("Reading from "+"20x20_temp_1.dat")
    for i in infile:
        data = i.split()
        temp1_20x20.append(data)

temp1_20x20 = np.array(temp1_20x20).astype(float).transpose()
temp24_20x20 = np.array(temp24_20x20).astype(float).transpose()


"""
maxtid = (max(gauss_legendre[:,2][-1], gauss_laguerre[:,2][-1],brute_force_monte_carlo
[:,2][-1], importance_sampling_monte_carlo [:,2][-1],
parallel_brute_force_monte_carlo[:,2][-1], parallel_importance_sampling_monte_carlo[:,2][-1]))
tidspan = np.linspace(0,float(maxtid),4)
_a = 5*(np.pi**2)/(16**2)
analytisk = np.zeros(4) +_a
"""
plt.plot(temp1_20x20[0], temp1_20x20[3] , label = "Temperatur ; 1\n<M>")#,marker ="x")
plt.plot(temp1_20x20[0], temp1_20x20[4] , label = "Temperatur ; 1\n<|M|>", marker ="x")
plt.plot(temp24_20x20[0], temp24_20x20[3] , label = "Temperatur ; 2.4\n<M>")#,marker ="x")
plt.plot(temp24_20x20[0], temp24_20x20[4] , label = "Temperatur ; 2.4\n<|M|>", marker ="x")
plt.title("Magnetisasjon\n <M> vs <|M|>",size=17)
plt.ylabel("Magnetisasjon ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
#plt.ylim(-1,1.5)
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()

plt.plot(temp1_20x20[0], temp1_20x20[1] , label = "t = 1    <E>")#,marker ="o")
plt.plot(temp24_20x20[0], temp24_20x20[1] , label = "t = 2.4  <E>")#,marker ="o")
plt.title("Energien per atom",size=17)
plt.ylabel("Energi ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()


plt.plot(temp1_20x20[0], (1/1)*temp1_20x20[2] , label = "t = 1 "+ r'$C_v $')#,marker ="o")
plt.plot(temp24_20x20[0],(1/(1*(2.4)**2))* temp24_20x20[2] , label = "t = 2.4 "+r'$C_v$')#,marker ="o")
plt.title("Varmen "r'$Cv$'+"\nRandom start",size=17)# Vi trenger kald start.
plt.ylabel("Energi utbytte ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()

plt.plot(temp1_20x20[0], (1/1)*temp1_20x20[5] , label = "t = 1 "+ r'$C_v $')#,marker ="o")
plt.plot(temp24_20x20[0],(1/(1*(2.4)))* temp24_20x20[5] , label = "t = 2.4 "+r'$C_v$')#,marker ="o")
plt.title("Suseptibilitet "+r'$\chi$'+"\nRandom start",size=17)# Vi trenger kald start.
plt.ylabel("Magnetismasjons utbytte ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()

plt.plot(temp1_20x20[0], temp1_20x20[6]/temp1_20x20[0] , label = "t = 1 ")#,marker ="o")
plt.plot(temp24_20x20[0], temp24_20x20[6]/temp24_20x20[0] , label = "t = 2.4 ")#,marker ="o")
plt.title("Godkjente konfigurasjons rate\nKald start",size=17)
plt.ylabel("Antall godkjente konfigurasjoner ; ",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')
plt.yscale('log')

#plt.ylim([0.15,0.20])
plt.show()

plt.plot(temp1_20x20[0], temp1_20x20[7] , label = r'$E_i$'+" = -j8 ")#,marker ="o")
plt.plot(temp1_20x20[0], temp1_20x20[8] , label = r'$E_i$'+" = -j4 ")#,marker ="o")
plt.plot(temp1_20x20[0], temp1_20x20[9] , label = r'$E_i$'+" = j0 ")#,marker ="o")
plt.plot(temp1_20x20[0], temp1_20x20[10] , label = r'$E_i$'+" = j4 ",marker ="o")
plt.plot(temp1_20x20[0], temp1_20x20[11] , label = r'$E_i$'+" = j8 ",marker ="o")
plt.title(r'$\Delta E $'+" T = 1\nRandom start",size=17) #Vi trenger kald start også
plt.ylabel("Prosent ;  %",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()

plt.plot(temp24_20x20[0], temp24_20x20[7] , label = r'$E_i$'+" = -j8 ")#,marker ="o")
plt.plot(temp24_20x20[0], temp24_20x20[8] , label = r'$E_i$'+" = -j4 ")#,marker ="o")
plt.plot(temp24_20x20[0], temp24_20x20[9] , label = r'$E_i$'+" = j0 ")#,marker ="o")
plt.plot(temp24_20x20[0], temp24_20x20[10] , label = r'$E_i$'+" = j4 ",marker ="o")
plt.plot(temp24_20x20[0], temp24_20x20[11] , label = r'$E_i$'+" = j8 ",marker ="o")
plt.title(r'$\Delta E $'+"T = 2.4\nRandom start",size=17)
plt.ylabel("Prosent ;  %",size=15)
plt.xlabel("Integrasjonspoeng ; ",size=15)
plt.grid()
plt.legend()
plt.xscale('log')

#plt.ylim([0.15,0.20])
plt.show()
