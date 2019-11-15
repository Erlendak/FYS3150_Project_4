import matplotlib.pyplot as plt
import numpy as np
import sys

simulation = []

with open("simulation_4d.dat", 'r') as infile:
    print("Reading from "+"simulation_4d.dat")
    for i in infile:
        data = i.split()
        simulation.append(data)

simulation = np.array(simulation).astype(float).transpose()
"""
print(simulation)
T = np.linspace(2,2.3,7)
NT = 7
simulation.reshape(NL, NT, -1)
for subarray in simulation:
    a,b,c,d,e,f,g = subarray
    plt.plot(a)
    plt.show()
"""

NL = 4
steps = int(simulation.shape[1]/NL)
print(steps)
elements = simulation.shape[0]
simulation = simulation.transpose()

latice40x40 = np.zeros([steps, elements])
latice60x60 = np.zeros([steps, elements])
latice80x80 = np.zeros([steps, elements])
latice100x100 = np.zeros([steps, elements])

for i in range(steps):
    latice40x40[i] = simulation[i]
    latice60x60[i] = simulation[steps+i]
    latice80x80[i] = simulation[(2*steps) + i]
    latice100x100[i] = simulation[(3*steps)+i]

print(latice40x40)
latice40x40 = latice40x40.transpose()
latice60x60 = latice60x60.transpose()
latice80x80 = latice80x80.transpose()
latice100x100 = latice100x100.transpose()

plt.plot(latice40x40[0], latice40x40[2] , label = "40x40    <E>")#,marker ="o")
plt.plot(latice60x60[0], latice60x60[2] , label = "60x60  <E>")#,marker ="o")
plt.plot(latice80x80[0], latice80x80[2] , label = "80x80    <E>")#,marker ="o")
plt.plot(latice100x100[0], latice100x100[2] , label = "100x100  <E>")#,marker ="o")
plt.title("Energien per atom",size=17)
plt.ylabel("Energi ; ",size=15)
plt.xlabel("Temperatur ; ",size=15)
plt.grid()
plt.legend()
#plt.yscale('log')

#plt.ylim([0.15,0.20])
plt.show()


plt.plot(latice40x40[0], latice40x40[5] , label = "40x40    <|M|>")#,marker ="o")
plt.plot(latice60x60[0], latice60x60[5] , label = "60x60  <|M|>")#,marker ="o")
plt.plot(latice80x80[0], latice80x80[5] , label = "80x80    <|M|>")#,marker ="o")
plt.plot(latice100x100[0], latice100x100[5] , label = "100x100  <|M|>")#,marker ="o")
plt.title("Magnetisme",size=17)
plt.ylabel("Magnetismasjon ; ",size=15)
plt.xlabel("Temperatur ; ",size=15)
plt.grid()
plt.legend()
#plt.yscale('log')

#plt.ylim([0.15,0.20])
plt.show()

plt.plot(latice40x40[0], latice40x40[3] , label = "40x40    "+r'$C_v$')#,marker ="o")
plt.plot(latice60x60[0], latice60x60[3] , label = "60x60    "+r'$C_v$')#,marker ="o")
plt.plot(latice80x80[0], latice80x80[3] , label = "80x80    "+r'$C_v$')#,marker ="o")
plt.plot(latice100x100[0], latice100x100[3] , label = "100x100    " +r'$C_v$')#,marker ="o")
plt.title("Varmen "+r'$Cv$',size=17)
plt.ylabel("Varmekapasitet ; ",size=15)
plt.xlabel("Temperatur ; ",size=15)
plt.grid()
plt.legend()
#plt.yscale('log')

#plt.ylim([0.15,0.20])
plt.show()


plt.plot(latice40x40[0], latice40x40[6] , label = "40x40    "+r'$\chi$')#,marker ="o")
plt.plot(latice60x60[0], latice60x60[6] , label = "60x60    "+r'$\chi$')#,marker ="o")
plt.plot(latice80x80[0], latice80x80[6] , label = "80x80    "+r'$\chi$')#,marker ="o")
plt.plot(latice100x100[0], latice100x100[6] , label = "100x100    "+r'$\chi$')#,marker ="o")
plt.title("Suseptibilitet"+r'$\chi$',size=17)
plt.ylabel("Suseptibilitet ; ",size=15)
plt.xlabel("Temperatur ; ",size=15)
plt.grid()
plt.legend()
#plt.yscale('log')

#plt.ylim([0.15,0.20])
plt.show()
