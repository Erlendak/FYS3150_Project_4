import matplotlib.pyplot as plt
import numpy as np
import sys




simulation = []

with open("simulation_4e.dat", 'r') as infile:
    print("Reading from "+"simulation_4e.dat")
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
#print(steps)
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

def sort_results(_ary):
    clone_ary = _ary.copy()
    _sort_list = _ary[0].argsort()
    #print(_ary.shape[1])
    #print(max(_sort_list))
    #clone_ary[:,10]
    for i in range(_ary.shape[1]):
        clone_ary[:,i] = _ary[:,(_sort_list[i]) ]
    return(clone_ary)

#print(latice40x40[:,0].argsort())
latice40x40 = sort_results(latice40x40.transpose())
#latice40x40 = sort_results(latice40x40)
#print(latice40x40[0].argsort())
#latice40x40 = np.matrix(latice40x40).sort(key=lambda row: row[0:], reverse = True)
latice60x60 = sort_results(latice60x60.transpose())
latice80x80 = sort_results(latice80x80.transpose())
latice100x100 = sort_results(latice100x100.transpose())

#print(latice40x)



plt.plot(latice40x40[0], latice40x40[2] , label = "40x40    <E>")#,marker ="o")
plt.plot(latice60x60[0], latice60x60[2] , label = "60x60  <E>")#,marker ="o")
plt.plot(latice80x80[0], latice80x80[2] , label = "80x80    <E>")#,marker ="o")
plt.plot(latice100x100[0], latice100x100[2] , label = "100x100  <E>")#,marker ="o")
plt.title("Energi Per Spin",size=17)
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
plt.title("Suseptibilitet "+r'$\chi$',size=17)
plt.ylabel("Suseptibilitet ; ",size=15)
plt.xlabel("Temperatur ; ",size=15)
plt.grid()
plt.legend()
#plt.yscale('log')

#plt.ylim([0.15,0.20])
plt.show()

T40x40 = latice40x40[0][(latice40x40[3].argmax())]
T60x60 = latice60x60[0][(latice60x60[3].argmax())]
T80x80 = latice80x80[0][(latice80x80[3].argmax())]
T100x100 = latice100x100[0][(latice100x100[3].argmax())]
print(T40x40)
print(T60x60)
print(T80x80)
print(T100x100)

x = np.linspace(0,150,500)

_y = np.array([T40x40, T60x60, T80x80, T100x100 ])
_x = np.array([1/40, 1/60, 1/80 , 1/100  ])
[a,b] = np.polyfit(_x, _y, deg = 1)
plt.subplot(2,1,1)
plt.plot(x,a*x +b , label = "L"+r'$\approx \infty$')#,marker ="o")
plt.gca().invert_xaxis()
plt.title("Fase trasisjon",size=17)
plt.ylabel("Temperatur ; ",size=15)
plt.xlabel("L ; ",size=15)
plt.grid()
plt.legend()

x = np.linspace(20,50000,500)
plt.subplot(2,1,2)
plt.plot(x,a*x +b , label = "L"+r'$\approx \infty$')#,marker ="o")
plt.gca().invert_xaxis()
plt.title("Fase trasisjon\n Zoomed out",size=17)
plt.ylabel("Temperatur ; ",size=15)
plt.xlabel("L ; ",size=15)
plt.grid()
plt.legend()
#plt.yscale('log')
plt.subplots_adjust(hspace=0.8)
#plt.ylim([0.15,0.20])
plt.show()
