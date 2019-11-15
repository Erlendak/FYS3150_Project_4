import matplotlib.pyplot as plt
import numpy as np
import sys

simulation = []

with open("simulation_4d.dat", 'r') as infile:
    print("Reading from "+"simulation_4d.dat")
    for i in infile:
        data = i.split()
        simulation.append(data)

simulation = np.array(simulation).astype(float)#.transpose()

print(simulation)
T = np.linspace(2,2.3,7)
NL = 7
NT = 7
simulation.reshape(NL, NT, -1)
for subarray in simulation:
    a,b,c,d,e,f,g = subarray
    plt.plot(T,a[0])
