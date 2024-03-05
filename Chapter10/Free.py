# Dani van Enk, 11823526
# Free.py simulates a quantum wave particle crossing a square potential of
#   V0 with energy E (> V0), it calculates |psi(x,t)|^2 for about 10s
# natrual constants are in effect
# It also calculates the reflection and transmission coefficients

# imports
from matplotlib import pyplot as ppl
import math
import cmath
import os
from matplotlib.backends.backend_pdf import PdfPages

# Path to current file
PATH = os.path.dirname(os.path.abspath(__file__))

# variables
hbar = 1#1.054571800e-34
m = 1#9.10938356e-31
tmax_steps = 100
pp = PdfPages(PATH + "/figures-free.pdf")

# lists
T_list = []
R_list = []
TR_list = []
time_list = []

# the Potential class creates a square potential given the width a and
#   the height V0
# it has a function to get the potential value
class Potential:
    # __int__() takes the width and height of the potential
    def __init__(self, a, V0):
        self.a = a
        self.V0 = V0

    # potential_value() gives the value of the potential given the  position
    def potential_value(self, x):
        if x < -self.a or x > self.a:
            return 0
        else:
            return self.V0

# the Wave_Particle class creates a quantum wave particle given the x position,
#   the energy, the mass, the potential (height, width and object itself)
# it also has a function wich gives the boundary wavefunction
# it also has a function wich gives the fourier transform of psi(x,0)
#   as function of k
# it also has a function to calculate the time dependent wavefunction
class Wave_Particle:
    # __init__() contains the begin parameters of the Wave_Particle class
    #   given the x position, energy (< V0), the potential (height,
    #   width and object itself)
    # it also calculates some initial variables
    def __init__(self, x, E, m, V0, a, Potential):
        self.x = x
        self.V = Potential
        self.E = E
        self.k10 = math.sqrt(2*m*self.E)
        self.k1 = self.k10
        self.k2 = cmath.sqrt(2*m*(self.E-self.V.V0))
        self.x0 = x
        self.deltax = 1
        self.eta = self.k1/self.k2 + self.k2/self.k1

    # wave_function_b() gives the boundary wavefunction
    def wave_function_b(self, x):
        A = 1
        G = 0
        F = cmath.exp(complex(0, -2*self.k1*self.V.a))/ \
            (cmath.cos(2*self.k2*self.V.a) - \
            complex(0, self.eta/2*cmath.sin(2*self.k2*self.V.a)))
        C = F*(self.k1 + self.k2)/(2*self.k2)* \
            cmath.exp(complex(0,(self.k1 - self.k2)*self.V.a))
        D = F*(self.k1 - self.k2)/(2*self.k2)* \
            cmath.exp(complex(0,(self.k1 + self.k2)*self.V.a))
        B = F*complex(0, (self.k2**2-self.k1**2) / \
            (2*self.k1*self.k2)*cmath.sin(2*self.k2*self.V.a))

        # return the correct parts of the boundary wavefunction
        if x < -self.V.a:
            return A*cmath.exp(complex(0, self.k1*x)) + \
                B*cmath.exp(complex(0, -self.k1*x))
        elif abs(x) <= self.V.a:
            return C*cmath.exp(complex(0, self.k2*x)) + \
                D*cmath.exp(complex(0, -self.k2*x))
        elif x > -self.V.a:
            return F*cmath.exp(complex(0, self.k1*x)) + \
                G*cmath.exp(complex(0, -self.k1*x))

    # phi_k() returns the analytical fourier transform of psi(x,0)
    #    with the analytical normalization factor
    def phi_k(self):
        N = (self.deltax/(2*math.pi**3))**(1/4)
        return N*math.exp(-self.deltax**2*(self.k1-self.k10)**2/8)* \
            cmath.exp(-complex(0, self.k1-self.k10)*self.x0)

    # integrate_psi() returns psi for (x,t) using steps in k (Riemann sum)
    def integrate_psi(self, steps,x,t):
        # variables
        int_value_psi = 0
        dk1 = 8/(self.deltax*steps)

        # for each step in k find the area, change k1, k2 and E
        for step in range(0, steps):
            int_value_psi += dk1*self.phi_k()*self.wave_function_b(x) * \
                cmath.exp(-complex(0, self.E*t)/hbar)
            self.k1 = self.k10 + step*dk1
            self.E = hbar**2*self.k1**2/(2*m)
            self.k2 = cmath.sqrt(self.k1**2-2*self.V.V0)

        return int_value_psi

    # wave_function() returns psi for a given (x,t)
    def wave_function(self, x, t):
        return self.integrate_psi(100, x, t)

# creating a potential object and a wave particle object 
V = Potential(.5,1)
WP = Wave_Particle(-10, 1.5, m, V.V0, V.a, V)

# plot setup
fig = ppl.figure("Wave function of a qyantum particle")

# get the |psi|^2 for different times
for time_step in range(0,tmax_steps+1):
    # variables
    time = time_step / 10
    area = 0
    area_R = 0
    area_T = 0

    # lists
    x_list = []
    psi_list = []

    # step plot setup
    ppl.clf()
    ppl.title("|\u03C8|^2 of a quantum particle at t = " + str(time))
    ppl.xlim(-10, 10)
    ppl.ylim(0, 1)
    ppl.xlabel("x position")
    ppl.ylabel("|\u03C8|^2")

    # loop over the x-axis
    for x_pos in range(-200,201):
        # change x to go in steps of .1
        x_pos /= 10

        # find psi for given (x, t)
        psi = WP.wave_function(x_pos,time)
        psi_star = psi.conjugate()

        # append x and |psi|^2
        x_list.append(x_pos)
        psi_list.append((psi*psi_star).real)

        # get area
        area += (psi*psi_star).real

        # find the reflected and transmitted area
        if x_pos < -V.a:
            area_R += (psi*psi_star).real
        elif x_pos > V.a:
            area_T += (psi*psi_star).real

    # calculate transmission and reflection coefficients
    T = area_T/area
    R = area_R/area

    # print transmission and reflection coefficients after time simulation
    if time == tmax_steps/10:
        print("Transmission", T)
        print("Reflection", R)
        print("Total", T+R)

    # plot |psi|^2 as function of x and show position of the potential
    ppl.plot(x_list, psi_list)
    ppl.plot([-V.a,-V.a,-V.a], [0,.5,1], color="black", linewidth=.5)
    ppl.plot([V.a,V.a,V.a], [0,.5,1], color="black", linewidth=.5)
    ppl.pause(0.01)

    # every five steps save the plot
    if time_step % 5 == 0:
        pp.savefig()

# plot show and pdf close
pp.close()
ppl.show()
