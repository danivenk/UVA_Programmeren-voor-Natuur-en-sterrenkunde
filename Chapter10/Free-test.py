# Dani van Enk, 11823526
# 

# imports
from matplotlib import pyplot as ppl
from matplotlib import animation
import math
import cmath
import os
import random as rd
from matplotlib.backends.backend_pdf import PdfPages
from fractions import Fraction as F
from mpl_toolkits.mplot3d import Axes3D
import imageio as imio

# Path to current file
PATH = os.path.dirname(os.path.abspath(__file__))

# variables
hbar = 1#1.054571800e-34
m = 1#9.10938356e-31
pp = PdfPages(PATH + "/figures-free.pdf")


# lists
x_list = []
psi_list = []

class Potential:
    def __init__(self, a, V0):
        self.a = a
        self.V0 = V0

    def potential_value(self, x):
        if x < -self.a or x > self.a:
            return 0
        else:
            return self.V0

class Wave_Particle:
    def __init__(self, x, E, m, V0, a, deltax, Potential):
        self.x = x
        self.V = Potential
        self.E = E
        self.k1 = math.sqrt(2*m*self.E)
        self.k2 = cmath.sqrt(2*m*(self.E-self.V.V0))
        self.eta = self.k1/self.k2 + self.k2/self.k1
        self.x0 = x
        self.deltax = deltax

    def wave_function_b(self, x):
        A = 1
        G = 0
        F = cmath.exp(complex(0, 2*self.k1*self.V.a))/(cmath.cos(2*self.k2*self.V.a) - complex(0, self.eta/2*cmath.sin(2*self.k2*self.V.a)))
        C = F*(self.k1 + self.k2)/(2*self.k2)*cmath.exp(complex(0,(self.k1 - self.k2)*self.V.a))
        D = F*(self.k1 - self.k2)/(2*self.k2)*cmath.exp(complex(0,(self.k1 + self.k2)*self.V.a))
        B = F*complex(0,(self.k2**2-self.k1**2)/(2*self.k1*self.k2)*cmath.sin(2*self.k2*self.V.a))
        # print(self.x < -self.V.a, abs(self.x) <= self.V.a, self.x > -self.V.a, self.x)
        if x < -self.V.a:
            return A*cmath.exp(complex(0, self.k1*x)) + B*cmath.exp(complex(0, -self.k1*x))
        elif abs(x) <= self.V.a:
            return C*cmath.exp(complex(0, self.k2*x)) + D*cmath.exp(complex(0, -self.k2*x))
        elif x > -self.V.a:
            return F*cmath.exp(complex(0, self.k1*x)) + G*cmath.exp(complex(0, -self.k1*x))

    def a_of_k(self, k):
        N = 1
        return N*math.exp(-self.deltax**2*(self.k1-k)**2/8)*cmath.exp(complex(0, self.k1-k)*self.x0)

    # def monte_carlo(self, a, b, fmax, steps):
    #     area = (b-a)*fmax

    #     x_list_corr = []
    #     y_list_corr = []
    #     x_list_incorr = []
    #     y_list_incorr = []

    #     for step in range(0, steps):
    #         x_rand = (b-a)*rd.random()-a
    #         y_rand = fmax*rd.random()

    #         if y_rand < self.a_of_k(x_rand).real:
    #             x_list_corr.append(x_rand)
    #             y_list_corr.append(y_rand)
    #         else:
    #             x_list_incorr.append(x_rand)
    #             y_list_incorr.append(y_rand)

    #     return area*len(x_list_corr)/(len(x_list_corr) + len(x_list_incorr))

    def integrate(self, a, b, steps,x,t):
        int_value = 0
        k = a
        self.dk = 8/(self.deltax*steps)
        self.dk = F((b-a), steps)

        for step in range(0, steps+1):
            if 0 < step < steps+1:
                factor = 2
            else:
                factor = 1
            self.k1 += step*self.dk
            self.E = hbar**2*self.k1**2/2*m
            self.k2 = cmath.sqrt(self.k1**2-2*self.V.V0)
            # self.wave_function_b(x)*self.a_of_k(k)*cmath.exp(complex(0,self.E*t)/hbar)*self.dk
            int_value += self.dk/2 * \
                factor*self.a_of_k(k)#.wave_function_b(x)#*self.a_of_k(k) #* \
                # cmath.exp(complex(0, self.E*t)/hbar)
            k += self.dk

        return int_value

    def wave_function(self, x, t):
        # print(self.wave_function_b(x), self.monte_carlo(-1,1,10, 100000), cmath.exp(complex(0,self.E*t)/hbar))
        # return self.wave_function_b(x)*self.monte_carlo(-10,10,5, 10000)*cmath.exp(complex(0,self.E*t)/hbar)
        return self.integrate(-10*self.deltax,10*self.deltax, 1000, x, t)

V = Potential(.5,1)
WP = Wave_Particle(-10, .75, m, V.V0, V.a, 1, V)

# print(WP.wave_function_b(-.5))

for x_pos in range(-1000,1001):
    x_pos /= 100
    psi = WP.wave_function(x_pos,0)
    psi_star = psi.conjugate()
    print(x_pos)
    x_list.append(x_pos)
    psi_list.append((psi*psi_star).real)


ppl.figure(1)
ppl.plot(x_list, psi_list)
# ppl.xlim(-1,1)
# ppl.ylim(0,)
pp.savefig()

# plot show and pdf close
pp.close()
ppl.show()
