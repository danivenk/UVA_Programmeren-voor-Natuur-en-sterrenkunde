# Dani van Enk, 11823526
# hspec.py calculates the descrete energy eigenstates of the Lennard-Jones
#  potential via the action integral.
#  It also test the integral for some test values

# imports
from matplotlib import pyplot as ppl
import math
import os
from matplotlib.backends.backend_pdf import PdfPages

# Path to current file
PATH = os.path.dirname(os.path.abspath(__file__))

# variables
gamma = 21.7
steps = 1000
e_count = -0.8
pp = PdfPages(PATH + "/figures-hspec.pdf")

# lists
E_pot = []
x_list = []
e_eigen = []
action_list = []
e_list = []
s_e_list = []

# Lennard_Jones_Potential() defines the Lennard-Jones Potential
def Lennard_Jones_Potential(x):
    return 4*((1/x)**12-(1/x)**6)

# f() defines the integrant of the action integral
def f(e,v):
    return gamma*math.sqrt(e-v)

# action() finds the upper and lower bounds of the action integral
#  and calulates the integral via the bisection method for a specified e-value
def action(e):
    # variables
    s_e = 0
    step = 0

    # upper and lower bounds of the integral and dx
    x_in = ((-4+math.sqrt(4**2-4*e*-4))/(2*e))**(1/6)
    x_out = ((-4-math.sqrt(4**2-4*e*-4))/(2*e))**(1/6)
    dx = abs(x_out - x_in)/steps

    # start x
    x = x_in

    # calculating the action integral with the bisection method
    while step < steps-1:
        # this try/except rules out any possible imaginary outcome of a step
        try:
            if x == x_in or x == x_out:
                s_e += dx/2*(f(e,Lennard_Jones_Potential(x)))
            else:
                s_e += dx*(f(e,Lennard_Jones_Potential(x)))
        except ValueError:
            pass

        # next step
        x += dx
        step += 1

    return x_in,x_out,s_e,e

# eigenvalues() finds the eigenstate values of the potential starting at estart
#  and for the descrete value (n+1/2)pi
def eigenvalues(estart, value):
    # initial values
    e_l = estart
    e_r = -.01
    iteration = 0

    # using the bisection method to find the values of e
    #  for e-v=0 with v=(n+1/2)pi
    while True:
        # variable
        e_m = (e_l + e_r)/2

        # bisection method
        if action(e_m)[2] - value < 0:
            e_l = e_m
        elif action(e_m)[2] - value > 0:
            e_r = e_m

        # breaking if the difference is really small
        #  and prevending the loop getting stuck on 1 value
        if abs(action(e_m)[2] - value) < 1e-10:
            break
        elif iteration > 100:
            return iteration

        # next step
        iteration += 1

    return e_m

# printing all the boundry x positions for the test e-values
#  compairing them to the descrete values
for index in range(1,10):
    # variables
    e = -index/10
    x_in,x_out,s_e,e_value = action(e)

    # listing the test e-values and their corresponding action integral values
    e_list.append(e)
    action_list.append(s_e)

    # printing out the boundry x positions and their corresponding
    #  action integal values
    print("x_in = " + str(x_in) + " and x_out = " + str(x_out) + " with s(" + \
        str(e) + ") = " + str(s_e))

    # printing all analytical eigen values
    for index2 in range(0,10):
        if (index2 + 1/2)*math.pi < s_e:
            print("  s(" + str(e) + ") = " + str(s_e) + \
                " is larger than:\n    (n + 1/2)*\u03C0 = " + \
                    str((index2 + 1/2)*math.pi) + " with n = " + str(index2))

# making note of the posible values for the eigenstates
print("the eigen value of the energy for n = 0 lies between e = -0.8 and \
    e = -0.7")
print("the eigen value of the energy for n = 1 lies between e = -0.5 and \
    e = -0.4")
print("the eigen value of the energy for n = 2 lies between e = -0.2 and \
    e = -0.1")

# finding all eigen values of the energy between -0.8 and -0.01
for index in range(0, 20):
    # variables
    e_m = eigenvalues(e_count, (index + 1/2)*math.pi)

    # excluding any energies outside the interval [-1,0]
    if e_m < 0 and e_m > -1:
        # variables
        e_count = e_m

        # printing the found values and appending all the data from action(e)
        print("Found an eigen value at " + str(e_count) + " this is the " + \
            str(index+1) + "th eigenvalue below e = -0.01")
        s_e_list.append(action(e_m))
    else:
        break

# finding the analytical values for the Lennard-Jones Potential
for index in range(96, 596):
    index /= 100
    x_list.append(index)
    E_pot.append(Lennard_Jones_Potential(index))

# making a 3D matrix for the plots of the found descrete values
plot_list = [[[],[]] for _ in range(len(s_e_list))]

# setting up the Lennard-Jones potential plot
ppl.figure("Atomic Lennard-Jones Potential and eigenvalues")
ppl.title("Atomic Lennard-Jones Potential and eigenvalues")
ppl.plot(x_list, E_pot, color="Blue")

# plotting all found descrete energystates
for index in range(0, len(s_e_list)):
    # variables
    xpos = s_e_list[index][0]
    dx = abs(s_e_list[index][1]-s_e_list[index][0])/20

    # ploting the horizontal lines for the energy levels
    while xpos < s_e_list[index][1]:
        plot_list[index][0].append(xpos)
        plot_list[index][1].append(s_e_list[index][3])

        # next point
        xpos += dx

    # plotting all the energy levels
    ppl.plot(plot_list[index][0],plot_list[index][1], color="Red")

# setting up the rest of plot-1 and saving it
ppl.legend(["Lennard-Jones Potential", "Eigenvalue"])
ppl.xlim(0,6)
ppl.ylim(-1,1)
ppl.xlabel("position in bohr-radius")
ppl.ylabel("energy in V0")
pp.savefig()

# ploting the test values of the action integral and saving it
ppl.figure("Test values of the action integral (9 test values for e)")
ppl.title("Test values of the action integral")
ppl.plot(e_list,action_list)
ppl.xlim(-.9,-.1)
ppl.ylim(0,10)
ppl.xlabel("energy in V0")
ppl.ylabel("action in Js/V0")
pp.savefig()

# plot show and pdf close
pp.close()
ppl.show()
