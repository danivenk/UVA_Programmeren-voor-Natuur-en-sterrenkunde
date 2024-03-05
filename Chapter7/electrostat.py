# Dani van Enk, 11823526
# electrostat.py calculates the potential according to Lap(V) = 0 numerically
#   using the Gauss-Siedel method and compares it with the analytical solution
# the problems calculated are: a plate allong the x-axis, a capacitor in the
#   y-direction, 2 poincharges at y = 0.5 and x = 0 & x = 1, the own
#   configuration is a capacitor with V0 on one end and 2V0 on the other, all
#   all analytical solutions are the best I could find

# imports
from matplotlib import pyplot as ppl
import math
import os
from matplotlib.backends.backend_pdf import PdfPages

# Path to current file
PATH = os.path.dirname(os.path.abspath(__file__))

# variables
a = 1
b = 1
V0 = 10
steps = 101
V_no = 1
V_capacitor_no = 1
V_pointcharge_no = 1
V_ownconf_no = 1
q = 50
dx = a/(steps-1)
dy = b/(steps-1)
pp = PdfPages(PATH + "/figures-electrostat.pdf")
ppextra = PdfPages(PATH + "/figures-electrostat-extra.pdf")

# lists
V_no_list = [0, 1]
V_capacitor_no_list = [0, 1]
V_pointcharge_no_list = [0, 1]
V_ownconf_no_list = [0, 1]
E_list = []
E_capacitor_list = []
E_pointcharge_list = []
E_ownconf_list = []
x_list = []
y_list = []
V_x_list = []
Va_x_list = []
V_capacitor_x_list = []
Va_capacitor_x_list = []
V_pointcharge_x_list = []
Va_pointcharge_x_list = []
V_ownconf_x_list = []
Va_ownconf_x_list = []
V_y_list = []
Va_y_list = []
V_capacitor_y_list = []
Va_capacitor_y_list = []
V_pointcharge_y_list = []
Va_pointcharge_y_list = []
V_ownconf_y_list = []
Va_ownconf_y_list = []

# this is a constant class from vana.py
# it creates a class with steps, top & right boundry, V0 and maxsteps
class constants:
    def __init__(self, nstep, a, b, V0, nmax):
        self.nstep = nstep
        self.a = a
        self.b = b
        self.V0 = V0
        self.nmax = nmax

# creates an object using the class constants from vana.py
myconstants = constants(100, 1., 1., 10., 101)

# F1() is part of the analytical solution according to vana.py
def F1(n):
    return (1./n)*(1-math.pow(-1., n))

# F2() is part of the analytical solution according to vana.py
def F2(n, x):
    return math.sin((math.pi*n*x)/myconstants.a)

# F3() is part of the analytical solution according to vana.py
def F3(n, y):
    return ((math.exp(((math.pi*n)/myconstants.a)*(myconstants.b-y)) - \
             math.exp(((-math.pi*n)/myconstants.a)*(myconstants.b-y)))/2.)

# F4() is part of the analytical solution according to vana.py
def F4(n):
    return (math.exp((math.pi*n*myconstants.b)/myconstants.a) - \
            math.exp((-math.pi*n*myconstants.b)/myconstants.a))/2.

# Va() is the analytical solution of the problem according to vana.py
def Va(x, y):
    V = 0
    Vcontribution = 0
    for n in range(1, myconstants.nmax, 2):
        Vcontribution = (2.*myconstants.V0/math.pi) * \
            F1(n)*F2(n, x)*F3(n, y)/F4(n)
        V += Vcontribution

        # I added this try except to prevent a Vcontribution/V = infinity
        try:
            if math.fabs(Vcontribution/V) < math.pow(10., -20):
                break
        except ZeroDivisionError:
            break
    return V

# Va_capacitor() is the analytical solution for the capacitor problem
def Va_capacitor(x, y):
    C = 2*V0
    V = C*y - V0
    return V

# Va_pointcharge() is the analytical solution for the pointcharge problem
def Va_pointcharge(x, y):
    k = 1.602e-19/(4*math.pi*8.85e-12)
    try:
        return k*(1/(math.sqrt(x**2 + (y - q/(steps-1))**2))**2 + \
            1/(math.sqrt((x - a)**2 + (y - q/(steps-1))**2))**2)
    except ZeroDivisionError:
        return V0

# Va_ownconf() is the analytical solution for the own configuration problem
def Va_ownconf(x, y):
    C = V0
    V = C*y + V0
    return V

# Gauss_Seidel() uses the Gauss-Seidel method to calculate the Potential values
def Gauss_Seidel(V,i,j):
    V[i][j] = 1/4*(V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1])
    return V[i][j]

# E() calculates the smoothness of a potential matrix
def E(V):
    E = 0
    for i in range(1, steps):
        for j in range(1, steps):
            E += (V[i][j]-V[i-1][j])**2 + (V[i][j]-V[i][j-1])**2
    return E
            

# creating a steps x steps matrix
V_matrix = [[] for _ in range(0, steps)]
Va_matrix = [[] for _ in range(0, steps)]
V_capacitor_matrix = [[] for _ in range(0, steps)]
Va_capacitor_matrix = [[] for _ in range(0, steps)]
V_pointcharge_matrix = [[] for _ in range(0, steps)]
Va_pointcharge_matrix = [[] for _ in range(0, steps)]
V_ownconf_matrix = [[] for _ in range(0, steps)]
Va_ownconf_matrix = [[] for _ in range(0, steps)]

# filling the matrix with zeros
for i in range(0,steps):
    j = 0
    while j < steps:
        V_matrix[i].append(0.)
        Va_matrix[i].append(0.)
        V_capacitor_matrix[i].append(0.)
        Va_capacitor_matrix[i].append(0.)
        V_pointcharge_matrix[i].append(0.)
        Va_pointcharge_matrix[i].append(0.)
        V_ownconf_matrix[i].append(0.)
        Va_ownconf_matrix[i].append(0.)
        j += 1

# setting boundary conditions along the x = 0 and x = a
for i in range(0, steps):
    V_matrix[i][0] = 0
    Va_matrix[i][0] = 0
    V_capacitor_matrix[i][0] = 0
    Va_capacitor_matrix[i][0] = 0
    V_pointcharge_matrix[i][0] = 0
    Va_pointcharge_matrix[i][0] = 0
    V_ownconf_matrix[i][0] = 0
    Va_ownconf_matrix[i][0] = 0
    V_matrix[i][steps - 1] = 0
    Va_matrix[i][steps - 1] = 0
    V_capacitor_matrix[i][steps - 1] = 0
    Va_capacitor_matrix[i][steps - 1] = 0
    V_pointcharge_matrix[i][steps - 1] = 0
    Va_pointcharge_matrix[i][steps - 1] = 0
    V_ownconf_matrix[i][steps - 1] = 0
    Va_ownconf_matrix[i][steps - 1] = 0

# setting boundary conditions along the y = 0 and y = b
for j in range(0, steps):
    V_matrix[0][j] = V0
    Va_matrix[0][j] = V0
    V_capacitor_matrix[0][j] = -V0
    Va_capacitor_matrix[0][j] = -V0
    V_pointcharge_matrix[0][j] = 0
    Va_pointcharge_matrix[0][j] = 0
    V_ownconf_matrix[0][j] = V0
    Va_ownconf_matrix[0][j] = V0
    V_matrix[steps - 1][j] = 0
    Va_matrix[steps - 1][j] = 0
    V_capacitor_matrix[steps - 1][j] = V0
    Va_capacitor_matrix[steps - 1][j] = V0
    V_pointcharge_matrix[steps - 1][j] = 0
    Va_pointcharge_matrix[steps - 1][j] = 0
    V_ownconf_matrix[steps - 1][j] = 2*V0
    Va_ownconf_matrix[steps - 1][j] = 2*V0

# setting a point charges at (x,y) -> (.5,0) & (.5,1)
V_pointcharge_matrix[q][0] = V0
Va_pointcharge_matrix[q][0] = V0
V_pointcharge_matrix[q][steps - 1] = V0
Va_pointcharge_matrix[q][steps - 1] = V0

# Calculating the smoothness factor for the first iteration
E_list.append(E(V_matrix))
E_capacitor_list.append(E(V_capacitor_matrix))
E_pointcharge_list.append(E(V_pointcharge_matrix))
E_ownconf_list.append(E(V_ownconf_matrix))

# Calculating the potential using the Gaus-Seidel method
for i in range(1, steps-1):
    for j in range(1, steps-1):
        V_matrix[i][j] = Gauss_Seidel(V_matrix, i, j)
        V_capacitor_matrix[i][j] = Gauss_Seidel(V_capacitor_matrix, i, j)
        V_pointcharge_matrix[i][j] = Gauss_Seidel(V_pointcharge_matrix, i, j)
        V_ownconf_matrix[i][j] = Gauss_Seidel(V_ownconf_matrix, i, j)

# Calculating the smoothness factor for the second iteration
E_list.append(E(V_matrix))
E_capacitor_list.append(E(V_capacitor_matrix))
E_pointcharge_list.append(E(V_pointcharge_matrix))
E_ownconf_list.append(E(V_ownconf_matrix))

# Iterating further untill the difference for the smoothness factor
# is smaller than a 1/10000 of a percent for the first problem
while abs(E_list[V_no] - E_list[V_no-1])/E_list[V_no] > 1e-6:
    # Calculating the potential using the Gaus-Seidel method
    for i in range(1, steps-1):
        for j in range(1, steps-1):
            V_matrix[i][j] = Gauss_Seidel(V_matrix, i, j)

    # Calculating the smoothness factor for the each iteration step
    # and counting the number of iterating steps   
    E_list.append(E(V_matrix))
    V_no_list.append(V_no)
    V_no += 1

    print("first", abs(E_list[V_no] - E_list[V_no-1])/E_list[V_no])

# Iterating further untill the difference for the smoothness factor
# is smaller than a 1/10000 of a percent for the capacitor problem
while abs(E_capacitor_list[V_capacitor_no] - \
    E_capacitor_list[V_capacitor_no-1])/E_capacitor_list[V_capacitor_no] > 1e-6:
    # Calculating the potential using the Gaus-Seidel method
    for i in range(1, steps-1):
        for j in range(1, steps-1):
            V_capacitor_matrix[i][j] = \
                Gauss_Seidel(V_capacitor_matrix, i, j)

    # Calculating the smoothness factor for the each iteration step
    # and counting the number of iterating steps
    E_capacitor_list.append(E(V_capacitor_matrix))
    V_capacitor_no_list.append(V_capacitor_no)
    V_capacitor_no += 1

    print("capacitor", abs(E_capacitor_list[V_capacitor_no] - \
        E_list[V_capacitor_no-1])/E_capacitor_list[V_capacitor_no])
    
# Iterating further untill the difference for the smoothness factor
# is smaller than a 1/10000 of a percent for the pointcharge problem
while abs(E_pointcharge_list[V_pointcharge_no] - E_pointcharge_list \
    [V_pointcharge_no-1])/E_pointcharge_list[V_pointcharge_no] > 1e-6:
    # Calculating the potential using the Gaus-Seidel method
    for i in range(1, steps-1):
        for j in range(1, steps-1):
            V_pointcharge_matrix[i][j] = \
                Gauss_Seidel(V_pointcharge_matrix, i, j)

    # Calculating the smoothness factor for the each iteration step
    # and counting the number of iterating steps
    E_pointcharge_list.append(E(V_pointcharge_matrix))
    V_pointcharge_no_list.append(V_pointcharge_no)
    V_pointcharge_no += 1

    print("pointcharge", abs(E_pointcharge_list[V_pointcharge_no] - \
        E_pointcharge_list[V_pointcharge_no-1])/E_pointcharge_list \
        [V_pointcharge_no])
    
# Iterating further untill the difference for the smoothness factor
# is smaller than a 1/10000 of a percent for the own configuration
while abs(E_ownconf_list[V_ownconf_no] - E_ownconf_list \
    [V_ownconf_no-1])/E_ownconf_list[V_ownconf_no] > 1e-6:
    # Calculating the potential using the Gaus-Seidel method
    for i in range(1, steps-1):
        for j in range(1, steps-1):
            V_ownconf_matrix[i][j] = \
                Gauss_Seidel(V_ownconf_matrix, i, j)

    # Calculating the smoothness factor for the each iteration step
    # and counting the number of iterating steps
    E_ownconf_list.append(E(V_ownconf_matrix))
    V_ownconf_no_list.append(V_ownconf_no)
    V_ownconf_no += 1

    print("ownconf", abs(E_ownconf_list[V_ownconf_no] - \
        E_ownconf_list[V_ownconf_no-1])/E_ownconf_list[V_ownconf_no])

# Calculating the potential using the analytical solution
for x_value in range(1, steps):
    for y_value in range(1, steps):
        Va_matrix[y_value][x_value] = Va(x_value*dx, y_value*dy)
        Va_capacitor_matrix[y_value][x_value] = \
            Va_capacitor(x_value*dx, y_value*dy)
        Va_pointcharge_matrix[y_value][x_value] = \
            Va_pointcharge(x_value*dx, y_value*dy)
        Va_ownconf_matrix[y_value][x_value] = \
            Va_ownconf(x_value*dx, y_value*dy)

# Listing the potential values at y = 0.5
for x_value in range(0, steps):
    x_list.append(x_value*dx)
    V_x_list.append(V_matrix[int(len(V_matrix)/2)][x_value])
    Va_x_list.append(Va_matrix[int(len(V_matrix)/2)][x_value])
    V_capacitor_x_list.append(V_capacitor_matrix \
        [int(len(V_capacitor_matrix)/2)][x_value])
    Va_capacitor_x_list.append(Va_capacitor_matrix \
        [int(len(V_capacitor_matrix)/2)][x_value])
    V_pointcharge_x_list.append(V_pointcharge_matrix \
        [int(len(V_pointcharge_matrix)/2)][x_value])
    Va_pointcharge_x_list.append(Va_pointcharge_matrix \
        [int(len(V_pointcharge_matrix)/2)][x_value])
    V_ownconf_x_list.append(V_ownconf_matrix[int(len(V_ownconf_matrix)/2)] \
        [x_value])
    Va_ownconf_x_list.append(Va_ownconf_matrix[int(len(V_ownconf_matrix)/2)] \
        [x_value])

# Listing the potential values at x = 0.5 (point charge x = 0)
for y_value in range(0, steps):
    y_list.append(y_value*dy)
    V_y_list.append(V_matrix[y_value][int(len(V_matrix)/2)])
    Va_y_list.append(Va_matrix[y_value][int(len(V_matrix)/2)])
    V_capacitor_y_list.append(V_capacitor_matrix[y_value] \
        [int(len(V_capacitor_matrix)/2)])
    Va_capacitor_y_list.append(Va_capacitor_matrix[y_value] \
        [int(len(V_capacitor_matrix)/2)])
    V_pointcharge_y_list.append(V_pointcharge_matrix[y_value][0])
    Va_pointcharge_y_list.append(Va_pointcharge_matrix[y_value][0])
    V_ownconf_y_list.append(V_ownconf_matrix[y_value] \
        [int(len(V_ownconf_matrix)/2)])
    Va_ownconf_y_list.append(Va_ownconf_matrix[y_value] \
        [int(len(V_ownconf_matrix)/2)])

# setting up contour plots 1&2 and saving them to figures-electrostat.pdf
ppl.figure("A contour plot of the numerical solution to the first problem")
ppl.title("A contour plot of the numerical solution to the first problem")
ppl.contourf(x_list, y_list, V_matrix)
ppl.xlim(0, a)
ppl.ylim(0, b)
ppl.xlabel("x position")
ppl.ylabel("y position")
pp.savefig()
ppl.figure("A contour plot of the analytical solution to the first problem")
ppl.title("A contour plot of the analytical solution to the first problem")
ppl.contourf(x_list, y_list, Va_matrix)
ppl.xlim(0, myconstants.a)
ppl.ylim(0, myconstants.b)
ppl.xlabel("x position")
ppl.ylabel("y position")
pp.savefig()

# setting up figure 3 and saving them to figures-electrostat.pdf
ppl.figure("Convergence monitor for the first problem")
ppl.title("Convergence monitor for the first problem")
ppl.plot(V_no_list,E_list)
ppl.xlim(0, V_no_list[len(V_no_list)-1])
ppl.ylim(0, E_list[0])
ppl.xlabel("Number iterations")
ppl.ylabel("Smoothness after iterating")
pp.savefig()

# setting up figure 4 and saving them to figures-electrostat.pdf
ppl.figure("Potential along the x-asis at y = 0.5 for the first problem")
ppl.title("Potential along the x-asis at y = 0.5 for the first problem")
ppl.plot(x_list, V_x_list)
ppl.plot(x_list, Va_x_list)
ppl.legend(["Numerical Potential", "Analytical Potential"])
ppl.xlim(0,a)
ppl.ylim(0,V0)
ppl.xlabel("x position")
ppl.ylabel("Potential in J")
pp.savefig()

# setting up figure 5 and saving them to figures-electrostat.pdf
ppl.figure("Potential along the y-asis at x = 0.5 for the first problem")
ppl.title("Potential along the y-asis at x = 0.5 for the first problem")
ppl.plot(y_list, V_y_list)
ppl.plot(y_list, Va_y_list)
ppl.legend(["Numerical Potential", "Analytical Potential"])
ppl.xlim(0,a)
ppl.ylim(0,V0)
ppl.xlabel("y position")
ppl.ylabel("Potential in J")
pp.savefig()

# setting up contour plots 1&2 and saving them to figures-electrostat-extra.pdf
ppl.figure("A contour plot of the numerical solution to the capacitor problem")
ppl.title("A contour plot of the numerical solution to the capacitor problem")
ppl.contourf(x_list, y_list, V_capacitor_matrix)
ppl.xlim(0, a)
ppl.ylim(0, b)
ppl.xlabel("x position")
ppl.ylabel("y position")
ppextra.savefig()
ppl.figure("A contour plot of the analytical solution to the capacitor problem")
ppl.title("A contour plot of the analytical solution to the capacitor problem")
ppl.contourf(x_list, y_list, Va_capacitor_matrix)
ppl.xlim(0, a)
ppl.ylim(0, b)
ppl.xlabel("x position")
ppl.ylabel("y position")
ppextra.savefig()

# setting up figure 3 and saving them to figures-electrostat-extra.pdf
ppl.figure("Convergence monitor for the capacitor problem")
ppl.title("Convergence monitor for the capacitor problem")
ppl.plot(V_capacitor_no_list, E_capacitor_list)
ppl.xlim(0, V_capacitor_no_list[len(V_capacitor_no_list)-1])
ppl.ylim(0, E_capacitor_list[0])
ppl.xlabel("Number iterations")
ppl.ylabel("Smoothness after iterating")
ppextra.savefig()

# setting up figure 4 and saving them to figures-electrostat-extra.pdf
ppl.figure("Potential along the x-asis at y = 0.5 for the capacitor problem")
ppl.title("Potential along the x-asis at y = 0.5 for the capacitor problem")
ppl.plot(x_list, V_capacitor_x_list)
ppl.plot(x_list, Va_capacitor_x_list)
ppl.legend(["Numerical Potential", "Analytical Potential"])
ppl.xlim(0, a)
ppl.ylim(-V0, V0)
ppl.xlabel("x position")
ppl.ylabel("Potential in J")
ppextra.savefig()

# setting up figure 5 and saving them to figures-electrostat-extra.pdf
ppl.figure("Potential along the y-asis at x = 0.5 for the capacitor problem")
ppl.title("Potential along the y-asis at x = 0.5 for the capacitor problem")
ppl.plot(y_list, V_capacitor_y_list)
ppl.plot(y_list, Va_capacitor_y_list)
ppl.legend(["Numerical Potential", "Analytical Potential"])
ppl.xlim(0, a)
ppl.ylim(-V0, V0)
ppl.xlabel("y position")
ppl.ylabel("Potential in J")
ppextra.savefig()

# setting up contour plots 6&7 and saving them to figures-electrostat-extra.pdf
ppl.figure("A contour plot of the numerical solution to the pointcharge problem")
ppl.title("A contour plot of the numerical solution to the pointcharge problem")
ppl.contourf(x_list, y_list, V_pointcharge_matrix)
ppl.xlim(0, a)
ppl.ylim(0, b)
ppl.xlabel("x position")
ppl.ylabel("y position")
ppextra.savefig()
ppl.figure("A contour plot of the analytical solution to the pointcharge problem")
ppl.title("A contour plot of the analytical solution to the pointcharge problem")
ppl.contourf(x_list, y_list, Va_pointcharge_matrix)
ppl.xlim(0, a)
ppl.ylim(0, b)
ppl.xlabel("x position")
ppl.ylabel("y position")
ppextra.savefig()

# setting up figure 8 and saving them to figures-electrostat-extra.pdf
ppl.figure("Convergence monitor for the pointcharge problem")
ppl.title("Convergence monitor for the pointcharge problem")
ppl.plot(V_pointcharge_no_list, E_pointcharge_list)
ppl.xlim(0, V_pointcharge_no_list[len(V_pointcharge_no_list)-1])
ppl.ylim(0, E_pointcharge_list[0])
ppl.xlabel("Number iterations")
ppl.ylabel("Smoothness after iterating")
ppextra.savefig()

# setting up figure 9 and saving them to figures-electrostat-extra.pdf
ppl.figure("Potential along the x-asis at y = 0.5 for the pointcharge problem")
ppl.title("Potential along the x-asis at y = 0.5 for the pointcharge problem")
ppl.plot(x_list, V_pointcharge_x_list)
ppl.plot(x_list, Va_pointcharge_x_list)
ppl.legend(["Numerical Potential", "Analytical Potential"])
ppl.xlim(0, a)
ppl.ylim(0, V0)
ppl.xlabel("x position")
ppl.ylabel("Potential in J")
ppextra.savefig()

# setting up figure 10 and saving them to figures-electrostat-extra.pdf
ppl.figure("Potential along the y-asis at x = 0 for the pointcharge problem")
ppl.title("Potential along the y-asis at x = 0 for the pointcharge problem")
ppl.plot(y_list, V_pointcharge_y_list)
ppl.plot(y_list, Va_pointcharge_y_list)
ppl.legend(["Numerical Potential", "Analytical Potential"])
ppl.xlim(0, a)
ppl.ylim(0, V0)
ppl.xlabel("y position")
ppl.ylabel("Potential in J")
ppextra.savefig()

# setting up contour plots 11&12 and saving them to figures-electrostat-extra.pdf
ppl.figure("A contour plot of the numerical solution to the own configuration")
ppl.title("A contour plot of the numerical solution to the own configuration")
ppl.contourf(x_list, y_list, V_ownconf_matrix)
ppl.xlim(0, a)
ppl.ylim(0, b)
ppl.xlabel("x position")
ppl.ylabel("y position")
ppextra.savefig()
ppl.figure("A contour plot of the analytical solution to the own configuration")
ppl.title("A contour plot of the analytical solution to the own configuration")
ppl.contourf(x_list, y_list, Va_ownconf_matrix)
ppl.xlim(0, a)
ppl.ylim(0, b)
ppl.xlabel("x position")
ppl.ylabel("y position")
ppextra.savefig()

# setting up figure 13 and saving them to figures-electrostat-extra.pdf
ppl.figure("Convergence monitor for the own configuration")
ppl.title("Convergence monitor for the own configuration")
ppl.plot(V_ownconf_no_list, E_ownconf_list)
ppl.xlim(0, V_ownconf_no_list[len(V_ownconf_no_list)-1])
ppl.ylim(0, E_ownconf_list[0])
ppl.xlabel("Number iterations")
ppl.ylabel("Smoothness after iterating")
ppextra.savefig()

# setting up figure 14 and saving them to figures-electrostat-extra.pdf
ppl.figure("Potential along the x-asis at y = 0.5 for the own configuration")
ppl.title("Potential along the x-asis at y = 0.5 for the own configuration")
ppl.plot(x_list, V_ownconf_x_list)
ppl.plot(x_list, Va_ownconf_x_list)
ppl.legend(["Numerical Potential", "Analytical Potential"])
ppl.xlim(0, a)
ppl.ylim(0, 1.5*V0)
ppl.xlabel("x position")
ppl.ylabel("Potential in J")
ppextra.savefig()

# setting up figure 15 and saving them to figures-electrostat-extra.pdf
ppl.figure("Potential along the y-asis at x = 0.5 for the own configuration")
ppl.title("Potential along the y-asis at x = 0.5 for the own configuration")
ppl.plot(y_list, V_ownconf_y_list)
ppl.plot(y_list, Va_ownconf_y_list)
ppl.legend(["Numerical Potential", "Analytical Potential"])
ppl.xlim(0, a)
ppl.ylim(0, 2*V0)
ppl.xlabel("y position")
ppl.ylabel("Potential in J")
ppextra.savefig()

# plot show and pdf close
pp.close()
ppextra.close()
ppl.show()
