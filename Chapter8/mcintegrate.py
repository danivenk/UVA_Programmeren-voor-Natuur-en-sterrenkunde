# Dani van Enk, 11823526
# mcintegrate.py calculates pi using the monte-carlo integration method
#   it also calculates the integral of a gaussian function with
#   a standard deviation of 1 and top at x = 0 using the monte-carlo method
# Lastly it calculates the the integration value of the gaussian as a function
#   of N and the error in that value

# imports
from matplotlib import pyplot as ppl
import math
import os
import random as rd
from matplotlib.backends.backend_pdf import PdfPages

# Path to current file
PATH = os.path.dirname(os.path.abspath(__file__))

# variables
steps = 100000 
std_steps = 100
N = 1
xmin = -1
xmax = 1
ymin = -1
ymax = 1
area = abs(xmax-xmin)*abs(ymax-ymin)
area_g = abs(xmax-xmin)*abs(ymax-0)
pp = PdfPages(PATH + "/figures-mcintegrate.pdf")

# lists
N_list = []
N_sq_inverse_list = []
int_gaussian = []
int_gaussian_error = []

# stepping_up() is a logarithmic step function which causes the steps
#   to get bigger with the power of 10
def stepping_up(N):
    for index in range(0, 10):
        if N < 10**(index+1) and N >= 10**(index):
            return 10**index

# Gaussian() is a gaussian function with standard deviation of 1
#   a center at x = 0
def Gaussian(x, std, mu):
    return math.exp(-((x-mu)**2/(2*std**2)))

# random_number() generates a number between a and b
def random_number(a, b):
    rd.seed()
    rd_nb = rd.random()
    diff = abs(b-a)
    return_value = rd_nb*diff + a
    return return_value

# monte_carlo() checks if a random point within the specified interval falls 
#   within the specified circle of radius 1 using the monte-carlo method
# it does this steps times and returns lists of the inside and outside points
def monte_carlo_circle(xmin, xmax, ymin, ymax,steps):
    # initial step value
    step = 0

    # lists
    x_list_inside = []
    x_list_outside = []
    y_list_inside = []
    y_list_outside = []

    # loop for each step
    while step < steps:
        # get random x, y and r values for a circle or radius 1
        rd_x = random_number(xmin, xmax)
        rd_y = random_number(ymin, ymax)
        rd_r = rd_x**2 + rd_y**2

        # check if the points are inside or outside the circle of radius 1
        if rd_r > 1:
            x_list_outside.append(rd_x)
            y_list_outside.append(rd_y)
        if rd_r < 1:
            x_list_inside.append(rd_x)
            y_list_inside.append(rd_y)
        
        # next step
        step += 1

    return x_list_inside, x_list_outside, y_list_inside, y_list_outside

# monte_carlo() checks if a random point within the specified interval falls
#   within the specified circle of radius 1 using the monte-carlo method
# it does this steps times and returns lists of the inside and outside points
def monte_carlo_gaussian(xmin, xmax, ymin, ymax, steps):
    # initial step value
    step = 0

    # lists
    x_list_inside = []
    x_list_outside = []
    y_list_inside = []
    y_list_outside = []

    # loop for each step
    while step < steps:
        # generate random x, y and f values for gaussian
        rd_x = random_number(xmin, xmax)
        rd_y = random_number(ymin, ymax)
        rd_f = Gaussian(rd_x, 1, 0)

        # check if the points are above (outside) or under (inside) the function
        if rd_f < rd_y:
            x_list_outside.append(rd_x)
            y_list_outside.append(rd_y)
        if rd_f > rd_y:
            x_list_inside.append(rd_x)
            y_list_inside.append(rd_y)

        # next step
        step += 1

    return x_list_inside, x_list_outside, y_list_inside, y_list_outside

# get the inside/outside lists for x and y for the circle using steps steps
x_list_ic, x_list_oc, y_list_ic, y_list_oc = \
    monte_carlo_circle(xmin, xmax, ymin, ymax, steps)

# Loop N for steps times (while here because the increment of N is variable)
while N <= steps:
    # initial values
    total = 0
    total_sq = 0

    # generate std_steps, N point integration values 
    for N_std in range(0,std_steps):
        x_list_ig_N, x_list_og_N, y_list_ig_N, y_list_og_N = \
            monte_carlo_gaussian(xmin, xmax, 0, ymax, N)  

        # calculate the integral 
        integralN = area_g*(len(x_list_ig_N)/(len(x_list_ig_N) + \
            len(x_list_og_N)))

        # add the calculated value to the total and squared total
        total += integralN
        total_sq += integralN**2

        # print current loop
        print("N", N, "STD_N", N_std)

    # calculate average and standard deviation
    average = total/std_steps
    average_sq = total_sq/std_steps
    std_dv = math.sqrt(average_sq - average**2)

    # append N, 1/sqrt(N),the integrated value and errors
    N_list.append(N)
    N_sq_inverse_list.append(1/math.sqrt(N))
    int_gaussian.append(average)
    int_gaussian_error.append(std_dv)

    # add to N according to stepping_up()
    N += stepping_up(N)

# calculate pi
pi = area*len(x_list_ic)/(len(x_list_ic) + len(x_list_oc))

# printing the values and relative difference to the real value for pi and
#   the intral
print("\u03C0 = " + str(pi))
print("The calculated value of \u03C0 is " + str(abs((pi-math.pi)/math.pi)*100) \
    + " percent from the real value")
print("The integral of a gausian function with standard deviation of 1 and " + \
    " center at 0 from -1 to 1 is " + str(int_gaussian[len(int_gaussian)-1]))
print("The calculated value of the integeral is " + \
    str(abs((int_gaussian[len(int_gaussian)-1]-math.sqrt(2*math.pi)* \
    math.erf(1/math.sqrt(2)))/math.sqrt(2*math.pi)* \
    math.erf(1/math.sqrt(2)))*100) + " percent from the real value")

# setup plot 1
ppl.figure("Random points inside or outside a circle " + \
    "using the monte-carlo method")
ppl.title("Random points inside or outside a circle " + \
    "using the monte-carlo method")
ppl.plot(x_list_ic, y_list_ic, 'go', markersize=1)
ppl.plot(x_list_oc, y_list_oc, 'ro', markersize=1)
ppl.legend(["points inside the circle", "points outside the circle"])
ppl.xlim(-1,1)
ppl.ylim(-1,1)
ppl.xlabel("x")
ppl.ylabel("y")
pp.savefig()

# setup plot 2
ppl.figure("Random points under or above a gaussian using the monte-carlo " + \
    "method using " + str(steps) + " points")
ppl.title("Random points under or above a gaussian using\n the monte-carlo " + \
    "method using " + str(steps) + " points")
ppl.plot(x_list_ig_N, y_list_ig_N, 'go', markersize=1)
ppl.plot(x_list_og_N, y_list_og_N, 'ro', markersize=1)
ppl.legend(["points under the function", "points above the function"])
ppl.xlim(-1,1)
ppl.ylim(0,1)
ppl.xlabel("x")
ppl.ylabel("y")
pp.savefig()

# setup plot 3 (N = 0 excluded)
ppl.figure("The integrated value as a function of the number of points " + \
    "using the monte carlo method, error bars included")
ppl.title("The integrated value as a function of\n the number of points " + \
    "using the monte carlo method, error bars included")
ppl.errorbar(N_list, int_gaussian, int_gaussian_error, ecolor="black", \
    elinewidth=1)
ppl.xscale("log")
ppl.xlim(1,steps)
ppl.ylim(int_gaussian[len(int_gaussian)-1]-int_gaussian_error[0], \
    int_gaussian[len(int_gaussian)-1]+int_gaussian_error[0])
ppl.xlabel("N")
ppl.ylabel("integrated value")
pp.savefig()

# setup plot 4 (N = 0 excluded)
ppl.figure("The size of the error as a function of the number of points " + \
    "using the monte carlo method")
ppl.title("The size of the error as a function of\n the number of points " + \
    "using the monte carlo method")
ppl.semilogx(N_list, int_gaussian_error)
ppl.semilogx(N_list, N_sq_inverse_list)
ppl.legend(["numerical error", "1/sqrt(N)"])
ppl.xlim(1, steps)
ppl.ylim(0, 1)
ppl.xlabel("N")
ppl.ylabel("N error")
pp.savefig()

# plot show and pdf close
pp.close()
ppl.show()
