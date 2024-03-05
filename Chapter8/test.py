# Dani van Enk, 11823526
# randomnumber.py generate random numbers using the flat distribution and
#   gaussian distributions, it also finds the average and standard deviation

# imports
import math
import os
import random as rd

# variables
steps = 1000
std_steps = 100
xmin = -1
xmax = 1
ymin = -1
ymax = 1
area = abs(xmax-xmin)*abs(ymax-ymin)
area_g = abs(xmax-xmin)*abs(ymax-0)

# lists
N_list = [0]
N_sq_inverse_list = [1]
int_gaussian = [0]
int_gaussian_error = [0]

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
def monte_carlo_gaussian(xmin, xmax, ymin, ymax, steps):
    step = 0

    x_list_inside = []
    x_list_outside = []
    y_list_inside = []
    y_list_outside = []

    while step < steps:
        rd_x = random_number(xmin, xmax)
        rd_y = random_number(ymin, ymax)
        rd_f = Gaussian(rd_x, 1, 0)

        # check if the points are inside or outside the circle of radius 1
        if rd_f < rd_y:
            x_list_outside.append(rd_x)
            y_list_outside.append(rd_y)
        if rd_f > rd_y:
            x_list_inside.append(rd_x)
            y_list_inside.append(rd_y)

        # print("gaussian step", step)

        step += 1

    return x_list_inside, x_list_outside, y_list_inside, y_list_outside

x_list_ig, x_list_og, y_list_ig, y_list_og = \
    monte_carlo_gaussian(xmin, xmax, 0, ymax, steps)

for N in range(1,steps, 1000):
    total = 0
    total_sq = 0
    for N_std in range(0,std_steps):
        x_list_ig_N, x_list_og_N, y_list_ig_N, y_list_og_N = \
            monte_carlo_gaussian(xmin, xmax, 0, ymax, N)
        integralN = area_g*(len(x_list_ig_N)/(len(x_list_ig_N) + len(x_list_og_N)))

        total += integralN
        total_sq += integralN**2

    # calculate average and standard deviation
    average = total/std_steps
    average_sq = total_sq/std_steps
    std_dv = math.sqrt(average_sq - average**2)

    N_list.append(N)
    N_sq_inverse_list.append(1/math.sqrt(N))
    int_gaussian.append(average)
    int_gaussian_error.append(std_dv)

integral = area_g*(len(x_list_ig)/(len(x_list_ig) + len(x_list_og)))

print("The integral of a gausian function with standard deviation of 1 and " + \
    " center at 0 from -1 to 1 is " + str(int_gaussian[len(int_gaussian)-1]))
print("The calculated value of the integeral is " + \
      str(abs((int_gaussian[len(int_gaussian)-1]-math.sqrt(2*math.pi)*math.erf(1/math.sqrt(2))) /
    math.sqrt(2*math.pi)*math.erf(1/math.sqrt(2)))*100) + \
    " percent from the real value")
