# Dani van Enk, 11823526
# randomnumber.py generate random numbers using the flat distribution and
#   gaussian distributions, it also finds the average and standard deviation

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
total_fd = 0
total_gd = 0
total_sq_fd = 0
total_sq_gd = 0
a = -5
b = 5
bins_per_tick = 10
pp = PdfPages(PATH + "/figures-randomnumber.pdf")

# lists
rd_list_fd = []
rd_list_gd = []

# Gaussian() is a gaussian function with standard deviation of 1
#   a center at x = 0
def Gaussian(x):
    return math.exp(-x**2/2)

# random_number() generates a number between a and b
def random_number(a,b):
    rd.seed()
    rd_nb = rd.random()
    diff = b-a
    return_value = rd_nb*diff + a
    return return_value

# monte_carlo() generates a random number between a and b using monte-carlo
def monte_carlo(a, b):
    rd_number = random_number(a, b)
    rd_f = random_number(0, Gaussian(0))
    f = Gaussian(rd_number)

    # if random function value is higher than function value, try again
    while rd_f > f:
        rd_number = random_number(a, b)
        rd_f = random_number(0, Gaussian(0))
        f = Gaussian(rd_number)
    
    return rd_number

# generate steps random numbers using the flat and gaussian distributions
#   and find total, squared total and append the random numbers
for i in range(0,steps):
    rd_number_fd = random_number(-5, 5)
    rd_number_gd = monte_carlo(-5, 5)

    print(i)
    total_fd += rd_number_fd
    total_gd += rd_number_gd
    total_sq_fd += (rd_number_fd)**2
    total_sq_gd += (rd_number_gd)**2
    rd_list_fd.append(rd_number_fd)
    rd_list_gd.append(rd_number_gd)

# calculate average and standard deviation
average_fd = total_fd/steps
average_gd = total_gd/steps
average_sq_fd = total_sq_fd/steps
average_sq_gd = total_sq_gd/steps
std_dv_fd = math.sqrt(average_sq_fd - average_fd**2)
std_dv_gd = math.sqrt(average_sq_gd - average_gd**2)

# print average and standard deviation of the flat and gaussian distributions
print("The average number after generating random numbers (flat distribution) " \
    "between " + str(a) + " and " + str(b) + " is " + str(average_fd) + \
      " with a standard deviation of " + str(std_dv_fd))
print("The average number after generating random numbers "+ \
    "(gassian distribution) between " + str(a) + " and " + str(b) + " is " + \
    str(average_gd) + " with a standard deviation of " + str(std_dv_gd))

# setup plot 1
ppl.figure("flat distribution using the random number method")
ppl.title("flat distribution using the random number method")
ppl.hist(rd_list_fd, bins=(b-a)*bins_per_tick, range=[a, b])
ppl.xlim(a, b)
ppl.xticks(range(a, b+1))
ppl.xlabel("random number value")
ppl.ylabel("count of random number value")
pp.savefig()

# setup plot 2
ppl.figure("gaussian distribution using the monte-carlo method")
ppl.title("gaussian distribution using the monte-carlo method")
ppl.hist(rd_list_gd, bins=(b-a)*bins_per_tick, range=[a, b])
ppl.xlim(a, b)
ppl.xticks(range(a, b+1))
ppl.xlabel("random number value")
ppl.ylabel("count of random number value")
pp.savefig()

# plot show and pdf close
pp.close()
ppl.show()
