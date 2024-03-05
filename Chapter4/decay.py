# Dani van Enk, 11823526
# decay.py calculates the exponential decay of Bromine numerically using
#   the Euler method and the Runge-Kutta method to the 4th order and finds
#   the error in the numerical values with respect tot the analytical values.

# imports
from matplotlib import pyplot as ppl
import math
import os
from matplotlib.backends.backend_pdf import PdfPages

# Path to current file
PATH = os.path.dirname(os.path.abspath(__file__))
# variables
N0 = 1
ENBr = N0
ENSe = 0
RKNBr = N0
RKNSe = 0
step = 0
Thalf = 16.1*3600
lamb = math.log(2,math.e)/(Thalf)
time_step = 0.01*Thalf
no_of_steps = 750
time = 0
pp = PdfPages(PATH + "/figures-decay.pdf")

# lists
time_list = [0]
EBromine, RKBromine = [ENBr], [RKNBr]
ESelenium, RKSelenium = [ENSe], [RKNSe]
AnalytBromine, AnalytSelenium = [N0], [0]
EDifferenceBromine, RKDifferenceBromine = [0], [0]

# formula for exponential decay
def N(t):
    return N0*math.exp(-lamb*t)

# use Euler method to calculate Bromine and Selenium left after decay
def Euler(NBr, NSe, N):
    NBr += time_step*-lamb*NBr
    NSe = (N0 - NBr)
    return NBr, NSe

# use Runge-Kutta method to calculate Bromine and Selenium left after decay
def Runge_Kutta(NBr, NSe, N):
    # assign the 4 orders for Runge-Kutta method
    k1 = time_step*-lamb*N
    k2 = time_step*-lamb*(N + k1/2)
    k3 = time_step*-lamb*(N + k2/2)
    k4 = time_step*-lamb*(N + k3)

    # actually do the Runge-Kutta method tot caclulate quantities
    NBr += (k1 + 2*k2 + 2*k3 + k4)/6
    NSe = (N0 - NBr)

    return NBr, NSe

# creating time axis
while step <= no_of_steps:
    # make time list
    time_list.append(time/3600)

    # make data list for the Euler method
    EBromine.append(ENBr)
    ESelenium.append(ENSe)

    # make data list for Runge-Kutta method
    RKBromine.append(RKNBr)
    RKSelenium.append(RKNSe)

    # make data list for analytical solution for decay
    AnalytBromine.append(N(time))
    AnalytSelenium.append(N0 - N(time))

    # make list of error in numerical
    EDifferenceBromine.append(abs(ENBr - N(time))/N(time))
    RKDifferenceBromine.append(abs(RKNBr - N(time))/N(time))

    # variables from functions
    ENBr, ENSe = Euler(ENBr, ENSe, N(time))
    RKNBr, RKNSe = Runge_Kutta(RKNBr, RKNSe, N(time))

    time += time_step
    step += 1

# plot-1 setup and save to pdf
ppl.figure("Numerical model of Bromine decay into Selenium (Euler method)")
ppl.title("Bromine decay to Selenium (Euler Method)")
ppl.plot(time_list, EBromine, color="Red")
ppl.plot(time_list, ESelenium, color="Blue")
ppl.xlim(0, no_of_steps*time_step/3600)
ppl.ylim(0, 1)
ppl.ylabel("Quantity/N0")
ppl.xlabel("time in hours")
ppl.legend(["Quantity of Bromine", "Quantity of Selenium"])
pp.savefig()

# plot-2 setup and save to pdf
ppl.figure("Error in numerical vs analytical for Euler method")
ppl.title("Error in numerical vs analytical for Euler method")
ppl.plot(time_list, EDifferenceBromine, color="Green")
ppl.xlim(0, no_of_steps*time_step/3600)
ppl.ylim(0)
ppl.ylabel("error in numerical")
ppl.xlabel("time in hours")
pp.savefig()

# plot-3 setup and save to pdf
ppl.figure("Error in numerical vs analytical for Runge-Kutta method")
ppl.title("Error in numerical vs analytical for Runge-Kutta method", y=1.05)
ppl.plot(time_list, RKDifferenceBromine, color="Yellow")
ppl.xlim(0, no_of_steps*time_step/3600)
ppl.ylim(0, 3.5e-9)
ppl.ylabel("error in numerical")
ppl.xlabel("time in hours")
pp.savefig()

# plot show and pdf close
pp.close()
ppl.show()
