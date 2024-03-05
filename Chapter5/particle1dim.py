# Dani van Enk, 11823526
# particle1dim.py uses the Adam-Bash method to calculate the position, momentum
# 	and energies of a 1 dimensional particle moving towards a potential
# two cases are shown, if E_kin < E_pot (**_min is used for such names) and
# 	if E_kin > E_pot (**_max is used for such names)

# imports
from matplotlib import pyplot as ppl
import math
import os
from matplotlib.backends.backend_pdf import PdfPages

# Path to current file
PATH = os.path.dirname(os.path.abspath(__file__))

# variables
V0 = 5
speed_min = 1
speed_max = 4
xpos = 5
mass = 1
time = 0
steps = 1000
time_max = 10
time_step = time_max/steps
momentum_min = mass*speed_min
momentum_max = mass*speed_max
x_min = 0
x_max = 0
pp = PdfPages(PATH + "/figures-particle1dim.pdf")

# lists
time_list = []
momentum_min_list = []
momentum_max_list = []
x_min_pos = []
x_max_pos = []
E_kin_min = []
E_kin_max = []
E_pot_min = []
E_pot_max = []
E_tot_min = []
E_tot_max = []

# Gauss() defines a gausian function
def Gauss(C,x,D):
	return C*math.exp(-(x)**2*D)

# Potential(x) defines the potential energy at x
def Potential(x):
	return Gauss(V0,x-xpos,1/2)

# Force(x) defines the force at x
def Force(x):
	return (x-xpos)*Potential(x)

# use the Adams-Bash method to calculate the next step
def Adams_Bash(number,fi,fi_1,h):
	return number + h*(3/2*fi-1/2*fi_1)

# calculate the time, position, momentum and energies at each timestep
for index in range(0, steps):
	# append time, position and momentum for speed_min and speed_max to their
	# 	respective lists
	time_list.append(time)

	x_min_pos.append(x_min)
	x_max_pos.append(x_max)

	momentum_min_list.append(momentum_min)
	momentum_max_list.append(momentum_max)

	# append total, kinetic and potential energies to their respective lists
	E_kin_min.append(momentum_min**2/2)
	E_kin_max.append(momentum_max**2/2)

	E_pot_min.append(Potential(x_min))
	E_pot_max.append(Potential(x_max))

	E_tot_min.append(momentum_min**2/2 + Potential(x_min))
	E_tot_max.append(momentum_max**2/2 + Potential(x_max))

	# next timestep
	time += time_step

	# calculate new positions and momentum
	x_min = Adams_Bash(x_min, momentum_min_list[index], \
		momentum_min_list[index-1], time_step)
	x_max = Adams_Bash(x_max, momentum_max_list[index], \
		momentum_max_list[index-1], time_step)

	momentum_min = Adams_Bash(momentum_min, Force(x_min_pos[index]), \
		Force(x_min_pos[index-1]), time_step)
	momentum_max = Adams_Bash(momentum_max, Force(x_max_pos[index]), \
		Force(x_max_pos[index-1]), time_step)

# plot-1 setup and save to pdf
ppl.figure("Position plot against the time for E_kin < E_pot")
ppl.title("Position plot against the time for E_kin < E_pot")
ppl.plot(time_list, x_min_pos, color="Green")
ppl.xlim(0, time_max)
ppl.ylim(-3, 3)
ppl.xlabel("time in seconds")
ppl.ylabel("position in meter")
pp.savefig()

# plot-2 setup and save to pdf
ppl.figure("Momentum plot against the time for E_kin < E_pot")
ppl.title("Momentum plot against the time for E_kin < E_pot")
ppl.plot(time_list, momentum_min_list, color="Red")
ppl.xlim(0, time_max)
ppl.ylim(-1, 1)
ppl.xlabel("time in seconds")
ppl.ylabel("momentum in kg*m/s")
pp.savefig()

# plot-3 setup and save to pdf
ppl.figure("Position plot against the time for E_kin > E_pot")
ppl.title("Position plot against the time for E_kin > E_pot")
ppl.plot(time_list, x_max_pos, color="Green")
ppl.xlim(0, time_max)
ppl.ylim(0, 40)
ppl.xlabel("time in s")
ppl.ylabel("position in m")
pp.savefig()

# plot-4 setup and save to pdf
ppl.figure("Momentum plot against the time for E_kin > E_pot")
ppl.title("Momentum plot against the time for E_kin > E_pot")
ppl.plot(time_list, momentum_max_list, color="Red")
ppl.xlim(0, time_max)
ppl.ylim(0, 4.5)
ppl.xlabel("time in s")
ppl.ylabel("momentum in kg*m/s")
pp.savefig()

# plot-5 setup and save to pdf
ppl.figure("E_kin, E_pot and E_tot in the time for E_kin < E_pot")
ppl.title("E_kin, E_pot and E_tot in the time for E_kin < E_pot")
ppl.plot(time_list, E_kin_min, color="Red")
ppl.plot(time_list, E_pot_min, color="Grey")
ppl.plot(time_list, E_tot_min, color="Yellow")
ppl.xlim(0, time_max)
ppl.ylim(0, .6)
ppl.legend(["E_kin", "E_pot", "E_tot"])
ppl.xlabel("time in s")
ppl.ylabel("energy in J")
pp.savefig()

# plot-6 setup and save to pdf
ppl.figure("E_kin, E_pot and E_tot in the time for E_kin > E_pot")
ppl.title("E_kin, E_pot and E_tot in the time for E_kin > E_pot")
ppl.plot(time_list, E_kin_max, color="Red")
ppl.plot(time_list, E_pot_max, color="Grey")
ppl.plot(time_list, E_tot_max, color="Yellow")
ppl.xlim(0, time_max)
ppl.ylim(0, 9)
ppl.legend(["E_kin", "E_pot", "E_tot"])
ppl.xlabel("time in s")
ppl.ylabel("energy in J")
pp.savefig()

# plot show and pdf close
pp.close()
ppl.show()
