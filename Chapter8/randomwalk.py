# Dani van Enk, 11823526
# randomwalk.py finds the paths of a dunken student, a particle in a box
# 	the average distance to (0,0) of 100 steps with/without a random stepsize,
# 	a wall, interaction between the particles
# It also calculates when the particles begin in the a left partition, what 
# 	the relaxation time is (left = right) and how long it takes for a particle
# 	to get to the right side, left and right are separated by the y-axis

# imports
from matplotlib import pyplot as ppl
import math
import os
import random as rd
from matplotlib.backends.backend_pdf import PdfPages

# Path to current file
PATH = os.path.dirname(os.path.abspath(__file__))

# variables
steps = 1000
relax_steps = 10000
sim_steps = 10
total_relaxation = 0
total_relaxation_sq = 0
total_first = 0
total_first_sq = 0
pp = PdfPages(PATH + "/figures-randomwalk.pdf")

# out_of_bounds() finds wall coords if (x,y) is out of bound,
# 	using (x,y) and the incoming angle 
def out_of_bounds(x,y, angle_in):
	# wall coords
	x_min = -10
	x_max = 10
	y_min = -10
	y_max = 10

	if x < x_min:
		y -= (x - x_min)*math.tan(angle_in)
		# print(x_min, y, "x_min")
		return True, x_min, y, 1
	elif x > x_max:
		y -= (x - x_max)*math.tan(angle_in)
		# print(x_max, y, "x_max")
		return True, x_max, y, 1
	elif y < y_min:
		x -= (y - y_min)/math.tan(angle_in)
		# print(x, y_min, "y_min")
		return True, x, y_min, 0
	elif y > y_max:
		x -= (y - y_max)/math.tan(angle_in)
		# print(x, y_max, "y_max")
		return True, x, y_max, 0
	else:
		return False, x, y, 0

# limits() finds minimum and maximum of a list and rounds it to whole given int
def limits(List, rounded_to):
	# variables
	maximum = 0
	minimum = 0

	# loop over all items in the list and finds min an max in list
	for items in List:
		if items > maximum:
			maximum = rounded_to*math.ceil(items/rounded_to)
		elif items < minimum:
			minimum = rounded_to*math.floor(items/rounded_to)
	
	return minimum, maximum

# random_step() finds the next step after (x,y) using given step size
def random_step(x, y, step_size):
	rd.seed()
	rd_angle = rd.random()*2*math.pi
	x += step_size*math.cos(rd_angle)
	y += step_size*math.sin(rd_angle)

	return x, y

# random_walk() makes a steps long random walk
# Wall, random stepsize and a LR partition of the box are enabled when True
def random_walk(steps, Wall, random_steps, LRpartition):
	# lists
	x_list = []
	y_list = []

	# setting initial coords for wall T/F and LRpartition T/F
	if Wall == False:
		x = 0
		y = 0
	elif Wall == True:
		x = 20*rd.random()-10
		y = 20*rd.random()-10
	if Wall == True and LRpartition == True:
		rd.seed()
		x = -10*rd.random()

	# appending initial coords to the lists
	x_list.append(x)
	y_list.append(y)

	# settig steps steps
	for step in range(0, steps):
		# step size setup according to input
		if random_steps == True:
			rd.seed()
			stepsize = rd.random()
		elif random_steps == False:
			stepsize = 1

		# set a step and find displacement and angle
		x, y = random_step(x, y, stepsize)
		dx = x - x_list[len(x_list)-1]
		dy = y - y_list[len(y_list)-1]
		angle_in = math.atan2(y_list[len(y_list)-1], x_list[len(x_list)-1])

		# check if still inside the wall, if wall is enabled
		if Wall == True:
			OoB, x_wall, y_wall, reflect_corr = out_of_bounds(x, y, angle_in)

		# while outside wall reflect from the wall if wall is enabled
		while Wall == True and OoB == True:
			# print("Out of bounds", x, y)
			dx = x - x_list[len(x_list)-1]
			dy = y - y_list[len(y_list)-1]
			dr = math.sqrt(dx**2 + dy**2)
			dr_wall =math.sqrt((x-x_wall)**2 + (y-y_wall)**2)
			angle_in = math.atan2(y_list[len(y_list)-1], x_list[len(x_list)-1])
			angle_out = math.pi*reflect_corr - angle_in
			x = x_wall + (dr-dr_wall)*math.cos(angle_out)
			y = y_wall + (dr-dr_wall)*math.sin(angle_out)
			OoB, x_wall, y_wall, reflect_corr = out_of_bounds(x, y, angle_in)

		# if position is not outside wall, append position to lists
		x_list.append(x)
		y_list.append(y)

	return x_list, y_list

# random_walk_with_interaction() makes a steps long random walk with
# 	interaction between the particles
# Wall, random stepsize and a LR partition of the box are enabled when True
def random_walk_with_interaction(particles, \
	steps, Wall, random_steps, LRpartition):
	# lists
	x_list_all = [[] for _ in range(0, particles)]
	y_list_all = [[] for _ in range(0, particles)]
	x_list = []
	y_list = []

	# create inital positions for the particles depending on input
	for particle in range(0, particles):
		if Wall == False:
			x = 0
			y = 0
		elif Wall == True:
			x = 20*rd.random()-10
			y = 20*rd.random()-10
		if Wall == True and LRpartition == True:
			rd.seed()
			x = -10*rd.random()

		# appending position to list with all position 
		# 	and list of current position per particle
		x_list_all[particle].append(x)
		y_list_all[particle].append(y)
		x_list.append(x)
		y_list.append(y)

	# find position of particles for steps steps
	for step in range(0, steps):
		for particle1 in range(0,particles):
			# setting stepsize according to input
			if random_steps == True:
				rd.seed()
				stepsize = rd.random()
			elif random_steps == False:
				stepsize = 1

			# generate new position for a particle
			# 	calculate displacement and angle
			x, y = random_step(x_list[particle1], y_list[particle1], stepsize)
			dx = x - x_list_all[particle1][len(x_list_all[particle1])-1]
			dy = y - y_list_all[particle1][len(y_list_all[particle1])-1]
			dr = math.sqrt(dx**2 + dy**2)
			angle_in = math.atan2(y_list_all[particle1] \
				[len(y_list_all[particle1])-1], x_list_all[particle1] \
				[len(x_list_all[particle1])-1])
		
			# check if still inside the wall, if wall is enabled
			if Wall == True:
				OoB, x_wall, y_wall, reflect_corr = \
					out_of_bounds(x, y, angle_in)

			# if particles hit each other during step bounce off
			for particle2 in range(0, len(x_list)):
				if x_list[particle1] == x_list[particle2] or \
					y_list[particle1] == y_list[particle2]:
					angle_out = math.pi - angle_in
					x += dr*math.cos(angle_out)
					y += dr*math.sin(angle_out)

			# while outside wall reflect from the wall if wall is enabled
			while Wall == True and OoB == True:
				# print("Out of bounds", x, y)
				dx = x - x_list_all[particle1][len(x_list_all[particle1])-1]
				dy = y - y_list_all[particle1][len(y_list_all[particle1])-1]
				dr = math.sqrt(dx**2 + dy**2)
				dr_wall = math.sqrt((x-x_wall)**2 + (y-y_wall)**2)
				angle_in = math.atan2(y_list_all[particle1] \
					[len(y_list_all[particle1])-1], x_list_all[particle1] \
					[len(x_list_all[particle1])-1])
				angle_out = math.pi*reflect_corr - angle_in
				x = x_wall + (dr-dr_wall)*math.cos(angle_out)
				y = y_wall + (dr-dr_wall)*math.sin(angle_out)
				OoB, x_wall, y_wall, reflect_corr = \
					out_of_bounds(x, y, angle_in)

			# if position is not outside wall, append position to lists
			x_list_all[particle1].append(x)
			y_list_all[particle1].append(y)
			x_list.append(x)
			y_list.append(y)

	return x_list_all, y_list_all

# av_distance_std() finds the average distance of particles particles 
# 	for steps steps from (0,0)
def av_distance_std(particles, steps, Wall, random_steps):
	# variabels and lists
	total_distance = 0
	total_distance_sq = 0
	position_list = []

	# generate the paths of particles particles
	for particle in range(0, particles):
		position_list.append(random_walk(steps, Wall, random_steps, False))

	# for each step find the distance and distance**2 for all particles
	for step in range(0, steps):
		for particle in range(0, particles):
			total_distance += math.sqrt(position_list[particle][0][step]**2 + \
				position_list[particle][1][step]**2)
			total_distance_sq += position_list[particle][0][step]**2 + \
				position_list[particle][1][step]**2

	# find the average distance and standard deviation
	average_distance = total_distance/(steps*particles)
	average_distance_sq = total_distance_sq/(steps*particles)
	distance_std = math.sqrt(average_distance_sq - average_distance**2)

	# print the result
	print("Using " + str(particles) + " non-interacting particles, Wall = " + \
		str(Wall) + " & random step size = " + str(random_steps))
	print("The average distance from (0,0) after " + str(steps) + " steps is " \
		+ str(average_distance) + " with a standard deviation of " + \
		str(distance_std))

# av_distance_std_with_interaction() finds the average distance of particles 
# 	particles for steps steps from (0,0) with interaction of the particles
def av_distance_std_with_interaction(particles, steps, Wall, random_steps):
	# variables and lists
	total_distance = 0
	total_distance_sq = 0
	position_list = []

	# creating a list of all positions of all the particles
	position_list = random_walk_with_interaction(particles, \
		steps, Wall, random_steps, False)

	# for each step find the distance and distance**2 for all particles
	for step in range(0, steps):
		for particle in range(0, particles):
			total_distance += math.sqrt(position_list[0][particle][step]**2 + \
				position_list[1][particle][step]**2)
			total_distance_sq += position_list[0][particle][step]**2 + \
				position_list[1][particle][step]**2

	# for each step find the distance and distance**2 for all particles
	average_distance = total_distance/(steps*particles)
	average_distance_sq = total_distance_sq/(steps*particles)
	distance_std = math.sqrt(average_distance_sq - average_distance**2)

	# print the result
	print("Using " + str(particles) + " interacting particles, Wall = " + \
		str(Wall) + " & random step size = " + str(random_steps))
	print("The average distance from (0,0) after " + str(steps) + " steps is " \
		+ str(average_distance) + " with a standard deviation of " + \
		str(distance_std))

# relaxation() finds the relaxation time for particles particles
#	it also finds the no of particles in a partitioned box in each partition
# all particles begin in the left partition
def relaxation(particles, steps, Wall, random_steps, LRpartition):
	# variables and lists
	left = 0
	right = 0
	first = 0
	position_list = []
	step_list = []
	left_list = []
	right_list = []

	# creating the paths of the particles
	for particle in range(0, particles):
		position_list.append( \
			random_walk(steps, Wall, random_steps, LRpartition))

	# for each step find how many particles are in the left
	# 	and right partition
	for step in range(0, steps):
		# initial right/left values
		left = 0
		right = 0

		# find the if a particle is left of right in the box
		for particle in range(0, particles):
			if position_list[particle][0][step] > 0:
				# find the first particle to get to the right
				if first == 0:
					first = step
				right += 1
			elif position_list[particle][0][step] < 0:
				left += 1

		# append the steps of iteration and the no of particles
		# 	in the left and right partition
		step_list.append(step)
		left_list.append(left)
		right_list.append(right)

		# find the relaxation time
		if left == right:
			relaxation_time = step
			break

	return step_list, left_list, right_list, relaxation_time, first

# get the paths of a drunken student (equal stepsize and no wall),
# 	a particle in a box and the limits of the student
x_list_ess, y_list_ess = random_walk(steps, False, False, False)
xmin_ess, xmax_ess = limits(x_list_ess, 10)
ymin_ess, ymax_ess = limits(y_list_ess, 10)
x_list_rss, y_list_rss = random_walk(steps, True, False, False)

# find the average distance from (0,0) for 100 particles with a random
# 	and equal stepsize after 10 and 1000 steps
av_distance_std(100, 10, False, False)
av_distance_std(100, 10, False, True)
av_distance_std(100, 1000, False, False)
av_distance_std(100, 1000, False, True)

# whiteline
print()

# find the average distance from (0,0) for 100 particles in a box with
# 	a random and equal stepsize after 10 and 1000 steps
av_distance_std(100, 10, True, False)
av_distance_std(100, 10, True, True)
av_distance_std(100, 1000, True, False)
av_distance_std(100, 1000, True, True)

# whiteline
print()

# repeat the relaxation function to rule out statistic fluctuations
for simulation in range(0, sim_steps):
	step_list, left_list, right_list, relaxation_time, first = \
		relaxation(100, relax_steps, True, False, True)
	total_relaxation += relaxation_time
	total_first += first
	total_relaxation_sq += relaxation_time**2
	total_first_sq += first**2

# whiteline
print()

# find the average distance from (0,0) for 100 particles in a box with
# 	a random and equal stepsize and interaction after 10 steps
av_distance_std_with_interaction(100, 10, True, False)
av_distance_std_with_interaction(100, 10, True, True)

# find the average relaxation and first right times
# 	their standard deviation
average_relaxation = total_relaxation/sim_steps
average_first = total_first/sim_steps
average_relaxation_sq = total_relaxation_sq/sim_steps
average_first_sq = total_first_sq/sim_steps
relaxation_std = math.sqrt(average_relaxation_sq - average_relaxation**2)
first_std = math.sqrt(average_first_sq - average_first**2)

# print the relaxation result
print("The relaxation time for 100 particles is " + str(average_relaxation) + \
	" with a standard deviation of " + str(relaxation_std) + \
	", using an average over " + str(sim_steps) + " iterations")
print("The time for 100 particles, at which the first particle arrives in " + \
	"the right partition is " + str(average_first) + \
	" with a standard deviation of " + str(first_std) + \
	", using an average over " + str(sim_steps) + " iterations")

# plot-1 setup
ppl.figure("Path of a single particle for " + str(steps))
ppl.title("Path of a single particle for " + str(steps))
ppl.plot(x_list_ess, y_list_ess, linewidth=.5)
ppl.xlim(xmin_ess, xmax_ess)
ppl.ylim(ymin_ess, ymax_ess)
ppl.xlabel("x")
ppl.ylabel("y")
pp.savefig()

# plot-2 setup
ppl.figure("Path of a single particle for " + str(steps) + \
	" with a wall & random steps size")
ppl.title("Path of a single particle for " +
          str(steps) + " with a wall & random steps size")
ppl.plot(x_list_rss, y_list_rss, linewidth=.5)
ppl.xlim(-10, 10)
ppl.ylim(-10, 10)
ppl.xlabel("x")
ppl.ylabel("y")
pp.savefig()

# plot-3 setup
ppl.figure("No. particles in the left & right partition of " + \
	"the box as a function of steps")
ppl.title("No. particles in the left & right partition of " + \
	"the box as a function of steps")
ppl.plot(step_list, left_list, color="Red")
ppl.plot(step_list, right_list, color="Blue")
ppl.legend(["Particles in Left partition", "Particles in Right partition"])
ppl.xlim(0, step_list[len(step_list)-1])
ppl.ylim(0, 100)
ppl.xlabel("iteration step")
ppl.ylabel("no. of particles")
pp.savefig()

# plot show and pdf close
pp.close()
ppl.show()
