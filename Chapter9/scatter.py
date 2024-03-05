# Dani van Enk, 11823526
# scatter.py calculates the paths and energies of particles using the
# 	using the rutherford scattering idea and using 2 moving particles
# it also uses the electrodynamic forces to find the path of 10 particles in
#	3D, as well as the path of inner the planets in the solar system using the
# 	gravitational forces
# because of the animations the code may take about 5-10 mins to run

# imports
from matplotlib import pyplot as ppl
from matplotlib import animation
import math
import os
import random as rd
from matplotlib.backends.backend_pdf import PdfPages
from fractions import Fraction as F
from mpl_toolkits.mplot3d import Axes3D
import imageio as imio

# Path to current file
PATH = os.path.dirname(os.path.abspath(__file__))

# variables
day = 3600*24
year = int(365.25*day)
t = 0
tmax = 25
tmax_dd = 5
tmax_3D = 10
tmax_CB3D = 3*year
steps = 100000
steps_3D = 100000
steps_CB3D = 10000
epsilon_0 = 8.854e-12
e = 1.60218e-19
c = 299792458
G = 6.67408e-11
mprot = 1.672623e-27
melec = 9.10938356e-31
rad_convert = math.pi/180
ly = 9.461e15
au = 1.496e11
dt = F(tmax,steps)
dt_dd = F(tmax_dd, steps)
dt_3D = F(tmax_3D, steps_3D)
dt_CB3D = F(tmax_CB3D, steps_CB3D)
no_of_particles3D = 10
correct_input = False
pp = PdfPages(PATH + "/figures-scatter.pdf")

# Get input parameters from user for the rutherford scattering
while correct_input == False:
	user_input = input("For a test particle an static pariticle please give:" + \
		"\nthe x-position of the dynamic particle xd in meters," + \
		"\nthe position of the static particle (xs,ys) in meters," + \
		"\nthe x-distance between the test and the static particle in m" + \
			"the speed of the dynamic particle in m/s"
		"\nthe mass in kg and charge in coulomb of both particles" + \
		"\nf.Ex: 0,5,5,0,1,1.672623e-27,1.67263e-27," + \
			"1.60218e-19,1.60218e-19\n")

	# Exclude any non number input or in the wrong format
	try:
		# separate parameters with ","
		# 	and make a list of numbers for the parametes
		user_input_list = [float(element) for element in user_input.split(",")]
		# make sure the user inputs 9 values
		if len(user_input_list) != 9:
			print("LengthError: the input has not enough values, " +
						"make sure to input 10 values in the given format")
			correct_input = False
		else:
			correct_input = True
	except ValueError:
		print("FormatError: Please try again in the correct format")
		correct_input = False

# the particle class defines a particle with inital position (x_init,y_init, z_init),
# 	speed, angle of that speed, mass and charge
#	it also has a function to describe the 2D motion of such charged particle
#	it also has a function to describe the 3D motion of such charged particle
class Particle:
	# __init() contains the begin parameters, those parameters are
	# 	number of particles, position (x_init,y_init,z_init), speed,
	# 	angle of the speed, mass and charge
	def __init__(self, i, x_init, y_init, z_init, speed, speed_polangle, \
		speed_aziangle, mass, charge):
		self.i = i
		self.x = x_init
		self.y = y_init
		self.z = z_init
		self.vx = speed*math.sin(speed_polangle*rad_convert)* \
			math.cos(speed_aziangle*rad_convert)
		self.vy = speed*math.sin(speed_polangle*rad_convert)* \
			math.sin(speed_aziangle*rad_convert)
		self.vz = speed*math.cos(speed_polangle*rad_convert)
		self.v = math.sqrt(self.vx**2 + self.vy**2 + self.vz**2)
		self.m = mass
		self.Q = charge

	# move2D() describes the motion of a charged particle in 2D given the
	# 	interacting	other particles, the timestep. There is also the
	# 	posibility to have 1 particle be static
	def move2D(self, particles, dt, static):
		# variables
		ax = 0
		ay = 0
		Ekin = 0
		Epot = 0

		# for rutherford scattering static = True
		# 	for moving particles static = False
		if static == True:
			Ekin += 1/2*self.m*(self.vx**2+self.vy**2)
		elif static == False:
			for particle in particles:
				Ekin += 1/2*particle.m*particle.v**2
			
		# calculate the acceleration taking the other particles into account
		# also calculating the potential Energy
		for particle in particles:
			if particle.i != self.i:
				dx = self.x - particle.x
				dy = self.y - particle.y
				angle = math.atan2(dy, dx)
				dr = math.sqrt((dx)**2 + (dy)**2)
				ax += 1/(4*math.pi*epsilon_0)*(particle.Q*self.Q)/ \
					(self.m*dr**2)*math.cos(angle)
				ay += 1/(4*math.pi*epsilon_0)*(particle.Q*self.Q)/ \
					(self.m*dr**2)*math.sin(angle)

				Epot += 1/(4*math.pi*epsilon_0)* \
					(particle.Q*self.Q)/dr

		# change speed and position according to the acelleration calculated above
		self.x += self.vx*dt + 1/2*ax*dt**2
		self.y += self.vy*dt + 1/2*ay*dt**2
		self.vx += ax*dt
		self.vy += ay*dt
		self.v = math.sqrt(self.vx**2 + self.vy**2)

		return self.x, self.y, Ekin, Epot

	# move3D() describes the motion of a charged particle in 3D given the
	# 	interacting other particles, the timestep.
	def move3D(self, particles, dt):
		# variables
		ax = 0
		ay = 0
		az = 0
		Epot = 0
		smallestr = 10

		# calculate the total kinetic enery for 1 step
		# for particle in particles:
		Ekin = 1/2*self.m*self.v**2

		# calculate the acceleration taking the other particles into account
		# also calculating the potential Energy
		for particle in particles:
			if self.i != particle.i:
				dx = self.x - particle.x
				dy = self.y - particle.y
				dz = self.z - particle.z
				dr = math.sqrt(dx**2 + dy**2 + dz**2)
				polangle = math.acos(dz/dr)
				aziangle = math.atan2(dy, dx)
				ax += 1/(4*math.pi*epsilon_0)*(particle.Q*self.Q) / \
					(self.m*dr**2)*math.sin(polangle)*math.cos(aziangle)
				ay += 1/(4*math.pi*epsilon_0)*(particle.Q*self.Q) / \
					(self.m*dr**2)*math.sin(polangle)*math.sin(aziangle)
				az += 1/(4*math.pi*epsilon_0)*(particle.Q*self.Q) / \
					(self.m*dr**2)*math.cos(polangle)

				# get the potential energy without counting double
				if self.i < particle.i:
					Epot += 1/(4*math.pi*epsilon_0) * \
						(particle.Q*self.Q)/dr

				# get the smallest distance between bodies
				if dr < smallestr:
					smallestr = dr

		# change speed and position according to the acelleration calculated above
		self.x += self.vx*dt + 1/2*ax*dt**2
		self.y += self.vy*dt + 1/2*ay*dt**2
		self.z += self.vz*dt + 1/2*az*dt**2
		self.vx += ax*dt
		self.vy += ay*dt
		self.vz += az*dt
		self.v = math.sqrt(self.vx**2 + self.vy**2 + self.vz**2)

		return self.x, self.y, self.z, Ekin, Epot, smallestr

# the Celestial_body class defines a Celestial_body with inital position
# 	(x_init,y_init, z_init), speed, angle of that speed, mass and charge
#	it also has a function to describe the 2D motion of such Celestial_body
#	it also has a function to describe the 3D motion of such Celestial_body
class Celestial_body:
	# __init() contains the begin parameters, those parameters are
	# 	number of particles, position (x_init,y_init,z_init), speed,
	# 	angle of the speed, mass, radius and  density
	def __init__(self, i, x_init, y_init, z_init, speed, speed_polangle, \
		speed_aziangle, mass, radius):
		self.i = i
		self.x = x_init
		self.y = y_init
		self.z = z_init
		self.vx = speed*math.sin(speed_polangle*rad_convert) * \
			math.cos(speed_aziangle*rad_convert)
		self.vy = speed*math.sin(speed_polangle*rad_convert) * \
			math.sin(speed_aziangle*rad_convert)
		self.vz = speed*math.cos(speed_polangle*rad_convert)
		self.v = math.sqrt(self.vx**2 + self.vy**2 + self.vz**2)
		self.m = mass
		self.R = radius
		self.dens = self.m/(4/3*math.pi*self.R**3)

	# Energy_Accel() calculates the energy and acceleration of the celestial
	# 	object according to the gravitational forces
	def Energy_Accel(self, bodies):
		# variables
		ax = 0
		ay = 0
		az = 0
		Ekin = 0
		Epot = 0
		smallestr = 100*au

		# calculate the total kinetic enery for 1 step
		Ekin += 1/2*self.m*self.v**2

		# calculate the acceleration taking the other bodies into account
		# also calculating the potential Energy
		for body in bodies:
			if self.i != body.i:
				dx = body.x - self.x
				dy = body.y - self.y
				dz = body.z - self.z
				dr = math.sqrt(dx**2 + dy**2 + dz**2)
				polangle = math.acos(dz/dr)
				aziangle = math.atan2(dy, dx)

				F = (G*self.m*body.m)/(dr**2)
				ax += F/self.m*math.sin(polangle) * \
					math.cos(aziangle)
				ay += F/self.m*math.sin(polangle) * \
					math.sin(aziangle)
				az += F/self.m*math.cos(polangle)

				# get the potential energy without counting double
				if self.i < body.i:
					Epot += -(G*self.m*body.m)/dr

				# get the smallest distance between bodies
				if dr < smallestr:
					smallestr = dr

		return ax, ay, az, Ekin, Epot, smallestr

	# move3D() describes the motion of a celestial body given the acceleration
	# 	and energy of the body
	def move3D(self, bodies, dt):
		# variables
		ax, ay, az, Ekin, Epot, smallestr = self.Energy_Accel(bodies)

		# change speed and position according to the acelleration calculated above
		self.x += self.vx*dt + 1/2*ax*dt**2
		self.y += self.vy*dt + 1/2*ay*dt**2
		self.z += self.vz*dt + 1/2*az*dt**2
		self.vx += ax*dt
		self.vy += ay*dt
		self.vz += az*dt
		self.v = math.sqrt(self.vx**2 + self.vy**2 + self.vz**2)

		return self.x, self.y, self.z, Ekin, Epot, smallestr

# deflection_angle() calculates the deflection angle according to the literature
def deflection_angle(md, ms, Qd, Qs, v0, b):
	Z1 = md/mprot
	Z2 = ms/mprot
	return 2*math.atan2(Z1*Z2*Qs*Qd, 4*math.pi*epsilon_0*md*v0**2*b)

# sd_2particle_numerical() finds the path of a dynamic particle 
#   at position (xd, yd) as it moves with a speed towards a static particle
#   at (xs, ys), both have a mass of m and an charge of Q. All variables
#   indicated with an s are for the static particle and all indicated with 
#   a d are for the dynamic particle
# 	the polar angle is 90 and z = 0 because this is 2D
def sd_2particle(xd, xs, ys, b, speed, md, ms, Qd, Qs):
	# initial kinetic and potential energy
	Ekin = 1/2*md*(speed)**2
	Epot = 1/(4*math.pi*epsilon_0)*(Qd*Qs)/math.sqrt((xd-xs)**2+(b)**2)

	# lists
	x_list = [xd]
	y_list = [ys + b]
	Ekin_list = [Ekin]
	Epot_list = [Epot]
	Etot_list = [Ekin+Epot]
	t_list = [0]

	# create a moving charged particle and a static particle
	particles_sd = [Particle(0,xd,ys + b,0,speed,90,0,md,Qd), \
		Particle(1,xs,ys,0,0,90,0,ms,Qs)]

	# calculate the position and kinetic, potential and total energy
	for step in range(0, steps):
		x, y, Ekin, Epot = particles_sd[0].move2D(particles_sd, dt, True)
		x_list.append(x)
		y_list.append(y)
		t_list.append(step*dt)
		Ekin_list.append(Ekin)
		Epot_list.append(Epot)
		Etot_list.append(Ekin+Epot)

	# printing the analytic and numerical values from the deflection angle
	print(str(deflection_angle(md,ms,Qd,Qs,speed,b)/rad_convert) + \
		" degrees is the angle from the analytics for the following parameters" + \
		"\nmd = " + str(md) + ", ms = " + str(ms) + ", Qd = " + str(Qd) + \
		", Qs = " + str(Qs) + ", speed = " + str(speed) + ", b = 0" + str(b))
	print(str(math.atan2(y_list[len(y_list)-1]-y_list[len(y_list)-2], \
		x_list[len(x_list)-1]-x_list[len(x_list)-2])/rad_convert) + \
		" degrees is the numerical found angle for the following parameters" +
		"\nmd = " + str(md) + ", ms = " + str(ms) + ", Qd = " + str(Qd) + \
		", Qs = " + str(Qs) + ", speed = " + str(speed) + ", b = 0" + str(b))

	"""
	Currently this function do any kind for rutherford scatering,
		for b=0 it is still not intirely correct. Also if the position
		between the charge and test charge at the begining and the end
		not the same is, a violation of conservation of energy occurs.
		So in that case this function is only good enough for an approximation
	"""

	return x_list, y_list, t_list, \
		Ekin_list, Epot_list, Etot_list, particles_sd[1]

# d_particles_2D() the path of 2 moving charged particles, particle 1
#   at position (x1, y1) as it moves with a speed towards particle 2
#   at (x2, y2), with mass of m1 and m2 and a charges of Q1 and Q2.
# 	All variables indicated with an 1 are for particle 1 and all indicated with
#   2 are for particle 2, the polar angel is 90 and z = 0 because this is 2D
def d_particles_2D(x1, y1, x2, y2, v1, angle1, v2, angle2, m1, m2, Q1, Q2):
	# lists
	x_list = [[], []]
	y_list = [[], []]
	Ekin_list = []
	Epot_list = []
	Etot_list = []
	t_list = []

	# creating 2 moving particles
	particles_d = [Particle(0,x1,y1,0,v1,90,angle1,m1,Q1), \
		Particle(1,x2,y2,0,v2,90,angle2,m2,Q2)]

	# looping for each step and calculating the paths
	# 	 energies of the particles
	for step in range(0, steps):
		# variables
		Ekin = 0
		Epot = 0

		# calculating the paths and energy parts for each particle
		# 	and appending those to their respective lists
		for particle in particles_d:
			x, y, Ekin_part, Epot_part = \
				particle.move2D(particles_d, dt_dd, False)
			x_list[particle.i].append(x)
			y_list[particle.i].append(y)
			Ekin += Ekin_part
			Epot += Epot_part

		# appending to lists
		t_list.append(step*dt_dd)
		Ekin_list.append(Ekin)
		Epot_list.append(Epot)
		Etot_list.append(Ekin+Epot)

	return x_list, y_list, t_list, Ekin_list, Epot_list, Etot_list

# moving_particle3D() creates a user difined number of parrtilces and calculates
# 	their paths and energies
def moving_particles3D(particles):
	# lists
	x_list = [[] for _ in range(0, particles)]
	y_list = [[] for _ in range(0, particles)]
	z_list = [[] for _ in range(0, particles)]
	Ekin_list = []
	Epot_list = []
	Etot_list = []
	t_list = []
	charges = [-e, e]
	particle_list = []

	# variables
	step_size = 1
	t_3D = 0

	# create particles with random inital parameters
	for particle in range(0, particles):
		# random variables
		x = rd.random()
		y = rd.random()
		z = rd.random()
		v = rd.random()
		polangle = 180*rd.random()
		aziangle = 360*rd.random()
		mass = mprot
		charge = rd.choice(charges)

		particle_list.append(Particle(particle, x, y, z, v, aziangle, \
			polangle, mass, charge))

	# calculate the position and energy for each particles
	for step in range(0, steps_3D):
		# variables
		Ekin = 0
		Epot = 0
		sr = 100

		# create the path of the particles,
		# 	it also calculates the total energies per step
		for particle in particle_list:
			x, y, z, Ekin_part, Epot_part, dr = \
				particle.move3D(particle_list, step_size*dt_3D)
			x_list[particle.i].append(x)
			y_list[particle.i].append(y)
			z_list[particle.i].append(z)

			# get the shortest distance between particles for this step
			if dr < sr:
				sr = dr

			# add the energies per particle
			Ekin += Ekin_part
			Epot += Epot_part

		# change the timestep according to the distance between particles
		if sr < .8:
			step_size = .1
		if sr < .5:
			step_size = .01
		if sr < .3:
			step_size = .001
		if sr < .1:
			step_size = .0001
		else:
			step_size = 1

		# get the time
		t_3D += step_size*step*dt_3D

		# appending to lists
		t_list.append(t_3D)
		Ekin_list.append(Ekin)
		Epot_list.append(Epot)
		Etot_list.append(Ekin+Epot)

	"""
	Sometimes the total energy is not constant this is because if the 
	particles come in really close proximity, the kinetic energy grows bigger
	but the potential energy not large enough to cancel this out. This
	causes the total energy to jump when 2 opposed charged particles 
	come in close proximity

	UPDATE: It should work now, but sometimes has a little hickup (small change)
	in the total energy
	"""

	return x_list, y_list, z_list, t_list, Ekin_list, Epot_list, Etot_list, \
		particle_list

# celestial_movement3D() describes the path of the planets in the solar system 
# 	around the sun
def celestial_movement3D():
	# object list
	bodies_list = [Celestial_body(0,0,0,0,0,0,0,2e30,69551e4), \
		Celestial_body(1,698169e5,0,0,47.4e3,7.005,90,3.285e23,2439.7e3), \
		Celestial_body(2,108939e6,0,0,35.02e3,3.39458,90,4.8675e24,6051.8e3), \
		Celestial_body(3,1521e8,0,0,29.78e3,7.155,90,5.97237e24,6371.0e3), \
		Celestial_body(4,2492e8,0,0,24.007e3,1.850,90,6.4171e23,3389.5e3)]

	# lists
	x_list = [[] for _ in range(0, len(bodies_list))]
	y_list = [[] for _ in range(0, len(bodies_list))]
	z_list = [[] for _ in range(0, len(bodies_list))]
	Ekin_list = []
	Epot_list = []
	Etot_list = []
	t_list = []

	# variables
	step_size = 1
	t_CB3D = 0

	# calculate the position and energy for each celestial body
	for step in range(0, steps_CB3D+1):
		# variables
		Ekin = 0
		Epot = 0
		sr = au

		# create the path of the celestial body as long at it is between 0 and a ly
		# 	for x, y and z. It also calculates the total energies per step
		for body in bodies_list:
			x, y, z, Ekin_part, Epot_part, dr = \
				body.move3D(bodies_list, dt_CB3D)
			x_list[body.i].append(x/au)
			y_list[body.i].append(y/au)
			z_list[body.i].append(z/au)

			# get the shortest distance between particles for this step
			if dr < sr:
				sr = dr

			# add the energies per body
			Ekin += Ekin_part
			Epot += Epot_part

		# change the timestep according to the distance between particles
		if sr < au:
			step_size = .1
		if sr < .05*au:
			step_size = .01
		if sr < .03*au:
			step_size = .001
		if sr < .01*au:
			step_size = .0001
		else:
			step_size = 1

		# get the time
		t_CB3D += step_size*step*dt_CB3D

		# appending to lists
		t_list.append(t_CB3D/day)
		Ekin_list.append(Ekin)
		Epot_list.append(Epot)
		Etot_list.append(Ekin+Epot)

	return x_list, y_list, z_list, t_list, Ekin_list, Epot_list, Etot_list, \
		bodies_list

# creating the paths and energies for the rutherford scattering using equal
#	charges and opposing charges and lastly for userdefined paramters
x_list_d_eqQ, y_list_d_eqQ, t_list_d_eqQ, Ekin_list_d_eqQ, Epot_list_d_eqQ, \
	Etot_list_d_eqQ, particle_s_eqQ = \
	sd_2particle(0,5,5,1,.5,mprot,mprot,e,e)
x_list_d_oppQ, y_list_d_oppQ, t_list_d_oppQ, Ekin_list_d_oppQ, Epot_list_d_oppQ, \
	Etot_list_d_oppQ, particle_s_oppQ = \
	sd_2particle(0,5,5,1,.5,mprot,mprot,e,-e)
x_list_d_user, y_list_d_user, t_list_d_user, Ekin_list_d_user, Epot_list_d_user, \
	Etot_list_d_user, particle_s_user = \
   sd_2particle(*user_input_list)

# creating the paths and energies for 2 moving particles using opposing charges
x_list_dd, y_list_dd, t_list_dd, Ekin_list_dd, Epot_list_dd, Etot_list_dd = \
	d_particles_2D(0,6,10,4,10,0,-10,0,melec,melec,e,-e)

# creating the paths and energies for a number of 3D particles
x_list3D, y_list3D, z_list3D, t_list3D, Ekin_list3D, Epot_list3D, \
	Etot_list3D, particles3D = moving_particles3D(no_of_particles3D)

# creating the paths and energies for a number of 3D celestial bodies
x_list_CB3D, y_list_CB3D, z_list_CB3D, t_list_CB3D, Ekin_list_CB3D, \
	Epot_list_CB3D, Etot_list_CB3D, bodies_CB3D = \
	celestial_movement3D()

# plot-1 setup
ppl.figure("Path of a dynamic particle in the direction of" + \
	"a charged particle of equal charge")
ppl.title("Path of a dynamic particle in the direction of\n" + \
	"a charged particle of equal charge")
ppl.plot(x_list_d_eqQ, y_list_d_eqQ)
ppl.xlim(0,10)
ppl.ylim(0,10)
ppl.xlabel("x in m")
ppl.ylabel("y in m")
ppl.plot([particle_s_eqQ.x], [particle_s_eqQ.y], "ro", markersize=1)
pp.savefig()

# plot-2 setup
ppl.figure("The total, kinetic and potential energy for a dynamic" + \
	"particle shot towards a charged particle with equal charge")
ppl.title("The total, kinetic and potential energy for a dynamic\n" + \
	"particle shot towards a charged particle with equal charge\n")
ppl.plot(t_list_d_eqQ, Ekin_list_d_eqQ, color="Red")
ppl.plot(t_list_d_eqQ, Epot_list_d_eqQ, color="Green")
ppl.plot(t_list_d_eqQ, Etot_list_d_eqQ, color="Blue")
ppl.legend(["Ekin", "Epot", "Etot"])
ppl.xlabel("time in s")
ppl.ylabel("Energy in J")
ppl.xlim(0, tmax)
pp.savefig()

# plot-3 setup
ppl.figure("Path of a dynamic particle in the direction of" + \
	"a charged particle of opposing charge")
ppl.title("Path of a dynamic particle in the direction of\n" + \
	"a charged particle of opposing charge")
ppl.plot(x_list_d_oppQ, y_list_d_oppQ)
ppl.xlim(0, 10)
ppl.ylim(0, 10)
ppl.xlabel("x in m")
ppl.ylabel("y in m")
ppl.plot([particle_s_oppQ.x], [particle_s_oppQ.y], "ro", markersize=1)
pp.savefig()

# plot-4 setup
ppl.figure("The total, kinetic and potential energy for a dynamic" + \
	"particle shot towards a charged particle with opposing charge")
ppl.title("The total, kinetic and potential energy for a dynamic\n" + \
	"particle shot towards a charged particle with opposing charge\n")
ppl.plot(t_list_d_oppQ, Ekin_list_d_oppQ, color="Red")
ppl.plot(t_list_d_oppQ, Epot_list_d_oppQ, color="Green")
ppl.plot(t_list_d_oppQ, Etot_list_d_oppQ, color="Blue")
ppl.legend(["Ekin", "Epot", "Etot"])
ppl.xlabel("time in s")
ppl.ylabel("Energy in J")
ppl.xlim(0, tmax)
pp.savefig()

# plot-5 setup
ppl.figure("Path of a dynamic particle in the direction of" + \
	"a charged particle with userdefined begin parameters")
ppl.title("Path of a dynamic particle in the direction of\n" + \
	"a charged particle with userdefined begin parameters")
ppl.plot(x_list_d_user, y_list_d_user)
ppl.xlim(0, 10)
ppl.ylim(0, 10)
ppl.xlabel("x in m")
ppl.ylabel("y in m")
ppl.plot([particle_s_user.x], [particle_s_user.y], "ro", markersize=1)
pp.savefig()

# plot-6 setup
ppl.figure("The total, kinetic and potential energy for a dynamic" + \
	"particle shot towards a charged particle, userdifined parameters")
ppl.title("The total, kinetic and potential energy for a dynamic\n" + \
	"particle shot towards a charged particle, userdifined parameters\n")
ppl.plot(t_list_d_user, Ekin_list_d_user, color="Red")
ppl.plot(t_list_d_user, Epot_list_d_user, color="Green")
ppl.plot(t_list_d_user, Etot_list_d_user, color="Blue")
ppl.legend(["Ekin", "Epot", "Etot"])
ppl.xlabel("time in s")
ppl.ylabel("Energy in J")
ppl.xlim(0, tmax)
pp.savefig()

# plot-7 setup
ppl.figure("Path of 2 moving electons with opposing charges")
ppl.title("Path of 2 moving electons with opposing charges")
ppl.plot(x_list_dd[0], y_list_dd[0], color="Red")
ppl.plot(x_list_dd[1], y_list_dd[1], color="Green")
ppl.legend(["Test-particle", "charged-particle"])
ppl.xlim(0, 10)
ppl.ylim(0, 10)
ppl.xlabel("x in m")
ppl.ylabel("y in m")
pp.savefig()

# plot-8 setup
ppl.figure("The kinetic, potential and total energy of " + \
	"2 moving electons with opposing charges")
ppl.title("The kinetic, potential and total energy of\n" + \
	"2 moving electons with opposing charges")
ppl.plot(t_list_dd, Ekin_list_dd, color="Red")
ppl.plot(t_list_dd, Epot_list_dd, color="Green")
ppl.plot(t_list_dd, Etot_list_dd, color="Blue")
ppl.legend(["Ekin", "Epot", "Etot"])
ppl.xlabel("time in s")
ppl.ylabel("Energy in J")
ppl.xlim(0, tmax_dd)
pp.savefig()

# plot-9 setup
fig3D = ppl.figure("The paths of " + str(no_of_particles3D) + \
	" charged particles in 3D")
ax3D = Axes3D(fig3D)
ax3D.set_title("The paths of " + str(no_of_particles3D) + \
	" charged particles in 3D")\

# plot all particles
for particle in particles3D:
	ax3D.plot(x_list3D[particle.i], y_list3D[particle.i], \
		z_list3D[particle.i], linewidth=.5)

ax3D.set_xlim(0, 1)
ax3D.set_ylim(0, 1)
ax3D.set_zlim(0, 1)
ax3D.set_xlabel("x in m")
ax3D.set_ylabel("y in m")
ax3D.set_zlabel("z in m")
pp.savefig()

# rotation function of the 3D plot
def animate3D(aziangle):
	ax3D.view_init(elev=10., azim=aziangle)

	return fig3D

# plot-10 setup
ppl.figure("The kinetic, potential and total energy of " + \
	str(no_of_particles3D) + " moving charged particles")
ppl.title("The kinetic, potential and total energy of \n" + \
	str(no_of_particles3D) + " moving charged particles")
ppl.plot(t_list3D, Ekin_list3D, color="Red")
ppl.plot(t_list3D, Epot_list3D, color="Green")
ppl.plot(t_list3D, Etot_list3D, color="Blue")
ppl.legend(["Ekin", "Epot", "Etot"])
ppl.xlabel("time in s")
ppl.ylabel("Energy in J")
ppl.xlim(0, t_list3D[len(t_list3D)-1])
pp.savefig()

# plot-11 setup
fig_CB3D = ppl.figure("The paths of " + str(len(bodies_CB3D) - 1) + \
	" planets and the sun in 3D for " + str(tmax_CB3D/year) + " years")
ax_CB3D = Axes3D(fig_CB3D)
ax_CB3D.set_title("The paths of " + str(len(bodies_CB3D) - 1) + \
	" planets and the sun in 3D\n for " + str(tmax_CB3D/year) + " years")

# plot all paths of the bodies
for body in bodies_CB3D:
	size = body.R/bodies_CB3D[0].R
	if body.i != 0:
		ax_CB3D.plot(x_list_CB3D[body.i], y_list_CB3D[body.i], \
			z_list_CB3D[body.i], linewidth=.3)
		ax_CB3D.scatter([body.x/au], [body.y/au], \
			[body.z/au], color="gray", s=20*size)
	else:
		ax_CB3D.plot(x_list_CB3D[body.i], y_list_CB3D[body.i], \
			z_list_CB3D[body.i], linewidth=.3)
		ax_CB3D.scatter([body.x/au], [body.y/au], [body.z/au], \
			color="yellow", s=10*size)

ax_CB3D.legend(["Sun", "Mercury", "Venus", "Earth", "Mars"])
ax_CB3D.set_xlim(-2, 2)
ax_CB3D.set_ylim(-2, 2)
ax_CB3D.set_zlim(-2, 2)
ax_CB3D.set_xlabel("x in au")
ax_CB3D.set_ylabel("y in au")
ax_CB3D.set_zlabel("z in au")
pp.savefig()

# rotation function of the 3D plot
def animate_CB3D(aziangle):
	ax_CB3D.view_init(elev=10., azim=aziangle)

	return fig_CB3D

# plot-12 setup
ppl.figure("The kinetic, potential and total energy of " + \
	str(len(bodies_CB3D)) + " moving charged particles" + \
	"for " + str(tmax_CB3D/year) + " years")
ppl.title("The kinetic, potential and total energy of \n" + \
	str(len(bodies_CB3D)) + " moving charged particles" + \
	" for " + str(tmax_CB3D/year) + " years")
ppl.plot(t_list_CB3D, Ekin_list_CB3D, color="Red")
ppl.plot(t_list_CB3D, Epot_list_CB3D, color="Green")
ppl.plot(t_list_CB3D, Etot_list_CB3D, color="Blue")
ppl.legend(["Ekin", "Epot", "Etot"])
ppl.xlabel("time in days")
ppl.ylabel("Energy in J")
ppl.xlim(0, tmax_CB3D/day)
pp.savefig()

# create the animation of the figures of each step
anim3D = animation.FuncAnimation(fig3D, animate3D, frames=360, interval=100)
anim3D.save(PATH + "/figure-3D.gif", dpi=80, writer="imagemagick")
anim_CB3D = animation.FuncAnimation(fig_CB3D, animate_CB3D, frames=360, interval=100)
anim_CB3D.save(PATH + "/figure-3D_CB.gif", dpi=80, writer="imagemagick")

# plot show and pdf close
pp.close()
ppl.show()
