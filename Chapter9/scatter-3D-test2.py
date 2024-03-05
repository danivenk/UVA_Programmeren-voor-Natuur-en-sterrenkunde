# Dani van Enk, 11823526
# scatter.py calculates the paths and energies of particles using the
# 	using the rutherford scattering idea and using 2 moving particles

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
t = 0
tmax = 25
tmax_dd = 5
tmax_3D = 1
tmax_CB3D = 1000*day
steps = 100000
steps_3D = 10000
steps_CB3D = 1
epsilon_0 = 8.854e-12
e = 1.60218e-19
c = 299792458
G = 6.67408e-11
mprot = 1.672623e-27
melec = 9.10938356e-31
mjup = 1.898e27
rad_convert = math.pi/180
ly = 9.461e15
au = 1.496e11
dt = F(tmax,steps)
dt_dd = F(tmax_dd, steps)
dt_3D = F(tmax_3D, steps_3D)
dt_CB3D = F(tmax_CB3D, steps_CB3D)
correct_input = False
pp = PdfPages(PATH + "/figures-scatter-3D-test2.pdf")

# Get input parameters from user for the rutherford scattering
while correct_input == False:
	user_input_3D = input("How many charged particles would you like" + \
		" to create in 3D?\nPlease give an integer answer\n")

	# Exclude any non number input or in the wrong format
	try:
		no_of_particles3D = int(user_input_3D)

		correct_input = True
	except ValueError:
		print("FormatError: Please try again in the correct format")
		correct_input = False

# the particle class defines a particle with inital position (x_init,y_init),
# 	speed, angle of that speed, mass and charge
#	it also has a function to describe the 2D motion of such charged particle
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

	# move3D() describes the motion of a charged particle in 3D given the
	# 	interacting other particles, the timestep. There is also the
	# 	posibility to have 1 particle be static
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

				if self.i < particle.i:				
					Epot += 1/(4*math.pi*epsilon_0) * \
						(particle.Q*self.Q)/dr

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
	# masses = [mprot, melec]
	particle_list = []
	particles_out_of_bound = []

	# create particles with random inital parameters
	for particle in range(0, particles):
		# random variables
		x = rd.random()
		y = rd.random()
		z = rd.random()
		v = rd.random()
		polangle = 180*rd.random()
		aziangle = 360*rd.random()
		mass = mprot # rd.choice(masses)
		charge = rd.choice(charges)

		particle_list.append(Particle(particle, x, y, z, v, aziangle, polangle, \
			mass, charge))

	# calculate the position and energy for each particles
	for step in range(0, steps_3D):
		# variables
		Ekin = 0
		Epot = 0
		step_size = 1

		# create the path of the particles,
		# 	it also calculates the total energies per step
		for particle in particle_list:
			x, y, z, Ekin_part, Epot_part, dr = \
				particle.move3D(particle_list, step_size*dt_3D)
			particles_out_of_bound.append(particle)
			x_list[particle.i].append(x)
			y_list[particle.i].append(y)
			z_list[particle.i].append(z)

			print(step, step_size, particle.i)

			# add the energies per particle
			Ekin += Ekin_part
			Epot += Epot_part

			if dr < step_size:
				step_size = dr

		# appending to lists
		t_list.append(step*dt_3D)
		Ekin_list.append(Ekin)
		Epot_list.append(Epot)
		Etot_list.append(Ekin+Epot)

		"""
		Sometimes the total energy is not constant this is because if the 
		particles come in really close proximity, the kinetic energy grows bigger
		but the potential energy not large enough to cancel this out. This
		causes the total energy to jump when 2 opposed charged particles 
		come in close proximity
		"""

	return x_list, y_list, z_list, t_list, Ekin_list, Epot_list, Etot_list, \
		particle_list

# creating the paths and energies for a number of 3D particles
x_list3D, y_list3D, z_list3D, t_list3D, Ekin_list3D, Epot_list3D, \
	Etot_list3D, particles3D = moving_particles3D(no_of_particles3D)

# plot-9 setup
fig3D = ppl.figure("The paths of " + str(no_of_particles3D) + \
	" charged particles in 3D")
ax3D = Axes3D(fig3D)
ax3D.set_title("The paths of " + str(no_of_particles3D) + \
	" charged particles in 3D")
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
ppl.xlim(0, tmax_3D)
pp.savefig()

# create the animation of the figures of each step
# anim3D = animation.FuncAnimation(fig3D, animate3D, frames=360, interval=100)
# anim3D.save(PATH + "/figure-3D.gif", dpi=80, writer="imagemagick")

# plot show and pdf close
pp.close()
# ppl.show()
