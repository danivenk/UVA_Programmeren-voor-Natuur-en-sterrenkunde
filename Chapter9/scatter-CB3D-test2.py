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
year = int(365.25*day)
t = 0
tmax = 25
tmax_dd = 5
tmax_3D = 1
tmax_CB3D = 3*year
steps = 100000
steps_3D = 10000
steps_CB3D = 10000
epsilon_0 = 8.854e-12
e = 1.60218e-19
c = 299792458
G = 6.67408e-11
mprot = 1.672623e-27
melec = 9.10938356e-31
mjup = 1.898e27
Rjup = 69911000
rad_convert = math.pi/180
ly = 9.461e15
au = 1.496e11
dt = F(tmax,steps)
dt_dd = F(tmax_dd, steps)
dt_3D = F(tmax_3D, steps_3D)
dt_CB3D = F(tmax_CB3D, steps_CB3D)
no_of_particles3D = 10
correct_input = False
pp = PdfPages(PATH + "/figures-scatter-CB3D-test2.pdf")

class Celestial_body:
	# __init() contains the begin parameters, those parameters are
	# 	number of particles, position (x_init,y_init,z_init), speed,
	# 	angle of the speed, mass and charge
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

	# use the Adams-Bash method to calculate the next step
	def Adams_Bash(self, number,fi,fi_1,h):
		return number + h*(3/2*fi-1/2*fi_1)

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
				a = math.sqrt(ax**2 + ay**2 + az**2)

				if self.i < body.i:
					Epot += -(G*self.m*body.m)/dr

				if dr < smallestr:
					smallestr = dr
					# print(dr/au)

		return ax, ay, az, Ekin, Epot, smallestr

	# move3D() describes the motion of a charged body in 3D given the
	# 	interacting other bodies, the timestep. There is also the
	# 	posibility to have 1 body be static
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

# creating the paths and energies for a number of 3D celestial bodies
x_list_CB3D, y_list_CB3D, z_list_CB3D, t_list_CB3D, Ekin_list_CB3D, \
	Epot_list_CB3D, Etot_list_CB3D, bodies_CB3D = \
	celestial_movement3D()

# plot-11 setup
fig_CB3D = ppl.figure("The paths of " + str(len(bodies_CB3D) - 1) + \
	" planets and the sun in 3D for " + str(tmax_CB3D/year) + " years")
ax_CB3D = Axes3D(fig_CB3D)
ax_CB3D.set_title("The paths of " + str(len(bodies_CB3D) - 1) + \
	" planets and the sun in 3D\n for " + str(tmax_CB3D/year) + " years")

for body in bodies_CB3D:
	size = body.R/bodies_CB3D[0].R
	if body.i != 0:
		ax_CB3D.plot(x_list_CB3D[body.i], y_list_CB3D[body.i], \
			z_list_CB3D[body.i], linewidth=.3)
		ax_CB3D.scatter([body.x/au], [body.y/au], [body.z/au], color="gray", s=20*size)
	else:
		ax_CB3D.plot(x_list_CB3D[body.i], y_list_CB3D[body.i], \
			z_list_CB3D[body.i], linewidth=.3)
		ax_CB3D.scatter([body.x/au], [body.y/au], [body.z/au], color="yellow", s=10*size)

ax_CB3D.legend(["Mercury", "Venus", "Earth", "Mars"])
ax_CB3D.set_xlim(-2, 2)
ax_CB3D.set_ylim(-2, 2)
ax_CB3D.set_zlim(-2, 2)
ax_CB3D.set_xlabel("x in au")
ax_CB3D.set_ylabel("y in au")
ax_CB3D.set_zlabel("z in au")
pp.savefig()

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
ppl.xlim(0, tmax_CB3D/day)
pp.savefig()

# create the animation of the figures of each step
# anim_CB3D = animation.FuncAnimation(fig_CB3D, animate_CB3D, frames=360, interval=100)
# anim_CB3D.save(PATH + "/figure-3D_CB.gif", dpi=80, writer="imagemagick")

# plot show and pdf close
pp.close()
# ppl.show()
