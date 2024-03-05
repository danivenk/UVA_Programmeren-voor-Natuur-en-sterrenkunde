# Dani van Enk, 11823526
# traffic.py in this program a traffic situation of N cars over M section
#   is made. The average speed per car per step and the average traffic-jam
#   speed per car per step is calulated.

# imports
from matplotlib import pyplot as ppl
from matplotlib import animation as ani
import math
import os
import random as rd
from matplotlib.backends.backend_pdf import PdfPages
from fractions import Fraction as F
import imageio as imio

# Path to current file
PATH = os.path.dirname(os.path.abspath(__file__))

# variables
vmax = 10
numeratormax = 101
r = 1
N = 50
M = 100
steps = 100
rad_convert = math.pi/180
section_size = F(str(360/M))
pp = PdfPages(PATH + "/figures-traffic.pdf")

# lists
x_list = [[] for step in range(0, steps)]
y_list = [[] for step in range(0, steps)]
d_list = [[] for step in range(0, steps)]
car_no_list = [car for car in range(0, N)]
density_list = []
av_v_list = []
av_v_tj_list = []

# Car_def is a class to make a car, it has the properties:
#   speed, car number, angle, x&y coords
# It also has a function to calculate the distance to the next car
# Lastely it has a function  to create a different step
class Car_def:
    # __init__() has the initial values of a car
    def __init__(self, i,N):
        rd.seed()
        self.v = int(vmax*rd.random())
        self.i = i
        self.angle = i*section_size*M/N
        self.x = r*math.cos(self.angle*rad_convert)
        self.y = r*math.sin(self.angle*rad_convert)

    # distance() calculates the distance to the next car
    def distance(self,car_in_front, N):
        # checks if the next car is at angle = 0 or not
        if self.i == N-1 and car_in_front.angle == 0 or \
            car_in_front.angle == self.angle:
            car_angle = car_in_front.angle + 360
        else:
            car_angle = car_in_front.angle
        self.d = abs((car_angle - self.angle)/section_size)

        return self.d
    
    # step() sets the next step
    def step(self):
        # if v under vmax accelerate
        if self.v < vmax:
            self.v += 1
        
        # if v is larger than distance to next car reduce speed
        if self.v >= self.d:
            self.v = self.d-1

        # if v is larger than 0 there is a random chance to decelerate
        if self.v > 0:
            rd.seed()
            p = rd.random()

            if p > .5:
                self.v -= 1

        # chance position
        self.angle += self.v*section_size
        self.x = r*math.cos(self.angle*rad_convert)
        self.y = r*math.sin(self.angle*rad_convert)

        return self.v, self.x, self.y

# average_v() calculates the average speed per step per car
def average_v(N, M, vmax):
    # initial variable
    vtot = 0

    # create N cars
    cars = [Car_def(car, N) for car in range(0,N)]

    # Loop for each step and car
    for step in range(0, steps):
        for car in range(0, len(cars)):
            # loop around if index is larger than length
            #   and find distance between cars
            try:
                cars[car].distance(cars[car+1], N)
            except IndexError:
                carI = 0
                cars[car].distance(cars[carI], N)

            # add speed and next step
            vtot += cars[car].v
            cars[car].step()
    
    return vtot/(steps*len(cars))

# average_v() calculates the average speed per step per car in a traffic-jam
def average_v_traffic_jam(N, M, vmax):
    # initial variable
    vtot = 0
    traffic_jam = 0

    # create N cars
    cars = [Car_def(car, N) for car in range(0, N)]

    # Loop for each step and car
    for step in range(0, steps):
        for car in range(0, len(cars)):
            # loop around if index is larger than length
            #   and find distance between cars
            try:
                cars[car].distance(cars[car+1], N)
            except IndexError:
                carI = 0
                cars[car].distance(cars[carI], N)

            # only add speed if there is a traffic jam
            if cars[car].d == 1:
                vtot += cars[car].v
                traffic_jam += 1

            # next step
            cars[car].step()

    # calculate average speed, if no traffic jam the average speed is 0
    try:
        v_av = vtot/traffic_jam
    except ZeroDivisionError:
        v_av = 0

    return v_av

# create an animation from the plots in the Animation/ directory
def animate(files):
    # list
    image = []

    # append the imread of each file
    for file in files:
        image.append(imio.imread(PATH + "/Animation/" + file))

    # create the gif-file
    imio.mimsave(PATH + "/figures-traffic.gif", image)

# for the densities make lists of the average speed overall
#   and in the trafic jam and append the densities
for N_dens in range(1,M):
    av_v_list.append(average_v(N_dens, M, vmax))
    av_v_tj_list.append(average_v_traffic_jam(N_dens, M, vmax))
    density_list.append(N_dens/M)

# create N cars from the variables
cars = [Car_def(car, N) for car in range(0, N)]


# loop for each step and car
for step in range(0, steps):
    for car in range(0, len(cars)):
        # append positions per step for each car
        x_list[step].append(cars[car].x)
        y_list[step].append(cars[car].y)
        # loop around if index is larger than length
        #   and find distance between cars
        try:
            cars[car].distance(cars[car+1], N)
        except IndexError:
            carI = 0
            cars[car].distance(cars[carI], N)

        # next step
        cars[car].step()

    # plot step and save as figure
    ppl.figure("Traffic for " + str(N) + " cars for " + str(steps) + " steps")
    ppl.title("Traffic for " + str(N) + " cars for " + str(steps) + " steps")
    ppl.plot(x_list[step], y_list[step], "ro", markersize=3)
    ppl.xlim(-1.5, 1.5)
    ppl.ylim(-1.5, 1.5)
    ppl.xlabel("x")
    ppl.ylabel("y")
    ppl.savefig(PATH + "/Animation/figure-" + str(step) + ".png")
    ppl.plot(x_list[step], y_list[step], "wo", markersize=4)

# create the animation of the figures of each step
animate(os.listdir(PATH + "/Animation/"))

# plot-2 setup
ppl.figure("Average speed per car per step as function of the car density")
ppl.title("Average speed per car per step as function of the car density")
ppl.semilogx(density_list, av_v_list)
ppl.xlim(1/M, 1)
ppl.ylim(0,vmax)
ppl.xlabel("Car density in N/M")
ppl.ylabel("Average speed per car per step")
pp.savefig()

# plot-3 setup
ppl.figure("Average speed per car per step in a traffic-jam as function of the car density")
ppl.title("Average speed per car per step in a traffic-jam \nas function of the car density")
ppl.semilogx(density_list, av_v_tj_list)
ppl.xlim(1/M, 1)
ppl.ylim(0,vmax)
ppl.xlabel("Car density in N/M")
ppl.ylabel("Average speed per car per step")
pp.savefig()

# plot show and pdf close
pp.close()
ppl.show()
