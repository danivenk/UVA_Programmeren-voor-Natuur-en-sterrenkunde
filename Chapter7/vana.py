import math
import matplotlib.pyplot as plt

class constants:
    def __init__(self,nstep, a,b,V0,nmax):
        self.nstep =nstep
        self.a = a
        self.b = b
        self.V0 = V0
        self.nmax = nmax

        
myconstants =  constants(100,1.,1.,10.,101)

def F1(n):
    return (1./n)*(1-math.pow(-1.,n))

def F2(n,x):
    return math.sin((math.pi*n*x)/myconstants.a)

def F3(n,y):
    return ((math.exp(((math.pi*n)/myconstants.a)*(myconstants.b -y)) -
             math.exp (((-math.pi*n)/myconstants.a)*(myconstants.b-y)))/2.)

def F4(n):
    return (math.exp((math.pi*n*myconstants.b)/myconstants.a)-
            math.exp((-math.pi*n*myconstants.b)/myconstants.a))/2.

def Va(x,y):
    V=0
    Vcontribution =0
    for n in range(1,myconstants.nmax,2):
        Vcontribution = (2.*myconstants.V0/math.pi)*F1(n)*F2(n,x)*F3(n,y)/F4(n)
        V += Vcontribution
        if math.fabs(Vcontribution/V)< math.pow(10.,-20):
            break
    return V

xmin = 0.
xmax = 1.
ymin = 0.
ymax = 1.

hx = (xmax-xmin)/myconstants.nstep
hy = (ymax-ymin)/myconstants.nstep

V = [[0]*myconstants.nstep for i in range (myconstants.nstep)]

x =[]
y =[]

for i in range(0,myconstants.nstep):
    x.append(i*hx)
    y.append(i*hy)

for i in range(1,myconstants.nstep):
    for j in range(1,myconstants.nstep):
        V[i][j]=Va(x[i],y[j])
    V[i][0]=myconstants.V0

Vxy1 = [0]*myconstants.nstep
xy1 = [0]*myconstants.nstep

for i in range(0,myconstants.nstep):
    Vxy1[i]=V[i][int(myconstants.nstep/2)]
    xy1[i]=x[i]

Vxy2 = [0]*myconstants.nstep
xy2 = [0]*myconstants.nstep

for i in range(0,myconstants.nstep):
    Vxy2[i]=V[int(myconstants.nstep/2)][i]
    xy2[i]=y[i]

plt.plot(xy1,Vxy1)
plt.plot(xy2,Vxy2)
plt.ylabel('Potential [a.u.]')
plt.xlabel('position [a.u.]')
# plt.show()
