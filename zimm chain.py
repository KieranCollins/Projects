import scipy
import numpy as np

import random 
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from random import gauss
import math

k = 1 #FENE spring constant /units
rmax= 1 #Lennard Jones 0 potential distance /meters
n = 2#Number of beads /integer
T =21000 # Total simulation time seconds
dt =1 #Time step /seconds
I = np.identity(3)
ns = 1 #some stokes constant check what real value is
a = 1 #diamiter or maybe radius must check
kb = 1
t = 1
ns = 1 #some stokes constant check what real value is
a = 1 #diamiter or maybe radius must check
k = 1
temp= 1

forcevalues = np.zeros((T + dt, n, 4))
forcevalues2 = np.zeros((T + dt, n+1, 4))
setup1 = np.zeros((n,4)) #Empty array for intial bead postion with extra coloumn for time
setup2 = np.zeros((n,4))

def torus_setup(setup1):    
    """
    Computes and returns the x,y,z coordinates of the intial polymer set up.
    
    Parameters
    ----------
    bond_angle : scalar
        Angle between beads.
        
    setup1 : numpy.ndarray
        Empty 2D numpy array of floats.
    Returns
    -------
    setup1 : numpy.ndarray of floats
        The x,y,z coordinates of the intial polymer set up.
    """
    for i in range(0, n):
       setup2[i,0] = i 
       setup2[i,1] = 0 # starting position of beads (chain along x axis)
       setup2[i,2] = 0 
    
    return (setup1)

setup1 = torus_setup(setup2)


def fene(r):
    
    rmax = 1
    r_abs = np.dot(r,r)
      
    r_abs = np.sqrt(r_abs)
    if r_abs == 1.5:
        print("Error")
    else:
        r_norm = r/r_abs
    y = (-r_norm*k*(r_abs-1)/(1-((r_abs-1)/rmax)**2))
    return y

# min to equate minimum and lenth of knot or something
min1 = 1

LJ_const = 0.5    
LJ_shift = 2**(1/6)
epsilon = 0.00000000000000000001
def LJ(a,b):
    
    r=b-a
    rr = np.dot(r,r)
    rr = math.sqrt(rr)
    LJ_shift = 2**(1/6)
    LJforce = (r/rr)*LJ_const*epsilon*(((48/(LJ_shift**12))/((rr)**13))-(((24/(LJ_shift**6))/((rr)**7))))
    #print(LJforce)
    if np.dot(LJforce,LJforce) > 0.1:
      print("LJ!!")
      print(LJforce)
    return(LJforce)

forcevalues[0,0:,0:] = setup1[:,:]


def Dij(kb, t, ns, a, forcevalues, I, n):
    """
    kb is boltzmans, t is temperature, ns is solvent viscosity,
    a is bead radius, r is seperation of beads, I is the identity tensor
    n is number of beads in the system
    
    This is based of eq 19 in Ermak
    """
    a = 1
    d_type = np.dtype((np.float, (3,3)))
    D = np.zeros((n,n), dtype=d_type)
    for i in range (n): #this bit only works for 2 beads as r is gonna be different for each set of beads
        for j in range(n):
            r = forcevalues[i,0:3] - forcevalues[j,0:3]
            if i == j:
                D[i][j] = I*(kb*t)/(6*np.pi*ns*a)
            else:
                r_abs = np.dot(r,r)
                r_abs = np.sqrt(r_abs)
                const = (kb*t)/((8*np.pi*ns)*(r_abs)) #creates constant part before the bracket in eq 19
                first = (np.outer(r,r)/(r_abs*r_abs)+ I) #first half of eq 19
                second = (2*a*a/r_abs*r_abs)*((I/3)-np.outer(r,r)/(r_abs**2)) #2nd part of eq 19
                D[i][j] = const*(first + second)
                
    return D



def motion_calculation():
  a = 1
  for i in range(dt, T+dt, dt): 
    setup1 = setup2
    Dijstep = Dij(kb, t, ns, a, setup1, I, n)
    D_reshape = np.reshape(Dijstep, (3*n*n,3))
    blank= np.zeros((3*n,3*n))
    for w in range(n):
        blank[:][3*w:3*w+3]= np.transpose(D_reshape[:][3*w*n:3*(w+1)*n])
    ##############################
    #print(blank)
    #j = 1####
    #for h in range(0, n):
    #   print(blank[0+3*h:3+3*h,0+3*j:3+3*j]) 
    ##############################  
    cov = 2*blank*dt
    ####print(np.shape(cov),"CHHHHHHHECK")
    mean = np.zeros(3*n)
    ####print(np.shape(mean),"CHHHHHHHECK")
    g = np.random.multivariate_normal(mean, cov)
    g = np.reshape(g,(n,3))
    ####print(g,"GGGGG")
    ####print(np.shape(g),"hi")
    xs=g[:,0]
    ys=g[:,1]
    zs=g[:,2]
    for j in range(0, n):
      if j==n-1:
         r_prev = setup1[j,0:] - setup1[j-1,0:] #<<This is for the previous bead
         spring_prev = fene(r_prev)/2
         spring_next = np.zeros(4)
          #<< this is for the next bead
      #spring_next = np.zeros(4) 
      elif j == 0:
         #<<This is for the previous bead
         spring_prev = np.zeros(4)
      
         r_next = -setup1[j,0:] + setup1[j+1,0:] #<< this is for the next bead
         spring_next = fene(r_next)/2
      else:
         r_prev = setup1[j,0:] - setup1[j-1,0:] #<<This is for the previous bead
         spring_prev = fene(r_prev)/2
      
         r_next = -setup1[j,0:] + setup1[j+1,0:] #<< this is for the next bead
         spring_next = fene(r_next)/2
         
      LJtotal = np.zeros((n,3))
      for h in range(0, n):
          if h!=j:
            LJtotal[h,:] = LJ(setup1[j,0:3],setup1[h,0:3])  # for all n here 
          if h ==j:
            LJtotal[h,:] = np.zeros(3)
              
      ran_kick = np.array([xs[j],ys[j],zs[j]])
      #################################################
      
      ran_kick= ran_kick/10000
      #print(ran_kick,"KICK",j)
      #print(LJtotal[0:3],"lJ",j)
      #print(spring_next[0:3],"spring next",j)
      #print(spring_prev[0:3],"spring prev",j)
      #print(np.sum(spring_next[0:3] + spring_prev[0:3] + LJtotal[0:3]),"forces")
      #D_all = np.ones((3,3))/(10*n)
      Fj = 0
      for h in range(0, n):
          D_all = blank[0+3*h:3+3*h,0+3*j:3+3*j]
          if h == j-1:
              spring_next2 = spring_next
              spring_prev2 = spring_prev
          elif h == j+1:
              spring_next2 = spring_next
              spring_prev2 = spring_prev
          else:
              spring_next2 = np.zeros(3)
              spring_prev2 = np.zeros(3)
          Fj += np.dot(D_all,(spring_next2[0:3] + spring_prev2[0:3] + LJtotal[h,:]))*t/(kb*temp)
      #print(Fj,"Fj")
      setup2[j,0:3] = setup1[j,0:3]+ Fj + ran_kick
    #print("NEXT!!!!!!!!!!!!!!!!!!!!!!!!!!!!")     
    forcevalues[i,0:,0:] = setup2[0:,0:] #updating the big array
  return(forcevalues)
  
motion_calculation()


def error_checking():
 
  for j in range(0, n-1):
    if abs(setup2[j,0] -setup2[j+1,0])>1.5:
       print(j,setup2[j,0],"large x error")
    if abs(setup2[j,1] -setup2[j+1,1])>1.5:
       print(j,setup2[j,1],"large y error")
    if abs(setup2[j,2] -setup2[j+1,2])>1.5:
       print(j,setup2[j,2],"large z error")
    dotter = setup2[j,0:] -setup2[j+1,0:]
    dotter2 = np.dot(dotter,dotter)
    dotter3 = math.sqrt( dotter2 )
    if abs(dotter3)>1.5:
       print(dotter3,setup2[j,0:], "check error")

error_checking()


forcevalues2 =  forcevalues


def plot_3D(forcevalues2):    

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    #ax.set_xlim3d([-3, 3.0])
    ax.set_xlabel('X')

    #ax.set_ylim3d([-3, 3.0])
    ax.set_ylabel('Y')

    #ax.set_zlim3d([-1, 1.0])
    ax.set_zlabel('Z')
    
    
    ##forcevalues2[i,0:,0]+=i/10 difusion hack term
    for i in range(0,T+dt,dt): # moves it in x axis 
        
       
       ax.plot3D(forcevalues2[i,0:,0],forcevalues2[i,0:,1],forcevalues2[i,0:,2])



diffsuion = np.zeros((T+1,4))
for i in range(0,T+dt,1): # moves it in x axis 
       diffsuion[i,0] = np.sum(np.abs(forcevalues[i,0:,0]))
       diffsuion[i,1] = np.sum(np.abs(forcevalues[i,0:,1]))
       diffsuion[i,2] = np.sum(np.abs(forcevalues[i,0:,2]))
       diffsuion[i,3] = i

def diff_inf(TTT,endddd):
    diffusion2 = np.zeros((T-TTT,2))
    for i in range(0,T-TTT,1):
        a = (diffsuion[i,0] - diffsuion[i+TTT,0])
        b = (diffsuion[i,1] - diffsuion[i+TTT,1])
        c = (diffsuion[i,2] - diffsuion[i+TTT,2])
        diffusion2[i,0] = np.sqrt(np.multiply(a,a)+ np.multiply(b,b) + np.multiply(c,c) )
        diffusion2[i,1] = i
        ####mask here with median
    diffusion2 = np.abs(diffusion2)
    median = np.mean(diffusion2[T-endddd:,0])
    std =np.ptp(diffusion2[T-endddd:,0])
    #percentile = np.percentile(diffusion2[T-endddd:,0], 16) 
    #mean = np.mean(diffusion2[T-endddd:,0])
    print(median/TTT,std/TTT,TTT) 
    
countbBb = 0
for b in [1,2,4,8,16,32,64,128,256,512,1024,2048,4096]:#,4096*8]: 
    diff_inf(b,5000)
    
for b in [4096*2,4096*4]:#,4096*8]: 
    diff_inf(b,20000)
    

RADIUS =np.zeros(T+1)
for i in range(dt,T ,dt):
    diffsuion = forcevalues[i,0:,0:]
    radi = 0 
    for j in range(0,n):
        a1 = (diffsuion[0,0] - diffsuion[j,0])
        b1 = (diffsuion[0,1] - diffsuion[j,1])
        c1 = (diffsuion[0,2] - diffsuion[j,2])########################
        radi += np.sqrt((1/(2*n**2)*(np.multiply(a1,a1)+ np.multiply(b1,b1) + np.multiply(c1,c1) )))
    RADIUS[i-1]= radi/(n**2)
#print(np.mean(RADIUS[0:-2]))
#np.std(RADIUS[0:-2])
#print(np.std(RADIUS[0::10]))
print(RADIUS[0::1000])
print(n)


