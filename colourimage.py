import requests
from bs4 import BeautifulSoup
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import statistics

n= 90 #beads
nn= n
def KnotCreation2(KnotNumber, n):
    """
    
    Parameters
    ----------
    KnotNumber :Between 2 - 251 for proper knots. knot index from the knotserver. http://www.colab.sfu.ca/KnotPlot/KnotServer/
    n : number of beads in the pollymer
    Returns
    -------
    Knot = 3d array of length n of knot positions, each with a spring length of 1 between the beads
    """
    def GetKnot(KnotNumber):
        url = 'http://www.colab.sfu.ca/KnotPlot/KnotServer/coord/'+ str(KnotNumber) + '.html'
        response = requests.get(url)
        soup = BeautifulSoup(response.text, "html.parser")
        #find = soup.findAll('pre')
        #found_string = str(find)
        
        Full_scrape = soup.findAll('html')
        Full_str = str(Full_scrape)
        
        def find_between( s, first, last ):
            try:
                start = s.index( first ) + len( first )
                end = s.index( last, start )
                return s[start:end]
            except ValueError:
                return "" 
            
        knot_string = find_between( Full_str, "<pre>", "</pre>")
            
        dt_object = np.zeros(3)
        dt = np.result_type(dt_object)
        array = np.fromstring(knot_string, dtype=dt, sep="\n")
        newlength = int(len(array)/3)
        array = np.reshape(array, (newlength, 3))
        return array
    
    return(GetKnot(KnotNumber))



interpolated = 1000001

#file_name = "knot number "+ str(50) + ".txt"
#text_file = open(file_name, "r")
#array = text_file.read()
array = KnotCreation2(2,n)
length = np.shape(array[0:,0])

length1 = int(length[0])
array2=np.zeros((interpolated,3))
array0=np.zeros((length1+1,3))
array0[0:length1,:]=array[:,:]
length1 = int(length[0])+1
array0[length,:]= array[0,:]


#print(array0)
xp = np.linspace(0,1,length1)
fp = array0[0:,0]
array2[0:,0] = np.interp(np.linspace(0,1,interpolated), xp, fp)

fp = array0[0:,1]
array2[0:,1] = np.interp(np.linspace(0,1,interpolated), xp, fp)

fp = array0[0:,2]
array2[0:,2] = np.interp(np.linspace(0,1,interpolated), xp, fp)

skip = int(np.floor(interpolated/(n-1)))
print(skip)
array3 = array2[0::skip]
print(np.shape(array))
print(np.shape(array3))

#array3 = array
#n= length1-1
lengths = np.zeros(n-1)
for i in np.arange(n-1):
    l = (array3[i,0]-array3[i+1,0])**2+(array3[i,1]-array3[i+1,1])**2+(array3[i,2]-array3[i+1,2])**2
    lengths[i] = np.sqrt(l)
    
    

mean = np.mean(lengths)
l = (array3[n-1,0]-array3[0,0])**2+(array3[n-1,1]-array3[0,1])**2+(array3[n-1,2]-array3[0,2])**2
print(np.sqrt(l)/mean)
print(mean,"mean")
first_range = np.ptp(lengths)/mean
print(first_range,"range")

array3 = np.divide(array3,mean)
print(array3)
print(lengths/mean)
"""
m = 0
r_list = []
while m < (len(knot)-1):
    r1 = knot[m]
    r2 = knot[m+1]
    r = r1 - r2
    r_abs = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
    r_list.append(r_abs)
    m = m + 1
mean_r = statistics.mean(r_list)
range_r = np.amax(r_list) - np.amin(r_list)
print("mean of r is ", mean_r)
print("Range of r is ", range_r)
"""


xvals = array3[0:, 0]
yvals = array3[0:, 1]
zvals = array3[0:, 2]
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(xvals, yvals, zvals, 'blue')

import numpy as np
import pylab 
import random 
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from random import gauss
import math
import time
import numba
from numba import jit


k = 1.5 #FENE spring constant /units
rmax= 1 #Lennard Jones 0 potential distance /meters
n = nn#Number of beads /integer
T =500 # Total simulation time seconds
dt =1 #Time step /seconds
setup1 = np.zeros((n,4)) #Empty array for intial bead postion with extra coloumn for time
setup2 = np.zeros((n,4))

P_no = 3 #Torus knot no. around cross section of torus /integer
Q_no = 2 #Torus knot no. around Z axis /integer
bond_angle = 2*math.pi/(n)

forcevalues = np.zeros((T + dt, n, 4))
forcevalues2 = np.zeros((T + dt, n+1, 4))
I = np.identity(3)

#creates empty matrix for mobility matrix
ns = 1 #some stokes constant check what real value is
a = 1 #diamiter or maybe radius must check
#print(V2)
kb = 1
t = 1



setup1[:,0:3] = array3



def fene(r ,r_max, k):
    """
    Computes and returns the source term (right-hand side)
    of the Poisson equation.
    
    Parameters
    ----------
    x : numpy.ndarray
        The gridline locations in the x direction
        as a 1D array of floats.
    y : numpy.ndarray
        The gridline locations in the y direction
        as a 1D array of floats.
    Lx : float
        Domain length in the x direction.
    Ly : float
        Domain length in the y direction.
    
    Returns
    -------
    b : numpy.ndarray of floats
        The forcing function as a 2D array.
    """
    r_abs = np.dot(r,r)
  
    r_abs = np.sqrt(r_abs)
    if r_abs == 1.5:
        print("Error")
    else:
        r_norm = r/r_abs
    r_final = (-r_norm*k*(r_abs-1)/(1-((r_abs-1)/rmax)**2))
    
    return(r_final)

     
print("hi")




def setups():
    """
    Computes and returns the source term (right-hand side)
    of the Poisson equation.
    
    Parameters
    ----------
    x : numpy ndarray
        The gridline locations in the x direction
        as a 1D array of floats.
    y : numpy.ndarray
        The gridline locations in the y direction
        as a 1D array of floats.
    Lx : float
        Domain length in the x direction.
    Ly : float
        Domain length in the y direction.
    
    Returns
    -------
    b : numpy.ndarray of floats
        The forcing function as a 2D array.
    """
    #all data over all time
    forcevalues[0,0:,0:]=setup1
    
forcevalues[0,0:,0:] = setup1



setup2 = setup1   
 
def motion_calculation():
  """
    Computes and returns the source term (right-hand side)
    of the Poisson equation.
    
    Parameters
    ----------
    x : numpy.ndarray
        The gridline locations in the x direction
        as a 1D array of floats.
    y : numpy.ndarray
        The gridline locations in the y direction
        as a 1D array of floats.
    Lx : float
        Domain length in the x direction.
    Ly : float
        Domain length in the y direction.
    
    Returns
    -------
    b : numpy.ndarray of floats
        The forcing function as a 2D array.
  """
  Fij2 = np.zeros(4)
  a = 1
  #setup2 = setup1
  for i in range(dt,T + dt,dt): # dt, T+dt,dt
    setup1 = setup2
    
    
    for j in range(0, n): #random kick
      #print("hi")
      a = 0#np.random.normal(0, 1)/1000
      b = 0#np.random.normal(0, 1)/1000
      c = 0#np.random.normal(0, 1)/1000
      #spring_prev = np.zeros(4) 
      if j==n-1:
         r_prev = setup1[j,0:] - setup1[j-1,0:] #<<This is for the previous bead
         spring_prev = fene(r_prev, rmax, k)/2
      
          #<< this is for the next bead
      #spring_next = np.zeros(4) 
      elif j == 0:
         #<<This is for the previous bead
         spring_prev = np.zeros(4)
      
         r_next = -setup1[j,0:] + setup1[j+1,0:] #<< this is for the next bead
         spring_next = fene(r_next, rmax, k)/2
      else:
         r_prev = setup1[j,0:] - setup1[j-1,0:] #<<This is for the previous bead
         spring_prev = fene(r_prev, rmax, k)/2
      
         r_next = -setup1[j,0:] + setup1[j+1,0:] #<< this is for the next bead
         spring_next = fene(r_next, rmax, k)/2
      
      # for all n here 
      
      Fij =  spring_next[0:4] + spring_prev[0:4] 
      
      setup2[j,0:] = setup1[j,0:] + [a,b,c,0]  + Fij # #add all forces
    setup2[0:,3] = i
    #print(setup2,i,"SETTTUP END")
    forcevalues[i,0:,0:] = setup2[0:,0:] #updating the big array
  return(forcevalues)
  
  
motion_calculation()

#print(forcevalues)
   
def error_checking():
  """
  Computes and returns the source term (right-hand side)
  of the Poisson equation.
    
  Parameters
  ----------
  x : numpy.ndarray
      The gridline locations in the x direction
      as a 1D array of floats.
  y : numpy.ndarray
      he gridline locations in the y direction
      as a 1D array of floats.
  Lx : float
       Domain length in the x direction.
  Ly : float
       Domain length in the y direction.
    
  Returns
  -------
  b : numpy.ndarray of floats
      The forcing function as a 2D array.
  """ 
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
  print("Checked",setup2[j,0])



def correction_for_plot():
    """
    Computes and returns the source term (right-hand side)
    of the Poisson equation.
    
    Parameters
    ----------
    x : numpy.ndarray
        The gridline locations in the x direction
        as a 1D array of floats.
    y : numpy.ndarray
        The gridline locations in the y direction
        as a 1D array of floats.
    Lx : float
        Domain length in the x direction.
    Ly : float
        Domain length in the y direction.
    
    Returns
    -------
    b : numpy.ndarray of floats
        The forcing function as a 2D array.
    """
    
    #print(forcevalues2)
    forcevalues2[:, :n, :] = forcevalues
    #print(forcevalues2)
    forcevalues2[:, n, :] = forcevalues[:, 0, :] #fixes something with matplot lib issue, puerly visual
    #print(forcevalues2,"end")
    #return(forcevalues2)
    
    
correction_for_plot()  

diffsuion = np.zeros((T+dt,4))
def plot_3D(forcevalues):    
    """
    Computes and returns the source term (right-hand side)
    of the Poisson equation.
    
    Parameters
    ----------
    x : numpy.ndarray
        The gridline locations in the x direction
        as a 1D array of floats.
    y : numpy.ndarray
        The gridline locations in the y direction
        as a 1D array of floats.
    Lx : float
        Domain length in the x direction.
    Ly : float
        Domain length in the y direction.
    
    Returns
    -------
    b : numpy.ndarray of floats
        The forcing function as a 2D array.
    """
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    #ax.set_xlim3d([-3, 3.0])
    ax.set_xlabel('X')

    # ax.set_ylim3d([-3, 3.0])
    ax.set_ylabel('Y')

    #ax.set_zlim3d([-1, 1.0])
    ax.set_zlabel('Z')
    for i in range(T,T+dt,dt):
      ax.plot3D(forcevalues[i,0:,0],forcevalues[i,0:,1],forcevalues[i,0:,2])
plot_3D(forcevalues)      

lengths = np.zeros(n-1)
array3 = forcevalues[i,0:,:3]
for i in np.arange(n-1):
    l = (array3[i,0]-array3[i+1,0])**2+(array3[i,1]-array3[i+1,1])**2+(array3[i,2]-array3[i+1,2])**2
    lengths[i] = np.sqrt(l)
    
    

mean = np.mean(lengths)
l = (array3[n-1,0]-array3[0,0])**2+(array3[n-1,1]-array3[0,1])**2+(array3[n-1,2]-array3[0,2])**2
print(np.sqrt(l)/mean)
print(mean,"mean")
print(first_range,"range before")
print(np.ptp(lengths)/mean,"range after")

array10= forcevalues[T,0:,:]
n = np.size(array10[:,0])
print(n)
#print(array10)

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


colours = np.zeros((90,3))
colours[:16,0] =   np.linspace(0,1,16)
colours[15:31,0] =   1- np.linspace(0,1,16)
colours[30:46,1] = np.linspace(0,1,16)
colours[45:61,1] = 1- np.linspace(0,1,16)
colours[60:76,2] = np.linspace(0,1,16)
colours[75:90,2] = 1- np.linspace(0,1,15)
# Make data
shade= 0
def plot(a):
  u = np.linspace(0, 2 * np.pi, 100)
  v = np.linspace(0, np.pi, 100)
  x = 0.5 * np.outer(np.cos(u), np.sin(v))-a[0]
  y = 0.5 * np.outer(np.sin(u), np.sin(v))-a[1]
  z = 0.5 * np.outer(np.ones(np.size(u)), np.cos(v))-a[2]
  # Plot the surface
  ax.plot_surface(x, y, z, color=(shade0, shade1, shade2))
  ax.set_xlim3d([-7, 7.0])
  ax.set_ylim3d([-7, 7.0])
  ax.set_zlim3d([-7, 7.0])
  rc('font',size=28)
  rc('font',family='serif')
  rc('axes',labelsize=26)
  ax.grid(False)
  ax.xaxis.pane.set_edgecolor('black')
  ax.yaxis.pane.set_edgecolor('black')
  ax.zaxis.pane.set_edgecolor('black')
  ax.xaxis.pane.fill = False
  ax.yaxis.pane.fill = False
  ax.zaxis.pane.fill = False
  ax.xaxis.set_ticklabels([])
  ax.yaxis.set_ticklabels([])
  ax.zaxis.set_ticklabels([])
  ax.xaxis._axinfo['tick']['inward_factor'] = 0
  ax.xaxis._axinfo['tick']['outward_factor'] = 0.0
  ax.yaxis._axinfo['tick']['inward_factor'] = 0
  ax.yaxis._axinfo['tick']['outward_factor'] = 0.0
  ax.zaxis._axinfo['tick']['inward_factor'] = 0
  ax.zaxis._axinfo['tick']['outward_factor'] = 0.0
  ax.zaxis._axinfo['tick']['outward_factor'] = 0.0
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  
  print(shade)

  
starter10 = 0
for i in range(n):
    shade0 = colours[i-1,0]
    shade1 = colours[i-1,1]
    shade2 = colours[i-1,2]
    plot(array10[i,:])

       
#ax.plot_surface(x, y, z, color=(0, 0, shade)) 

31
#ax.set_title(r'$3_1$ Knot 90 Bead Polymer')
ax.view_init(elev=45, azim=160)
fig.savefig("colour.png")       
#plt.show()