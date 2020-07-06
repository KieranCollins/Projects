#import requests
#from bs4 import BeautifulSoup
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import requests
from bs4 import BeautifulSoup

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
for KnotNumber in range(2, 252): # 2,252
    for n in [20,25,30,35,40,45,50]: 
        array = KnotCreation2(KnotNumber, n)
        print(array)
        interpolated = 1000001
        length = np.size(array[0:,0])
        length1 = length
        array2=np.zeros((interpolated,3))
        array0=np.zeros((length1+1,3))
        array0[0:length1,:]=array[:,:]
        length1 = length+1
        array0[length,:]= array[0,:]
        xp = np.linspace(0,1,length1)
        fp = array0[0:,0]
        array2[0:,0] = np.interp(np.linspace(0,1,interpolated), xp, fp)
        fp = array0[0:,1]
        array2[0:,1] = np.interp(np.linspace(0,1,interpolated), xp, fp)
        fp = array0[0:,2]
        array2[0:,2] = np.interp(np.linspace(0,1,interpolated), xp, fp)
        skip = int(np.floor(interpolated/(n-1)))
        array3 = array2[0::skip]
        lengths = np.zeros(n-1)
        for i in np.arange(n-1):
            l = (array3[i,0]-array3[i+1,0])**2+(array3[i,1]-array3[i+1,1])**2+(array3[i,2]-array3[i+1,2])**2
            lengths[i] = np.sqrt(l)
        mean = np.mean(lengths)
        l = (array3[n-1,0]-array3[0,0])**2+(array3[n-1,1]-array3[0,1])**2+(array3[n-1,2]-array3[0,2])**2
        #print(np.sqrt(l)/mean)
        #print(mean,"mean")
        first_range = np.ptp(lengths)/mean
        #print(first_range,"range")
        array3 = np.divide(array3,mean)
        #print(array3)
        #print(lengths/mean)
        import random 
        from random import gauss
        import math
        k = 1.5 #FENE spring constant /units
        rmax= 1 #Lennard Jones 0 potential distance /meters
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
        forcevalues[0,0:,0:] = setup1
        setup2 = setup1     
        motion_calculation()
        #diffsuion = np.zeros((T+dt,4))
        lengths = np.zeros(n-1)
        array3 = forcevalues[i,0:,:3]
        for i in np.arange(n-1):
            l = (array3[i,0]-array3[i+1,0])**2+(array3[i,1]-array3[i+1,1])**2+(array3[i,2]-array3[i+1,2])**2
            lengths[i] = np.sqrt(l)
        mean = np.mean(lengths)
        l = (array3[n-1,0]-array3[0,0])**2+(array3[n-1,1]-array3[0,1])**2+(array3[n-1,2]-array3[0,2])**2
        print("###############################")
        print("n=",n,"knot=",KnotNumber)
        #print(np.sqrt(l)/mean)
        print(mean,"mean")
        print(first_range,"range before")
        rangerrrr = np.ptp(lengths)/mean
        print(rangerrrr,"range after")
        print("###############################")
        array = forcevalues[-1,0:,:3]
        fname = "knot number"+ str(KnotNumber) + "b" +str(n) + ".txt"
        np.savetxt(fname, array, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
        #plot_3D(forcevalues)