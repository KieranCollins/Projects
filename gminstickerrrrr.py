import requests
from bs4 import BeautifulSoup
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import statistics
blank = np.zeros((250,2))
def GetKnot(KnotNumber):
    """
    Parameters
    ----------
    KnotNumber :Between 2 - 251 for proper knots. knot index from the knotserver. http://www.colab.sfu.ca/KnotPlot/KnotServer/
    n : number of beads in the pollymer
    Returns
    -------
    Knot = 3d array of length n of knot positions, each with a spring length of 1 between the beads
    """
    url = 'http://www.colab.sfu.ca/KnotPlot/KnotServer/mseq-coord/'+ str(KnotNumber) + '.html'
    response = requests.get(url)
    soup = BeautifulSoup(response.text, "html.parser")   
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
    file_name = "knot number "+ str(KnotNumber) + ".txt"
    array = np.size(array[:,0])
    return array




for m in range(2, 252):#2
    print(m)
    print(GetKnot(m))
    blank[m-2,0] = GetKnot(m)
    blank[m-2,1] = m
    
fname = "Minstick.txt"
np.savetxt(fname, blank, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)    









