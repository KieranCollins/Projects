import requests
from bs4 import BeautifulSoup
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import statistics

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
    url = 'http://www.colab.sfu.ca/KnotPlot/KnotServer/mseq-coord/+ str(KnotNumber) + '.html'
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
    array = np.shape(array[0:])
    return array,file_name




for m in range(2, 252):
    print(m)
    knot,file_name = GetKnot(m)
    knot_string = str(knot)
    text_file = open(file_name, "w")
    text_file.write(str(knot))
    text_file.close()