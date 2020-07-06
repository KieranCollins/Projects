
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

for KnotNumber in range(251, 252): # 2,252
    print(KnotNumber)
    file_name = "knot number "+ str(KnotNumber) + ".txt"
    #knot_string = str(knot)
    text_file = open(file_name, "r")
    array = text_file.read()
    text_file.close()
    #print(array)
    for n in [10,20,30,40]:
        print(array)