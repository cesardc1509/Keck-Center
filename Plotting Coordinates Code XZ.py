# -*- coding: utf-8 -*-
"""
Created on Fri Dec 31 11:40:36 2021

@author: Cesar Diaz
cdiazcarav@miners.utep.edu
 
W.M. Keck Center for 3D Innovation
The University of Texas at El Paso
"""

import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np

dataset = pd.read_excel('name.xlsx')
dataset = pd.DataFrame(dataset)
dataset = np.array(dataset)

x_coords_index = 2 # For indexes starting at 0
z_coords_index = 4

plt.figure(figsize=(11, 10))
plt.scatter(dataset[:, x_coords_index], dataset[:, z_coords_index], s=10)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim([-22, 22])
plt.ylim([-15, 25])
plt.show()