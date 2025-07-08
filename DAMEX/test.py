

# Import some files and libraries 
import sys
sys.path.insert(1, 'C:/Users/mourahib/Desktop/github/Estimation-mixture-model/DAMEX')

import Functions
import extreme_data 
import pandas as pd
import numpy as np
import sys
from extreme_data import rank_transformation
from Functions import damex



df_danube = pd.read_csv("C:/Users/mourahib/Desktop/github/Estimation-mixture-model/DAMEX/Danube4.csv")
print(df_danube)




x_rank_0_danube = rank_transformation(df_danube)
#x_rank_0 is a numpy array where each value is itself an array
# kth extremer points for the sum-norm
k_0 = int(2e3)
ind_extr_0_danube = np.argsort(np.sum(x_rank_0_danube, axis=1))[::-1][:k_0]
x_extr_0_danube = x_rank_0_danube[ind_extr_0_danube]
R_danube=np.min(np.max(x_extr_0_danube , axis=1 ) )-1
eps_dmx = 0.0005
mu_min = 0.0001
alphas_0_danube, mass_danube = damex(x_extr_0_danube, R_danube, eps_dmx, mu_min=mu_min)
print("With threshold", mu_min, "the extreme directions resulting from stations {1 , 2 , 11 , 12}", alphas_0_danube , "with mass" , mass_danube)

