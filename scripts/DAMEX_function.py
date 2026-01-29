# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 14:10:36 2026

@author: Anas Mourahib 

Function for the DAMEX algorithm 
"""
import sys
sys.path.insert(1, 'C:/Users/20254817/Desktop/Githib/Penalized_least-squares_estimator/extern/DAMEX')

from extreme_data import rank_transformation
from damex_algo import damex


 
import numpy as np

def DAMEX_test(data : np.ndarray  , mu : np.array ,  epsilon = 0.01  ) -> list:
    kn = round( data.shape[0]**(1/2) )
    #kn = round(data.shape[0] * 0.08)
    x_rank_0 = rank_transformation(data)
    ind_extr_0 = np.argsort(np.sum(x_rank_0, axis=1))[::-1][:kn]
    x_extr_0 = x_rank_0[ind_extr_0]
    R = np.min(np.max(x_extr_0 , axis=1 ) )-1
    list_extreme_direction = []
    for threshold in mu:
        extreme_directions, mass_danube = damex(x_extr_0, R, epsilon, mu_min=threshold)
        list_extreme_direction.append(extreme_directions)
    for i, sublist in enumerate(list_extreme_direction):
        for j, subsublist in enumerate(sublist):
            list_extreme_direction[i][j] = [int(x) for x in subsublist]    
    return(list_extreme_direction)


##################This is a  test for a simulation from a mixture logistic model 

#import pyreadr

#X = pyreadr.read_r('C:/Users/20254817/Desktop/Githib/Penalized_least-squares_estimator/Data/simu.Rds')[None]
#X = np.array(X)
#a = 0.01     # start
#b = 0.1    # end (exclusive)
#gap = 0.01     # step size
#mu = np.arange(a, b, gap)


#extr_dir = DAMEX_test(data = X , mu = mu)


    