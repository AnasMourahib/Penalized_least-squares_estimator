import numpy as np 
import pandas as pd
import sys 
import sys
sys.path.append('./likelihood_free_estimator')
from simulate_model import simulate_model

##Number of saples 

num_samples = 200

N = 1000 

A = np.array([[1/3 , 1/3 , 1/3 , 0],
              [1/2 , 0 ,  0 , 1/2] , 
              [0, 3/4, 1/4, 0],
              [0, 1/2, 0, 1/2] 
              ])

alphas = np.random.uniform(0.0 , 1.0 , size= num_samples)

for i in range(num_samples):
    print(f"i = {i}")
    print(f"alpha = {alphas[i]}")


simulations = []
parameters = []

for i in range(num_samples):
    print("i=", i)
    print(f"alpha = {alphas[i]}")
    result = simulate_model(N, A, alphas[i])    
    result_np = np.array(result)
    print(result_np)
    simulations.append(result_np)
    parameters.append([alphas[i]])

    