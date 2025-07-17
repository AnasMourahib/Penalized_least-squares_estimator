import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

#Import the necessary packages
mev = importr('mev')

# Activate automatic conversion globally (optional)
numpy2ri.activate()

# Load the R script
ro.r.source("C:/Users/mourahib/Desktop/github/PenalizedLeastSquaresEstimator/R/Simulation_mixture_model.R")

# Get the R function
simulate_mixture_model = ro.globalenv['N_generate_Mix_log']

# Prepare the inputs
N = 100
A = np.array([
    [1/3 , 1/3 , 1/3 , 0 ],
    [1/2 , 0   , 0   , 1/2],
    [0   , 3/4 , 1/4 , 0],
    [0   , 1/2 , 0   , 1/2]
])
alpha = 0.25

# Ensure NumPy-to-R conversion works
with localconverter(ro.default_converter + numpy2ri.converter):
    result = simulate_mixture_model(N, A, alpha)

# Print the result
print(result)
