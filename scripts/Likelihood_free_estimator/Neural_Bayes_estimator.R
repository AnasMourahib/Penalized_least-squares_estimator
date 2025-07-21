# Install packages if not already installed
install.packages("JuliaCall")
install.packages("JuliaConnectoR")

# Load libraries
library(JuliaConnectoR)  # optional if you use JuliaCall primarily
library(JuliaCall)

# Setup Julia with the folder containing julia.exe
julia_setup(JULIA_HOME = "C:/Users/mourahib/AppData/Local/Programs/Julia-1.11.6/bin")

# Install packages in Julia if not already installed
julia_command('using Pkg; Pkg.add("NeuralEstimators")')
julia_command('using Pkg; Pkg.add("Flux")')

# Load the packages in the current Julia session
julia_library("NeuralEstimators")
julia_library("Flux")

# Alternatively, you can combine loading packages in one command
# juliaEval("using NeuralEstimators, Flux")



#> Starting Julia ...
#>
