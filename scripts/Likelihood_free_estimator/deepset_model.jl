# deepset_model.jl
using Flux

d = 2    # dimension of each replicate
w = 32   # number of neurons in each hidden layer

# Final layer with clamp instead of sigmoid
final_layer = Dense(w, 1, x -> clamp.(x, 0.0, 1.0))

psi = Chain(Dense(d, w, relu), Dense(w, w, relu), Dense(w, w, relu))
phi = Chain(Dense(w, w, relu), Dense(w, w, relu), final_layer)

deepset = DeepSet(psi, phi)
estimator = PointEstimator(deepset)
