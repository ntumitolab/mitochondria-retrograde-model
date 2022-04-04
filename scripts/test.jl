using ModelingToolkit
using DifferentialEquations
using GalacticOptim
using Optim

rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
p = [1.0, 100.0]
f = OptimizationFunction(rosenbrock)
prob = GalacticOptim.OptimizationProblem(f, x0, p, lb=[-1.0, -1.0], ub=[1.0, 1.0])
sol = solve(prob, Optim.ParticleSwarm(lower=prob.lb, upper=prob.ub, n_particles=100))
