# DDM
Estimating Parameters of the Drift Diffusion Model Using Particle Swarm Optimization

Firstly, the goal is to fit a DDM to reaction time and accuracy data. Afterwards, the likelihood of the observed data, given the model parameters, is maximized. This process is equivalent to the minimization of the negative log-likelihood of the observed data, given the model parameter.


The following parameters of the DDM are estimated:

    v = drift rate (speed of evidence accumulation)
    s_v = variability in drift rate across trials
    a = decision boundary separation (threshold)
    s_z = variability in starting point (expressed relative to $a$)
    t_0 = non-decision time (encoding + motor)
    s_t0 = variability in non-decision time
    z = starting point bias (relative to $a$)


In one folder (Output_Function) the objective function to be minimized, vectors for reaction times, accuracies, and condition levels are specified. Further, a vector for the parameters per condition has to be specified, serving as initial parameters.


In another folder the functions to estimate the DDM parameters based on the negative Log-Likelihood are stored (negLL_DDM_final).

For the PSO algorithm (PSO_algorithm_simple), first the number of instruction conditions and DDM parameters have to be specified. Secondly, the total number of parameters (m) is defined. Assuming two experimental conditions, this results in 14 parameters which are estimated. 

Afterwards, the parameter bounds have to be specified. In our example, for the drift rate v the boundary is set between -5 and 5, for its variability s_v the boundaries are 0 and 2. For boundary separation a, boundaries are set between 0.3 and 3. For variability of starting point s_z the boundaries are between 0 and 0.5. For non-decision time t0 they are between 0.1 and 1 and for variability of non-decision time s_t0 between 0 and 0.3. For starting point z the values can only range between 0 and 1. 

Then, number of particles, inertia weights, cognitive and social acceleration coefficients, the maximum number of iterations as well as the maximum number of PSO runs have to be specified. 

     n = 100 (number of particles)
     w_max = 0.9 (inertia weight)
     w_min = 0.4 (inertia weight)
     k1 = 2 (acceleration factor)
     k2 = 2 (acceleration factor)
     max_iteration = 1000
     max_run = 1

This form could be adapted to include anarchy (PSO_algorithm_anarchy, PSO_algorithm_adapted anarchy), rebounding mechanisms (PSO_algorithm_rebounded; PSO_algorithm_anarchy_rebounded, PSO_algorithm_adapted_anarchy_rebounded) or updated weights (PSO_algorithm_simple_weights).
