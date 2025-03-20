function tiv_call(conc_c3b_sites)
    tend = 10*60 #seconds
    tspan = (0.0, tend)
    y0 = zeros(3)
    y0[1] = conc_c3b_sites*1e6
    y0[2] = 1e-5 
    prob = ODEProblem(TIVfunc!, y0, tspan)
    sol = solve(prob, AutoTsit5(Rosenbrock23()))
    return sol
end
