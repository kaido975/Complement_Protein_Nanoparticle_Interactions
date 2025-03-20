function TIVfunc!(dydt, y, p, t)
    dydt[1] = (-k_tick * y[1] - kf_pC3b * y[1] * y[3] )
    dydt[2] = (k_tick * y[1] + kf_pC3b*y[1]*y[3] - kf_C3bBbC3b*y[2])
    dydt[3] = (k_cat * y[2] - div_fac*kf_C3bBbC3b*y[3]) 
end
