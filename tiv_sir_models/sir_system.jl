function SIRfunc!(dydt, y, p, t)
    dydt[1] = ( -K_mod * y[1] * y[2] )
    dydt[2] = (K_mod * y[1] * y[2]  - kf_C3bBbC3b * y[2]  )
end
