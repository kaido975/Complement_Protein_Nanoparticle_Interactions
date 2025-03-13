function critical_spacing(distance, concentration)
    m(t, p) =  p[1]*(abs.((t .- p[3]))./p[3]).^p[2] .* exp.( p[4] * ((t .- p[3])./p[3]) ) 
    p0 = [0.5, 2.0, 40.0, 2.0]
    idx = findall((35.0 .> distance .> 25.0))
    tdata = distance[idx]
    ydata =  concentration[idx] 
    fit = curve_fit(m, tdata, ydata, p0)
    params = fit.param
    params = round.(params, digits=1)
    exponent = params[2]
    B = params[4]
    crit = params[3]

    #Plot and equation
    equation = L"\mathbf{A \|d_r\|^{%$exponent} exp(%$B d_r); d_r=\dfrac{d - %$crit}{%$crit}}"
    y_pred = m(tdata, params)
    Plots.plot(tdata, y_pred, c = "blue", label="fit")
    display(Plots.scatter!(tdata, ydata, xflip=true, c = "red", label="data",legend=:topleft, labelweight=:bold,
    xtickfontsize=12,ytickfontsize=12,yguidefontsize=12,xguidefontsize=12,legendfontsize=12, title="Model Fit: "*equation,titlefontsize=14,
    xlabel="Average IgG-IgG Distance (nm)", ylabel = "C3a (μM)"))
    return
    
end