function critical_spacing(distance, concentration)
    m(t, p) =  p[1]*(abs.((t .- p[3]))./p[3]).^p[2] .* exp.( p[4] * ((t .- p[3])./p[3]) ) 
    p0 = [0.8, 2.3, 33.0, -4.5]
    idx = findall((37 .> distance .> 25.0))
    tdata = distance[idx]
    ydata =  concentration[idx] 
    lower_p = [0.5; -Inf; -Inf; -Inf]
    upper_p = [1.0; Inf; Inf; Inf]
    fit = curve_fit(m, tdata, ydata, p0; lower = lower_p, upper = upper_p)
    params = fit.param
    params = round.(params, digits=1)
    exponent = params[2]
    B = params[4]
    crit = params[3]
    print(params, "\n")
    #Plot and equation
    equation = L"\mathbf{A \|d_r\|^{%$exponent} exp(%$B d_r); d_r=\dfrac{d - %$crit}{%$crit}}"
    y_pred = m(tdata, params)
    Plots.plot(tdata, y_pred, c = "blue", label="fit")
    display(Plots.scatter!(tdata, ydata, xflip=true, c = "red", label="data",legend=:topleft, labelweight=:bold,
    xtickfontsize=12,ytickfontsize=12,yguidefontsize=12,xguidefontsize=12,legendfontsize=12, title="Model Fit: "*equation,titlefontsize=14,
    xlabel="Average IgG-IgG Distance (nm)", ylabel = "C3b/Binding Site"))
    return
end
