using DifferentialEquations
using Plots
using LsqFit
using LaTeXStrings

include("parameters.jl")
include("sir_system.jl")
include("protein_spacing.jl")
include("sir_call.jl")
include("tiv_call.jl")
include("tiv_system.jl")
include("criticality.jl")
#Number of binding sites / particle
nps = collect(range(5, 250, 250-4))

values_tiv = [] #Store relative C3b concentration
values_sir = [] #Store relative C3b concentration

for sites_np in nps
    conc_c3b_sites = sites_np*conc_np
    #Normalized C3b Concentration TIV model
    sol_tiv = tiv_call(conc_c3b_sites)
    push!(values_tiv, maximum(sol_tiv[2, :])/ sol_tiv[1, 1])

    #Normalized C3b concentration from SIR model    
    sol_sir = sir_call(conc_c3b_sites)
    push!(values_sir, maximum(sol_sir[2, :])/ sol_sir[1, 1])
end

Plots.scatter(distance, values_tiv, xflip=true, label=false, c = "blue")
display(Plots.plot!(distance, values_tiv, xlabel="Protein-Protein spacing (nm)", ylabel="C3b (μM)/Binding Site (μM)", title="TIV model",
                    xtickfontsize=14,ytickfontsize=14,yguidefontsize=14,xguidefontsize=14,legendfontsize=14, xflip=true, label=false, c = "blue"))

Plots.scatter(distance, values_sir, xflip=true, label=false, c = "blue")
display(Plots.plot!(distance, values_sir, xlabel="Protein-Protein spacing (nm)", ylabel="C3b (μM)",title="SIR model",
                    xtickfontsize=14,ytickfontsize=14,yguidefontsize=14,xguidefontsize=14,legendfontsize=14, xflip=true, label=false, c = "blue"))


#Estimating critical spacing for TIV model
critical_spacing(distance, values_tiv)

#Estimating critical spacing for SIR model
critical_spacing(distance, values_sir)
