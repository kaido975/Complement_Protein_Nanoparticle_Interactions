# this script installs all required packages
using Pkg
Pkg.update()
if ! in("DifferentialEquations",keys(Pkg.dependencies())) Pkg.add("DifferentialEquations") end
if ! in("Plots",keys(Pkg.dependencies())) Pkg.add("Plots") end
if ! in("LsqFit",keys(Pkg.dependencies())) Pkg.add("LsqFit") end
if ! in("LaTeXStrings",keys(Pkg.dependencies())) Pkg.add("LaTeXStrings") end
@warn("Packages installed!")