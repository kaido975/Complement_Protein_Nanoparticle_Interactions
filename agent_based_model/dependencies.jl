# this script installs all required packages
using Pkg
Pkg.update()
if ! in("Agents",keys(Pkg.installed())) Pkg.add(Pkg.PackageSpec(;name="Agents", version="5.6.3")) end
if ! in("Plots",keys(Pkg.installed())) Pkg.add("Plots") end
@warn("Packages installed!")