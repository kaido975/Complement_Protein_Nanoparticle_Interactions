include("angle_compute.jl")
n_sites = collect(range(5, 250, 250-4))
# n_sites = [5, 6, 10, 20, 22, 25, 30, 35,40, 50, 60, 70, 100, 120, 150, 200, 250]
n_sites = Int.(n_sites)
distance = []

r = 0.5 #Radius of the nanopartcile (75 nanometers)
cx = r
cy = r
cz = r
# 1 unit = 150 nanometers

for k = 1:size(n_sites, 1)
    i = range(0, n_sites[k]-1, length = n_sites[k]) .+ 0.5
    ϕ = acos.(1 .- 2*i/n_sites[k])
    ω = (1 + 5^0.5)/2   
    θ = 2*pi * i / ω
    local x, y, z = r*cos.(θ) .* sin.(ϕ) .+ cx, r*sin.(θ) .* sin.(ϕ) .+ cy, r*cos.(ϕ) .+ cz

    d = 100

    for i =1:n_sites[k]
        for j=1:n_sites[k]
            p1 = (x[i], y[i], z[i])
            p2 = (x[j], y[j], z[j])
            d_new = r*Main.AngleBetweenVectors.angle(p1 .- (0.5, 0.5, 0.5), p2 .- (0.5, 0.5, 0.5))
            if d_new < d && d_new!=0
                d = d_new
            end
        end
    end
    push!(distance, d*150) #Protein-Protein spacing in nanometers
end
