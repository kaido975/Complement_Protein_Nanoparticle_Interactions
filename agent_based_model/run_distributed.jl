@everywhere using Agents, Random, Agents.Pathfinding
@everywhere using Distributions
@everywhere using LinearAlgebra
@everywhere using Distributed
@everywhere using Plots, StatsPlots
@everywhere using Statistics
@everywhere ENV["GKSwstype"]="nul"


@everywhere include("PATH/angle_compute.jl") ## Include angle_compute.jl 
cd("PATH") ## Change the path to the folder containing run_distributed.jl

@everywhere n_agents = 200
@everywhere @agent complement_agent ContinuousAgent{3} begin
    type::Symbol          
    engaged::Int  
    state::Int
    age::Float32
end
@everywhere const v0 = (0.0, 0.0, 0.0)
@everywhere c3(id, pos, engaged, state, age) = complement_agent(id, pos, v0, :c3, engaged, state, age)
@everywhere lp(id, pos, engaged, state, age) = complement_agent(id, pos, v0, :lp, engaged, state, age)

@everywhere x = collect(range(0, 1, length = 60))
@everywhere y = collect(range(0, 1, length = 60))
@everywhere z = collect(range(0, 1, length = 60))
@everywhere z = reshape(z, (1, 1, size(z, 1)))
@everywhere y = reshape(y, (1, size(y, 1)))
@everywhere cx = 0.5
@everywhere cy = 0.5
@everywhere cz = 0.5
@everywhere r = 0.5

#Boundary layer for diffusion
@everywhere ϵ = 0.02
@everywhere mask = r^2 + ϵ/10 .<= (x .- cx).^2 .+ (y .- cy).^2 .+ (z .- cz).^2 .<= r^2 + ϵ

@everywhere n = 50
@everywhere u = range(-π, π; length = n)
@everywhere v = range(0, π; length = n)
@everywhere x = r*cos.(u) * sin.(v)'
@everywhere y = r*sin.(u) * sin.(v)'
@everywhere z = r*ones(n) * cos.(v)'

#Distribution of agents on the grid
@everywhere avg_step = 20
@everywhere i = range(0, n_agents-1, length = n_agents) .+ 0.5
@everywhere ϕ = acos.(1 .- 2*i/n_agents)
@everywhere ω = (1 + 5^0.5)/2   
@everywhere θ = 2*pi * i / ω
@everywhere x_lp, y_lp, z_lp = r*cos.(θ) .* sin.(ϕ) .+ cx, r*sin.(θ) .* sin.(ϕ) .+ cy, r*cos.(ϕ) .+ cz


@everywhere cus_scheduler(model) = (a.id for a in allagents(model) if a.state==0)

@everywhere function initialize_model(;
    seed = 42,  ## seed for random number generator,
    walk = mask, ## Boundary layer for agent movement
    surface_diffusion = 0.35, #Surface Diffusion Rate
    dt = 0.1, #Time step
    k_tick = 4e-4, #Tickover rate
    k_cat = 8e-2, #Autocatalysis rate
    decay_time = 5.0 #Bound C3b time to decay
)
   
    lps = Int.(collect(range(1, n_agents, length = n_agents)))  #Binding Sites

    rng = MersenneTwister(seed)
    space = ContinuousSpace((1, 1, 1); periodic = false)

    properties = (
        walkable_path = AStar(space; walkmap = mask),
        surface_diffusion = surface_diffusion,
        dt = dt,
        k_tick  =k_tick,
        k_cat = k_cat, 
        decay_time = decay_time,
        lps = lps
    )

    model = ABM(complement_agent, space; rng, properties)
    for i in 1:n_agents
        id = nextid(model)
        agent = lp(id, (x_lp[i], y_lp[i], z_lp[i]), 0, 1, 0)
        add_agent_pos!(agent, model)  
    end

    return model
end

@everywhere function agent_step!(agent::complement_agent, model)
    if agent.type == :c3
        c3_step!(agent, model)
    else
        lp_step!(agent, model)
    end
end

@everywhere function c3_step!(c3_agent, model)

    c3_agent.age += model.dt
    if c3_agent.engaged != 0
        if c3_agent.age > model.decay_time
            c3_agent.state = 1
            return 
        end
    end
    
    if c3_agent.engaged == 0

        if c3_agent.age > model.decay_time/3
            kill_agent!(c3_agent, model, model.walkable_path)
            return
        end  

        if c3_agent.state != 1 
            lp_neigh = [x for x in nearby_agents(c3_agent, model, 0.0005) if x.type == :lp]  
            if !isempty(lp_neigh)
                for lp_agent in lp_neigh  
                    if lp_agent.engaged == 0
                        c3_agent.engaged = 1
                        lp_agent.engaged = 1
                        move_agent!(c3_agent, lp_agent.pos, model)
                        deleteat!(model.lps, findall(x->x==lp_agent.id, model.lps))
                        return
                    end
                end
            end
        end
        
        if is_stationary(c3_agent, model.walkable_path)
            plan_route!(c3_agent, random_walkable(model, model.walkable_path), model.walkable_path)     
        end
        move_along_route!(c3_agent, model, model.walkable_path, model.surface_diffusion, model.dt)        
        return  
    else 
        if length(model.lps)==0
            return
        end
        # Autocatalysis
        if c3_agent.state != 1
            r = rand(Uniform(0, 1))    
            if (r < model.k_cat)
                u = 2π * rand(model.rng) - π
                v = π * rand(model.rng)
                x = 0.05*cos(u) * sin(v)
                y = 0.05*sin(u) * sin(v)
                z = 0.05*cos(v)
                positions = c3_agent.pos 
                if (0.0 <= positions[1] <= 1.0) && 0.0 <= positions[2] <=1.0 && 0.0 <= positions[3] <=1.0
                    id = nextid(model)
                    add_agent_pos!( c3(id, positions, 0, 0, 0) , model)
                    return
                end
            end
        end
    
    end

    return
end

@everywhere function lp_step!(lp_agent, model)
    #Binding sites are static
    return
end

@everywhere function model_step!(model)
    #Recruite c3b:: tickover
    r = rand(Uniform(0, 1))
    if length(model.lps)!=0
        id = model.lps[rand(1:end)]
    else 
        return        
    end

    if (r < model.k_tick)
        x, y, z = model.agents[id].pos
        model.agents[id].engaged = 1
        deleteat!(model.lps, findall(x->x==id, model.lps))
        id = nextid(model)
        vel = sincos(2π * rand(model.rng)) .* model.surface_diffusion
        vel = (0,0,0)
        add_agent_pos!( c3(id, (x, y, z), 1, 0, 0) , model)
        return
    end

end

@everywhere nsteps = 1400 #Total steps
@everywhere runs = 500 #Number of simulations

@everywhere models = [initialize_model(seed = x) for x in rand(UInt16, runs)]
@everywhere unbound_sites(model) = length(model.lps)
@everywhere active_agents(model) = length([a.id for a in allagents(model) if a.state==0])
@everywhere mdata = [unbound_sites, active_agents]
@everywhere adf, mdf, mds = ensemblerun!(models, agent_step!, model_step!, nsteps; mdata, parallel = true)

using CSV
CSV.write(string(n_agents)*"agents_data.csv", mdf) 
