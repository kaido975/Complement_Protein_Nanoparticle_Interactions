using Agents, Random, Agents.Pathfinding
using Distributions
using InteractiveDynamics
using GLMakie
GLMakie.activate!()
using LinearAlgebra
import Plots

include("angle_compute.jl")
steps = 1400

@agent complement_agent ContinuousAgent{3} begin
    type::Symbol
    engaged::Int
    state::Int
    age::Float32
end
const v0 = (0.0, 0.0, 0.0)
c3(id, pos, engaged, state, age) = complement_agent(id, pos, v0, :c3, engaged, state, age)
lp(id, pos, engaged, state, age) = complement_agent(id, pos, v0, :lp, engaged, state, age)

r = 0.5

x = collect(range(0, 1.0, length = 80))
y = collect(range(0, 1.0, length = 80))
z = collect(range(0, 1.0, length = 80))
z = reshape(z, (1, 1, size(z, 1)))
y = reshape(y, (1, size(y, 1)))

cx = r
cy = r
cz = r


# The two lines below could be merged, but I stored the mask
# for code clarity.
ϵ = 0.02
mask = r^2 + ϵ/12  .<= (x .- cx).^2 .+ (y .- cy).^2 .+ (z .- cz).^2 .<= r^2 + ϵ
findall(mask.==1)
# mask = BitArray(trues(10,10,10))
n = 50
u = range(-π, π; length = n)
v = range(0, π; length = n)
x = r*cos.(u) * sin.(v)'
y = r*sin.(u) * sin.(v)'
z = r*ones(n) * cos.(v)'

#Distribution of agents on the grid
n_agents = 200

i = range(0, n_agents-1, length = n_agents) .+ 0.5
ϕ = acos.(1 .- 2*i/n_agents)
ω = (1 + 5^0.5)/2   
θ = 2*pi * i / ω

x_lp, y_lp, z_lp = r*cos.(θ) .* sin.(ϕ) .+ cx, r*sin.(θ) .* sin.(ϕ) .+ cy, r*cos.(ϕ) .+ cz

function initialize_model(;
    seed = 42,  ## seed for random number generator,
    walk = mask,
    speed = 0.35,
    dt = 0.1, 
    k_tick = 4e-3, 
    k_cat = 6e-2,
    decay_time = 5.0
)

    lps = Int.(collect(range(1, n_agents, length = n_agents)))

    rng = MersenneTwister(seed)
    space = ContinuousSpace((1.0, 1.0, 1.0); periodic = false)

    properties = (
        walkable_path = AStar(space; walkmap = mask),
        speed = speed,
        dt = dt,
        k_tick  =k_tick,
        k_cat = k_cat, 
        decay_time = decay_time,
        lps = lps,
        counter = [0]
    )

    model = ABM(complement_agent, space; rng, properties)
    for i in 1:n_agents
        id = nextid(model)
        agent = lp(id, (x_lp[i], y_lp[i], z_lp[i]), 0, 1, 0)
        add_agent_pos!(agent, model)  
    end

    return model
end


function agent_step!(agent::complement_agent, model)
    if agent.type == :c3
        c3_step!(agent, model)
    else
        lp_step!(agent, model)
    end
end

function c3_step!(c3_agent, model)

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
        move_along_route!(c3_agent, model, model.walkable_path, model.speed, model.dt)        
        return  
    else 
        if length(model.lps)==0
            return
        end
        if c3_agent.state != 1
            r_num = rand(Uniform(0, 1))    
            if (r_num < model.k_cat)
                positions = c3_agent.pos 
                id = nextid(model)
                add_agent_pos!( c3(id, positions, 0, 0, 0) , model)
                return
            end
        end
    
    end

    return
end

function lp_step!(lp_agent, model)
    return
end

function model_step!(model)

    #Recruite c3:: tickover
    model.counter[1] += 1

    r_num = rand(Uniform(0, 1))
    if length(model.lps)!=0
        id = model.lps[rand(1:end)]
    else 
        return        
    end

    if (r_num < model.k_tick)
        x, y, z = model.agents[id].pos
        model.agents[id].engaged = 1
        deleteat!(model.lps, findall(x->x==id, model.lps))
        id = nextid(model)
        vel = sincos(2π * rand(model.rng)) .* model.speed
        vel = (0,0,0)
        add_agent_pos!( c3(id, (x, y, z), 1, 0, 0) , model)
        return
    end

end

#Scheduler to find active agents
cus_scheduler(model) = (a.id for a in allagents(model) if a.state==0)

function static_preplot!(ax, model)
    GLMakie.surface!(ax,
        x .+ r,
        y .+ r,
        z .+ r;
        colorrange = (0, 0),
        colormap = :vermeer,

        ambient = Vec3f(0., 0, 0),
        diffuse = Vec3f(0., 0., 0),
        specular = Vec3f(0, 0, 0),
        shininess = 0
    )
end

model = initialize_model()
groupcolor(a) = a.type == :c3 ? a.state==1 ? :black : a.engaged ==1 ? :red : :blue : :lightgreen 
groupsize(a) =  a.type == :c3 ? 0.01 : 0.012
position_agent(a) = a.type == :c3 ? a.engaged ==1 ?  (a.pos .- (r,r,r)).*0.02 : 0.0 : 0.0

f = abmvideo(
    "quasi_2d.mp4",
    model, agent_step!, model_step!;
    figure = (resolution = (1600, 1600),),
    frames = steps,
    framerate = 20,
    ac = groupcolor,
    as = groupsize,
    offset = position_agent,
    static_preplot!,
    title = "C3b spreading"
)

