__precompile__()
module DDAAModule

export run_simulation #, initialize_agents, Edge, Agent # fill in exported names here

using CSV, DataFrames
using LinearAlgebra
using Tables
using Random, Distributions  # hide+
using Distances
using Plots
using MosekTools
using JLD2
# using Plots


include("C2step.jl")
include("D2step.jl")
include("CJointMod.jl")
include("DJointMod.jl")
include("GenAgentsMod.jl")
include("neighborhoodmod.jl")
include("logtofile.jl")
using .GenAgentsMod, .logtofile, .neighborhoodmod
using .C2stepMod, .D2stepMod, .CJointMod, .DJointMod 


# filename = "results_joint_alpha1.jld2"
# noiselevels = [Sigma1, Sigma2];

mutable struct Edge
    k:: Int
    j:: Int
    r:: Float64
    w:: Vector{Float64}
    e:: Float64
    z:: Vector{Float64}
    lambdak::Vector{Float64}
    lambdaj::Vector{Float64}
    yk:: Vector{Float64}
    yj:: Vector{Float64}
    function Edge(k::Int, j::Int, r::Float64, w::Vector{Float64}, e::Float64)
        return new(k,j,r,w, e, Vector{Float64}(undef,0),Vector{Float64}(undef,0),Vector{Float64}(undef,0), Vector{Float64}(undef,0),Vector{Float64}(undef,0))
    end
end

mutable struct Agent
    id::Int
    pos::Vector{Float64}
    vel::Vector{Float64}
    goal::Vector{Float64}
    p0::Vector{Float64}
    cooperative::Int
    est_pos::Vector{Float64}
    est_vel::Vector{Float64}
    neighborhood:: Vector{Edge}
    q::Matrix{Float64}
    nom_traj::Matrix{Float64}
    ref_traj::Matrix{Float64}
    function Agent(id::Int, pos::Vector{Float64},vel::Vector{Float64},goal::Vector{Float64},p0::Vector{Float64},cooperative::Int,traj::Matrix{Float64},ref_traj::Matrix{Float64})
        return new(id,pos,vel,goal,p0,cooperative,pos,vel,Vector{Edge}(undef, 0),traj,traj,ref_traj)
    end
end

function initialize_agents(var_set,n_cooperative = 4, n_noncooperative = 2, nsteps = 100)
    a = Vector{Agent}(undef, n_cooperative+n_noncooperative)
    # positions = GenAgentsMod.generate_agents_random(n_cooperative+n_noncooperative, Float64(var_set.r), Float64(var_set.r))
    # goals = GenAgentsMod.generate_goals_random(n_cooperative+n_noncooperative, Float64(var_set.r), Float64(var_set.r))
    positions = GenAgentsMod.generate_agents_crossing(n_cooperative+n_noncooperative, Float64(var_set.r))
    goals = GenAgentsMod.generate_goals_crossing(n_cooperative+n_noncooperative, Float64(var_set.r))
    if n_noncooperative != 0
        pos_nc = reduce(vcat,positions[n_cooperative+1:end])
        vel_nc = zeros(var_set.dim*n_noncooperative,)
    end
    b=1
    for i in n_cooperative+1:+n_cooperative+n_noncooperative
        id = i;
        pos = positions[i] 
        goal = goals[i] 
        vel = (goal-pos)/norm(goal-pos)*var_set.umax;
        vel_nc[b:b+1] = vel;
        nom_traj = zeros(Float64,(2,var_set.Thor+1))
        for t = 0:var_set.Thor
            nom_traj[:,t+1] = pos+vel*t
        end
        b=b+2
        a[i] = Agent(id, pos,vel, goal,pos,2,nom_traj,nom_traj)
    end
    
    for i in 1:n_cooperative
        id = i;
        pos = positions[i]  
        goal = goals[i]; 
        vel = [0.0, 0.0];
        nom_traj = zeros(Float64,(2,var_set.Thor+1))
        for t = 0:var_set.Thor
            nom_traj[:,t+1] = pos+vel*t
        end
        ref_traj = zeros(Float64, (2,nsteps+var_set.Thor))
        ref_traj[:,1] = pos
        for t = 1:nsteps-1+var_set.Thor
            if norm(ref_traj[:,t]-goal) < var_set.umax
                ref_traj[:,t+1] = goal
            else
                ref_traj[:,t+1] = ref_traj[:,t]+(goal-ref_traj[:,t])/norm(goal-ref_traj[:,t])*var_set.umax
            end
        end
        a[i] = Agent(id, pos,vel, goal,pos,1,nom_traj, ref_traj)
    end
    return a
end

function move_agents(all_agents,Distr1)
    for agent in all_agents
        setfield!(agent, :pos, agent.pos+agent.vel+rand(Distr1,(2,)));
        if agent.cooperative == 2
            setfield!(agent, :nom_traj, agent.nom_traj.+agent.vel);
        end
    end

    return all_agents
end

function generate_measurements(r,all_agents, Distr2)
    cur_pos = reduce(vcat,transpose.((p->p.pos).(all_agents)))
    pairs = pairwise(Euclidean(), cur_pos, dims=1)
    rd = deepcopy(pairs)
    w = zeros(Float64, (length(all_agents),length(all_agents),2))
    e = zeros(Float64, (length(all_agents),length(all_agents)))
    concat = Vector{Edge}(undef,0);
    for a in eachindex(all_agents)
        for b in a+1:eachindex(all_agents)[end]
            radius = pairs[a,b]
            if radius <= r
                # edge exists
                if (all_agents[a].cooperative == 1) && (all_agents[b].cooperative == 1)
                    w[a,b,:] = cur_pos[b,:] - cur_pos[a,:] + rand(Distr2,(2,1));
                    w[b,a,:] = cur_pos[a,:] - cur_pos[b,:] + rand(Distr2,(2,1));
                    e[a,b] = norm(cur_pos[b,:] - cur_pos[a,:])+rand(Distr2);
                    e[b,a] = norm(cur_pos[a,:] - cur_pos[b,:])+rand(Distr2);
                    concat = vcat(concat,Edge(a, b, rd[a,b], w[a,b,:],e[a,b]));
                    concat = vcat(concat,Edge(b, a, rd[b,a], w[b,a,:],e[b,a]));
                elseif (all_agents[a].cooperative == 1) && (all_agents[b].cooperative == 2)
                    w[a,b,:] = cur_pos[b,:]-cur_pos[a,:] + rand(Distr2,(2,1))
                    e[a,b] = norm(cur_pos[b,:] - cur_pos[a,:])+rand(Distr2);
                    concat = vcat(concat,Edge(a, b, rd[a,b], w[a,b,:],e[a,b]));
                elseif (all_agents[a].cooperative == 2) && (all_agents[b].cooperative == 1)
                    w[b,a,:] = cur_pos[a,:]-cur_pos[b,:] + rand(Distr2,(2,1))  
                    e[b,a] = norm(cur_pos[a,:] - cur_pos[b,:])+rand(Distr2);
                    concat = vcat(concat,Edge(b, a, rd[b,a], w[b,a,:],e[b,a]));
                end
            end
        end
    end
    return concat
end

function set_updated_state(all_agents, results)
    for agent in all_agents
        agent.est_pos = agent.pos;
        agent.est_vel = agent.vel;
    end
    return all_agents
end

function run_simulation(filename, r,nsteps, var_set)
    # var_set = constant_variables(r, r_coll, umax, Thor, dim, Sigma1, Distr1, Sigma2, Distr2,k_c,k_nc,rho,iter,alpha);
    all_agents_coll = Matrix{Agent}(undef, var_set.k_c+var_set.k_nc,nsteps)
    all_agents = initialize_agents(var_set,var_set.k_c,var_set.k_nc,nsteps);
    results_pos = zeros(Float64, (nsteps, var_set.iter))
    results_ctrl = zeros(Float64, (nsteps, var_set.iter))
    
    for n in 1:nsteps
        println("step $n")
        edges = generate_measurements(r, all_agents,var_set.Distr2); # returns all existing edges and meas [k, j, dkj, e1kj, e2kj]    
        # Agents_cent = C2stepMod.centralized_2_step(all_agents, edges,var_set,n)
        # Agents_dist = D2stepMod.distributed_2_step(all_agents, edges, var_set)
        # Agents_c_joint = CJointMod.centralized_joint(all_agents, edges, var_set)
        Agents_d_joint, results_pos_n, results_ctrl_n = DJointMod.distributed_joint(all_agents, edges, var_set)
        all_agents = copy(Agents_d_joint)
        results_pos[n,:] = results_pos_n;
        results_ctrl[n,:] = results_ctrl_n;
        all_agents_coll[:,n] = deepcopy(all_agents)
        move_agents(all_agents,var_set.Distr1);  # moves all_agents
    end
    @save(filename, var_set, all_agents_coll, results_pos, results_ctrl)
    return all_agents_coll, results_pos, results_ctrl
end


end # module
