module GenAgentsMod
using LinearAlgebra
using Distances

function generate_agents_random(K::Int, r_comm::Float64, r_coll::Float64)
    agents = []
    new_agent = [10.0, 25.0];
    push!(agents, new_agent);
    new_agent = [40.0, 90.0];
    push!(agents, new_agent);
    new_agent = [75.0, 50.0];
    push!(agents, new_agent);
    new_agent = [5.0, 85.0];
    push!(agents, new_agent);
    new_agent = [60.0, 15.0];
    push!(agents, new_agent);
    new_agent = [20.0, 70.0];
    push!(agents, new_agent);
    new_agent = [35.0, 40.0];
    push!(agents, new_agent);
    new_agent = [90.0, 10.0];
    push!(agents, new_agent);
    return agents
end

function generate_goals_random(K::Int, r_comm::Float64, r_coll::Float64)
    agents = []
    new_agent = [55.0, 80.0];
    push!(agents, new_agent);
    new_agent = [20.0, 60.0];
    push!(agents, new_agent);
    new_agent = [85.0, 30.0];
    push!(agents, new_agent);
    new_agent = [10.0, 50.0];
    push!(agents, new_agent);
    new_agent = [45.0, 95.0];
    push!(agents, new_agent);
    new_agent = [30.0, 25.0];
    push!(agents, new_agent);
    new_agent = [65.0, 70.0];
    push!(agents, new_agent);
    new_agent = [5.0, 40.0];
    push!(agents, new_agent); 
    return agents
end

function generate_agents_head_on(K::Int,r_comm::Float64)
    agents = []
    new_agent = [0.0 50.0];
    push!(agents, new_agent)
    new_agent = [100.0,50.0];
    push!(agents, new_agent)  
    new_agent = [50.0,0.0];
    push!(agents, new_agent)  
    new_agent = [100.0,50.0];
    push!(agents, new_agent)
    new_agent = [14.65,14.65];
    push!(agents, new_agent) 
    new_agent = [14.65,85.35];
    push!(agents, new_agent) 
    new_agent = [85.35,85.35];
    push!(agents, new_agent) 
    new_agent = [85.35,14.65];
    push!(agents, new_agent) 
    return agents
end

function generate_agents_crossing(K::Int, r_comm::Float64)
    agents = []
    new_agent = [100.0, 50.0];
    push!(agents, new_agent)
    new_agent = [0.0, 50.0];
    push!(agents, new_agent)  
    new_agent = [50.0, 100.0];
    push!(agents, new_agent) 
    new_agent = [50.0, 0.0];
    push!(agents, new_agent) 
    new_agent = [14.65,14.65];
    push!(agents, new_agent) 
    new_agent = [14.65,85.35];
    push!(agents, new_agent) 
    new_agent = [85.35,85.35];
    push!(agents, new_agent) 
    new_agent = [85.35,14.65];
    push!(agents, new_agent)  
    return agents
end

function generate_goals_crossing(K::Int, r_comm::Float64)
    agents = []
    new_agent = [0.0, 50.0];
    push!(agents, new_agent)
    new_agent = [100.0, 50.0];
    push!(agents, new_agent)  
    new_agent = [50.0, 0.0];
    push!(agents, new_agent)  
    new_agent = [50.0, 100.0];
    push!(agents, new_agent)  
    new_agent = [85.35,85.35];
    push!(agents, new_agent) 
    new_agent = [85.35,14.65];
    push!(agents, new_agent) 
    new_agent = [14.65,14.65];
    push!(agents, new_agent) 
    new_agent = [14.65,85.35];
    push!(agents, new_agent)
    return agents
end
end