module DJointMod
using LinearAlgebra
using JuMP,SCS, Ipopt #, GLPK, MosekTools

include("neighborhoodmod.jl")
include("logtofile.jl")

using .logtofile
using .neighborhoodmod

mutable struct Results
    positions::Vector{Float64}
    controls::Vector{Float64}
    function Results(positions::Vector{Float64}, controls::Vector{Float64})
        return new(positions, controls)
    end
end

function index_z(k,j,t,Thor)
    if k < j 
        idxk1 = t
        idxk2 = Thor+1+t
        idxj1 = 2*(Thor+1)+t
        idxj2 = 3*(Thor+1)+t
    else
        idxj1 = t
        idxj2 = Thor+1+t
        idxk1 = 2*(Thor+1)+t
        idxk2 = 3*(Thor+1)+t
    end
    return idxk1, idxk2, idxj1, idxj2
end





function solve_for_z(all_agents,vars)
    res = 0
    conv = false
    for agent in all_agents
        for edge in agent.neighborhood
            edge.z = 0.5*(edge.yk+edge.yj);
            # compute lambda
            edge.lambdak = edge.lambdak +vars.rho * (edge.yk - edge.z)
            edge.lambdaj = edge.lambdaj +vars.rho * (edge.yj - edge.z)
            res = res + norm(edge.yk-edge.z)
        end
    end
    if res < 10^-6
        conv = true
    end
    return all_agents, conv
end

function z_collision(all_agents, edge, vars)
    neighbor = all_agents[edge.j]
    agent = all_agents[edge.k]

    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, z[1:size(edge.yk,1)]>=-1000)
    @objective(model, Min, edge.lambdak'*(edge.yk-z)+edge.lambdaj'*(edge.yj-z) +vars.rho/2*((edge.yk - z)'*(edge.yk - z)+(edge.yj - z)'*(edge.yj - z)))

    deadlock = false
    #detect head on collision
    for t in 2:vars.Thor+1
        idxk1, idxk2, idxj1, idxj2 = index_z(edge.k, edge.j, t,vars.Thor)
        deltap = [edge.yk[idxk1], edge.yk[idxk2],0]-[edge.yk[idxk1-1], edge.yk[idxk2-1],0]
        eta = (all_agents[edge.k].nom_traj[:,t]-all_agents[edge.j].nom_traj[:,t])/norm(all_agents[edge.k].nom_traj[:,t]-all_agents[edge.j].nom_traj[:,t])
        test = cross(vcat(eta,0), deltap)
        if sum(abs.(test))<10^-6
            deadlock = true
        end
    end

        
    for t in 2:vars.Thor+1
        idxk1, idxk2, idxj1, idxj2 = index_z(edge.k, edge.j, t,vars.Thor)
        deltabar = agent.nom_traj[:,t]-neighbor.nom_traj[:,t]
        if neighbor.cooperative == 1
            delta = [z[idxk1]-z[idxj1], z[idxk2]-z[idxj2]]
        else 
            delta = [z[idxk1]-neighbor.nom_traj[1,t], z[idxk2]-neighbor.nom_traj[2,t]]
        end
        if deadlock
            Rmat = [cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]
            eta = Rmat*(deltabar/norm(deltabar))
        else
            eta = (deltabar/norm(deltabar))
        end
        @constraint(model, norm(deltabar)+eta'*(delta-deltabar)>=vars.r_coll)
    end


    # # solve for z
    optimize!(model)
    return vec(value.(z))
end


function solve_for_z_collision(all_agents, vars)
    res = 0
    conv = false
    for agent in all_agents
        for edge in agent.neighborhood
            temp = edge.z
            coll = 0
            for t in 1:vars.Thor
                idxk1, idxk2, idxj1, idxj2 = index_z(edge.k,edge.j,t,vars.Thor)
                if any(norm([edge.yk[idxk1] edge.yk[idxk2]]-[edge.yk[idxj1] edge.yk[idxj2]]) <= vars.r_coll)
                    coll = coll+1;
                end
            end
            if coll != 0 # collision in traj
                edge.z = z_collision(all_agents, edge, vars)
            else
                edge.z = 0.5*(edge.yk+edge.yj);
            end
            # compute lambda
            edge.lambdak = edge.lambdak +vars.rho * (edge.yk - edge.z)
            edge.lambdaj = edge.lambdaj +vars.rho * (edge.yj - edge.z)
            res = res + norm(temp-edge.z)
        end
    end
    if res < 10^-4
        conv = true
    end
    return all_agents, conv
end

function set_variables_y(agent, model, vars)
    n_neighbors = length(agent.neighborhood)
    @variable(model, p[1:vars.dim, 1:n_neighbors+1, 1:vars.Thor+1])
    @variable(model, u[1:vars.dim, 1:vars.Thor])
    return model, p, u
end

function localization_cost(agent, vars, p, model)
    Sigma1e = inv(vars.Sigma1*I(2));
    Sigma2e = inv(vars.Sigma2*I(2));
    Sigma2 = Sigma2e[1];
    @expression(model, j1, (p[:,1,1]-(agent.est_pos+agent.est_vel))'*Sigma1e*(p[:,1,1]-(agent.est_pos+agent.est_vel)))
   
    @NLexpression(model, j2, 
        Sigma2*sum(
        (agent.neighborhood[b].e^2- 
        ((p[1,b+1,1]-p[1,1,1])^2 + (p[2,b+1,1]-p[2,1,1])^2))^2
        for b in eachindex(agent.neighborhood)
        ));

    return model, j1, j2
end

function MPC_cost(agent, vars, p, u, model)
    weight = 1; #norm(agent.est_pos-agent.goal)
    if weight <= 10^-5
        weight = 10^-5
    end
    @expression(model, j3, sum(1/(vars.Thor*weight)* (p[:,1,t]-agent.goal)'*I(vars.dim)*(p[:,1,t]-agent.goal) for t in 2:vars.Thor+1))
    return model, j3
end

function lagrangian_cost(agent,vars, p, model)
    @expression(model, jf1, sum(agent.neighborhood[b].lambdak'*[vec(p[:,1,:]'); vec(p[:,b+1,:]')]+vars.rho/2*([vec(p[:,1,:]'); vec(p[:,b+1,:]')]-agent.neighborhood[b].z)'*([vec(p[:,1,:]'); vec(p[:,b+1,:]')]-agent.neighborhood[b].z) for b in eachindex(agent.neighborhood) if agent.neighborhood[b].k < agent.neighborhood[b].j))
    @expression(model, jf2, sum(agent.neighborhood[b].lambdak'*[vec(p[:,b+1,:]'); vec(p[:,1,:]')]+vars.rho/2*([vec(p[:,b+1,:]'); vec(p[:,1,:]')]-agent.neighborhood[b].z)'*([vec(p[:,b+1,:]'); vec(p[:,1,:]')]-agent.neighborhood[b].z) for b in eachindex(agent.neighborhood) if agent.neighborhood[b].k > agent.neighborhood[b].j))
    return model, jf1, jf2
end


function dynamic_constraint(vars, model, p, u)
    @constraint(model, [t in 1:vars.Thor], p[:,1,t+1] == p[:,1,t] + u[:,t])
    @constraint(model, [t in 1:vars.Thor], u[:,t]'*u[:,t]<=vars.umax^2)
    return model
end

function collision_constraint_y(all_agents, agent,vars, p, model)
    for (j,edge) in enumerate(agent.neighborhood)
        neighbor = all_agents[edge.j]
        for t in 2:vars.Thor+1
            deltabar = agent.nom_traj[:,t]-neighbor.nom_traj[:,t]
            delta = p[:,1,t]-p[:,j+1,t]
            @constraint(model, norm(deltabar)+(deltabar/norm(deltabar))'*(delta-deltabar)>=vars.r_coll)
        end
    end
    return model
end

function solve_for_y(all_agents, vars)
    for agent in all_agents
        model = Model(Ipopt.Optimizer)
        set_silent(model)
        
        model, p, u = set_variables_y(agent, model, vars)
        model, j1, j2 = localization_cost(agent, vars, p, model)
        model, j3 = MPC_cost(agent, vars, p, u, model)
        model, jf1, jf2 = lagrangian_cost(agent, vars, p, model)

        # model = localization_constraint(agent, model, p, e)
        model = dynamic_constraint(vars, model, p, u)
        model = collision_constraint_y(all_agents, agent, vars, p, model)
        
        @NLobjective(model, Min, vars.alpha*(j1+j2)+(1-vars.alpha)*j3+jf1+jf2)
        optimize!(model)

        for (b,edge) in pairs(IndexLinear(),agent.neighborhood)
            if edge.k < edge.j
                edge.yk = [vec(value.(p[:,1,:])'); vec(value.(p[:,b+1,:])')]
            else
                edge.yk = [vec(value.(p[:,b+1,:])'); vec(value.(p[:,1,:])')]
            end
        end
        agent.q = value.(p[:,1,:]);
        agent.vel = vec(value.(u[:,1]));
    end
    return all_agents
end

function save_values(all_agents, i, results_coll)
    for (k,agent) in enumerate(all_agents)
        results_coll[k,i] = Results(agent.q[:,1], agent.vel)
    end
    return results_coll
end

function joint_algorithm(all_agents, edges, vars, results_coll)
    dimz = vars.dim*(vars.Thor+1)*2
    edges = neighborhoodmod.initialize_consensus_variables(edges, dimz)
    all_agents = neighborhoodmod.set_agent_neighborhood(all_agents, edges)
    final_iter = vars.iter;
    for i in 1:vars.iter
        agents_y = solve_for_y(all_agents, vars)
        agents_y_u = neighborhoodmod.communicate_y_neighborhood(agents_y)
        all_agents, converged = solve_for_z(agents_y_u, vars)
        results_coll = save_values(all_agents, i, results_coll)
        for agent in all_agents
            agent.nom_traj = agent.q
        end
        if converged
            final_iter = i
            break; 
        end   
    end

    for agent in all_agents
        agent.est_pos = agent.q[:,1]
        agent.est_vel = agent.vel
    end
    return all_agents, results_coll, final_iter
end

function calculate_error(results_coll, final_iter, all_agents)
    errors_pos = fill(NaN, size(results_coll,2));
    errors_ctrl = fill(NaN, size(results_coll,2));
    for i in 1:final_iter
        errors_pos[i] = 0.00;
        errors_ctrl[i] = 0.00;
        for k in eachindex(results_coll[:,i])
            errors_pos[i] += norm(results_coll[k,i].positions-all_agents[k].pos)
            errors_ctrl[i] += norm(results_coll[k,i].controls-results_coll[k,final_iter].controls)
        end
    end
    return errors_pos, errors_ctrl
end

function distributed_joint(all_agents, edges, vars)
    results_coll = Array{Results}(undef, size(all_agents,1),vars.iter)
    all_agents_update, results_coll, final_iter = joint_algorithm(all_agents, edges, vars,results_coll)
    results_pos, results_ctrl = calculate_error(results_coll, final_iter, all_agents_update)
    return all_agents_update,results_pos, results_ctrl
end
end

