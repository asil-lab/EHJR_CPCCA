module CJointMod
using LinearAlgebra
using JuMP, SCS
using MosekTools

function joint_algorithm(all_agents, edges, vars)
    Sigma1e = inv(vars.Sigma1*I(2));
    Sigma2e = inv(vars.Sigma2*I(2));
    K = vars.k_c+vars.k_nc;

    model = Model(Mosek.Optimizer)
    @variable(model, e[1:vars.dim,1:length(edges)])
    @variable(model, p[1:vars.dim,1:K,1:vars.Thor+1])
    @variable(model, u[1:vars.dim, 1:K, 1:vars.Thor])

    @expression(model, j1, sum(0.5*(p[:,a,1]-(all_agents[a].est_pos+all_agents[a].est_vel))'*Sigma1e*(p[:,a,1]-(all_agents[a].est_pos+all_agents[a].est_vel)) for a in eachindex(all_agents)) )
    @expression(model, j2, sum(0.5*(edges[b].w-e[:,b])'*Sigma2e*(edges[b].w-e[:,b]) for b in eachindex(edges)))
    @constraint(model, [b in eachindex(edges)], e[:,b]==p[:,edges[b].j,1]-p[:,edges[b].k,1])
    
    weight = [norm(agent.est_pos-agent.goal) for agent in all_agents]
    for w in weight
        if w <= 10^-10
            w = 10^-10
        end
    end

    @expression(model, j3, sum(0.5*1/(vars.Thor*weight[a])* (p[:,a,t]-all_agents[a].goal)'*I(vars.dim)*(p[:,a,t]-all_agents[a].goal) for a in eachindex(all_agents[1:vars.k_c]), t in 2:vars.Thor+1))
    
    @constraint(model, [a in eachindex(all_agents[1:vars.k_c]), t in 1:vars.Thor], p[:,a,t+1] == p[:,a,t] + u[:,a,t])
    @constraint(model, [a in eachindex(all_agents[1:vars.k_c]), t in 1:vars.Thor], u[:,a,t]'*u[:,a,t]<=vars.umax^2)
    
    # collision constraint
    @constraint(model, [k in eachindex(all_agents[1:vars.k_c]), j in eachindex(all_agents), t in 2:vars.Thor+1; k!=j],
    norm(all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t])+
    ((all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t])/norm(all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t]))'*
    ((p[:,k,t]-p[:,j,t])-(all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t]))-vars.r_coll>=0)

    # nc agents on nom traj
    @constraint(model, [a in eachindex(all_agents[vars.k_c+1:end]), t in 1:vars.Thor+1], p[:,a,t] == all_agents[a].nom_traj[:,t])

    @objective(model, Min, vars.alpha*(j1+j2)+(1-vars.alpha)*j3)
    set_silent(model)

    optimize!(model)

    # update our estimation
    
    p_prev = [all_agents[a].pos for a in eachindex(all_agents)]
 
    pos = value.(p);
    vel = [pos[:,a,1] - p_prev[a] for a in eachindex(all_agents)];
    # error = 0;
    for a in eachindex(all_agents)
    #     error = error + norm(pos[:,a]-all_agents[a].pos)
        all_agents[a].est_vel = vec(vel[a])
        all_agents[a].est_pos = vec(pos[:,a,1])
    end
    for a in eachindex(all_agents[1:vars.k_c])
        all_agents[a].vel = value.(u[:,a,1])
        all_agents[a].est_vel = value.(u[:,a,1])
        all_agents[a].nom_traj[:,1:end-1] = value.(p[:,a,2:end])
        all_agents[a].nom_traj[:,end] = value.(p[:,a,end])
    end

    return all_agents
end

function centralized_joint(all_agents, edges, vars)
    all_agents_update = joint_algorithm(all_agents, edges, vars)
    return all_agents_update
end

end