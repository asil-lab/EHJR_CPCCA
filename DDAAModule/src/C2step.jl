module C2stepMod
using LinearAlgebra
using JuMP, SCS, MosekTools #, GLPK, MosekTools

function localization(all_agents,edges,vars)
    Sigma1e = inv(vars.Sigma1*I(2));
    Sigma2e = inv(vars.Sigma2*I(2));
    K = vars.k_c+vars.k_nc;

    # model = Model(SCS.Optimizer)
    model = Model(Mosek.Optimizer)
    @variable(model, e[1:vars.dim, 1:length(edges)])
    @variable(model, p[1:vars.dim, 1:K]>=-1000)

    @expression(model, j1, sum(0.5*(p[:,a]-(all_agents[a].pos+all_agents[a].est_vel))'*Sigma1e*(p[:,a]-(all_agents[a].est_pos+all_agents[a].est_vel)) for a in eachindex(all_agents)) )
    @expression(model, j2, sum(0.5*(edges[b].w-e[:,b])'*Sigma2e*(edges[b].w-e[:,b]) for b in eachindex(edges)))
    @constraint(model, [b in eachindex(edges)], e[:,b]==p[:,edges[b].j]-p[:,edges[b].k])
    @objective(model, Min, j1+j2)
    set_silent(model)

    optimize!(model)

    p_prev = [all_agents[a].est_pos for a in eachindex(all_agents)]

    # # update our estimation
    pos = value.(p);
    vel = [pos[:,a] - p_prev[a] for a in eachindex(all_agents)];
    # error = 0;
    for a in eachindex(all_agents)
    #     error = error + norm(pos[:,a]-all_agents[a].pos)
        all_agents[a].est_vel = vec(vel[a])
        all_agents[a].est_pos = vec(pos[:,a])
        # all_agents[a].vel = vec(vel[a])
    end
    return all_agents
end

function control_strat_centralized(all_agents, vars)
    K = size(all_agents,1)
    model = Model(SCS.Optimizer)
    # model = Model(Mosek.Optimizer)

    @variable(model, u[1:vars.dim, 1:K, 1:vars.Thor])
    @variable(model, p[1:vars.dim, 1:K, 1:vars.Thor+1])
    # @expression(model, J4, )

    @objective(model, Min, sum(1/(norm(all_agents[a].est_pos-all_agents[a].goal))* (p[:,a,t]-all_agents[a].goal)'*I(vars.dim)*(p[:,a,t]-all_agents[a].goal) for a in eachindex(all_agents[1:vars.k_c]), t in 2:vars.Thor+1))
    # +sum((u[:,a,t]-u[:,a,t-1])'*I(vars.dim)*(u[:,a,t]-u[:,a,t-1]) for a in eachindex(all_agents[1:vars.k_c]), t in 2:vars.Thor)

    @constraint(model, [a in eachindex(all_agents[1:vars.k_c]), t in 1:vars.Thor], p[:,a,t+1] == p[:,a,t] + u[:,a,t])
    @constraint(model, [a in eachindex(all_agents[1:vars.k_c])], p[:,a,1] == all_agents[a].est_pos)
    @constraint(model, [a in eachindex(all_agents[1:vars.k_c]), t in 1:vars.Thor], u[:,a,t]'*u[:,a,t]<=vars.umax^2)
    
    # collision constraint
    @constraint(model, [k in eachindex(all_agents), j in eachindex(all_agents), t in 2:vars.Thor+1; k!=j],
    norm(all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t])+
    ((all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t])/norm(all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t]))'*
    ((p[:,k,t]-p[:,j,t])-(all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t]))-vars.r_coll>=0)

    # nc agents on nom traj
    @constraint(model, [b in eachindex(all_agents[vars.k_c+1:K]), t in 2:vars.Thor+1], p[:,b,t] == all_agents[b].nom_traj[:,t])
    set_silent(model)

    # print(model)
    optimize!(model)

    # Update velocity, nominal trajectory
    for a in eachindex(all_agents[1:vars.k_c])
        all_agents[a].vel = value.(u[:,a,1])
        all_agents[a].est_vel = value.(u[:,a,1])
        all_agents[a].nom_traj[:,1:end-1] = value.(p[:,a,2:end])
        all_agents[a].nom_traj[:,end] = value.(p[:,a,end])
    end

    return all_agents
end

function control_strat_relaxed(all_agents, vars)
    K = size(all_agents,1)
    model = Model(Mosek.Optimizer)

    @variable(model, u[1:vars.dim, 1:K, 1:vars.Thor])
    @variable(model, p[1:vars.dim, 1:K, 1:vars.Thor+1])
    @variable(model, D[1:vars.dim, 1:vars.dim,1:2, 1:vars.Thor])
    # @expression(model, J4, )

    @objective(model, Min, sum(1/(norm(all_agents[a].pos-all_agents[a].goal))* (p[:,a,t]-all_agents[a].goal)'*I(vars.dim)*(p[:,a,t]-all_agents[a].goal) for a in eachindex(all_agents[1:vars.k_c]), t in 2:vars.Thor+1))
    # +sum((u[:,a,t]-u[:,a,t-1])'*I(vars.dim)*(u[:,a,t]-u[:,a,t-1]) for a in eachindex(all_agents[1:vars.k_c]), t in 2:vars.Thor)

    @constraint(model, [a in eachindex(all_agents[1:vars.k_c]), t in 1:vars.Thor], p[:,a,t+1] == p[:,a,t] + u[:,a,t])
    @constraint(model, [a in eachindex(all_agents[1:vars.k_c])], p[:,a,1] == all_agents[a].pos)
    @constraint(model, [a in eachindex(all_agents[1:vars.k_c]), t in 1:vars.Thor], u[:,a,t]'*u[:,a,t]<=vars.umax^2)
    
    # collision constraint
    # @constraint(model, [k in eachindex(all_agents[1:vars.k_c]), j in eachindex(all_agents), t in 2:vars.Thor+1; k!=j],)
    a=1
    for k in eachindex(all_agents[1:vars.k_c])
        for j in eachindex(all_agents[1:vars.k_c])
            if k!=j
                for t = 1:vars.Thor
                    delta = p[:,k,t+1]-p[:,j,t+1]
                    @constraint(model, [D[:,:,a,t] delta; delta' 1]>=0, PSDCone())

                    @constraint(model, tr(D[:,:,a,t])>=vars.r_coll)
                end
            end
        end
        a=a+1;
    end
    # nc agents on nom traj
    set_silent(model)

    # print(model)
    optimize!(model)

    # Update velocity, nominal trajectory
    for a in eachindex(all_agents[1:vars.k_c])
        all_agents[a].vel = value.(u[:,a,1])
        all_agents[a].est_vel = value.(u[:,a,1])
        all_agents[a].nom_traj[:,1:end-1] = value.(p[:,a,2:end])
        all_agents[a].nom_traj[:,end] = value.(p[:,a,end])
    end

    return all_agents
end

function control_strat_ref_traj(all_agents, vars, T) 
    K = size(all_agents,1)
    # model = Model(SCS.Optimizer)
    model = Model(Mosek.Optimizer)

    @variable(model, u[1:vars.dim, 1:K, 1:vars.Thor])
    @variable(model, p[1:vars.dim, 1:K, 1:vars.Thor+1])
    @expression(model, J4, sum((p[:,a,t]-all_agents[a].ref_traj[:,T+t-1])'*I(vars.dim)*(p[:,a,t]-all_agents[a].ref_traj[:,T+t-1]) for a in eachindex(all_agents[1:vars.k_c]), t in 2:vars.Thor+1))
    @objective(model, Min, J4)

    @constraint(model, [a in eachindex(all_agents[1:vars.k_c]), t in 1:vars.Thor], p[:,a,t+1] == p[:,a,t] + u[:,a,t])
    @constraint(model, [a in eachindex(all_agents[1:vars.k_c])], p[:,a,1] == all_agents[a].pos)
    @constraint(model, [a in eachindex(all_agents[1:vars.k_c]), t in 1:vars.Thor], u[:,a,t]'*u[:,a,t]<=vars.umax^2)
    
    # collision constraint
    @constraint(model, [k in eachindex(all_agents[1:vars.k_c]), j in eachindex(all_agents), t in 2:vars.Thor+1; k!=j],
    norm(all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t])+
    ((all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t])/norm(all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t]))'*
    ((p[:,k,t]-p[:,j,t])-(all_agents[k].nom_traj[:,t]-all_agents[j].nom_traj[:,t]))-vars.r_coll>=0)

    # nc agents on nom traj
    set_silent(model)

    # print(model)
    optimize!(model)

    # Update velocity, nominal trajectory
    for a in eachindex(all_agents[1:vars.k_c])
        all_agents[a].vel = value.(u[:,a,1])
        all_agents[a].est_vel = value.(u[:,a,1])
        all_agents[a].nom_traj[:,1:end-1] = value.(p[:,a,2:end])
        all_agents[a].nom_traj[:,end] = value.(p[:,a,end])
    end

    return all_agents
end

function centralized_2_step(all_agents,edges,var_set,n)
    all_agents_update = localization(all_agents, edges,var_set)
    all_agents_new = control_strat_centralized(all_agents_update, var_set)
    return all_agents_new
end

end