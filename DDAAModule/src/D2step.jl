module D2stepMod
using LinearAlgebra
using JuMP, SCS
using MosekTools

include("neighborhoodmod.jl")
using .neighborhoodmod
include("logtofile.jl")
using .logtofile

function estimate_state_distributed(all_agents, edges, vars)
    Sigma1e = inv(vars.Sigma1*I(2));
    Sigma2e = inv(vars.Sigma2*I(2));
    dimz = vars.dim*3

    # error = zeros(1,iter);
    edges1 = copy(neighborhoodmod.initialize_consensus_variables(edges, dimz))
    all_agents = copy(neighborhoodmod.set_agent_neighborhood(all_agents, edges1))

    # Iterations
    for i = 1:vars.iter
        # Minimizing for q and y
        res = 0;
        for agent in all_agents
            n_neighbors = size(agent.neighborhood,1)

            # model = Model(SCS.Optimizer)
            model = Model(Mosek.Optimizer)
            @variable(model, p[1:vars.dim, 1:n_neighbors+1])
            @variable(model, e[1:vars.dim, 1:n_neighbors])

            @expression(model, j1, 0.5*(p[:,1]-(agent.est_pos+agent.est_vel))'*Sigma1e*(p[:,1]-(agent.est_pos+agent.est_vel)))
            @expression(model, j2, sum(0.5*(agent.neighborhood[b].w-e[:,b])'*Sigma2e*(agent.neighborhood[b].w-e[:,b]) for b in eachindex(agent.neighborhood)))
            @constraint(model, [b in eachindex(agent.neighborhood)], e[:,b]==p[:,b+1]-p[:,1])
            
            @expression(model, jf1, sum(agent.neighborhood[b].lambdak'*[e[:,b]' p[:,1]' p[:,b+1]']'+vars.rho/2*([e[:,b]' p[:,1]' p[:,b+1]']'-agent.neighborhood[b].z)'*([e[:,b]' p[:,1]' p[:,b+1]']'-agent.neighborhood[b].z) for b in eachindex(agent.neighborhood) if agent.neighborhood[b].k < agent.neighborhood[b].j))
            @expression(model, jf2, sum(agent.neighborhood[b].lambdak'*[-e[:,b]' p[:,b+1]' p[:,1]']'+vars.rho/2*([-e[:,b]' p[:,b+1]' p[:,1]']'-agent.neighborhood[b].z)'*([-e[:,b]' p[:,b+1]' p[:,1]']'-agent.neighborhood[b].z) for b in eachindex(agent.neighborhood) if agent.neighborhood[b].k > agent.neighborhood[b].j))

            @objective(model, Min, j1+j2+jf1+jf2)
            set_silent(model)

            optimize!(model)

            # Save for this agent
            for (b,edge) in pairs(IndexLinear(),agent.neighborhood)
                if edge.k < edge.j
                    temp = copy(edge.yk)
                    edge.yk = vec([value.(e[:,b])' value.(p[:,1])' value.(p[:,b+1])'])
                    res = res + norm(edge.yk-temp);
                else
                    temp = copy(edge.yk)
                    edge.yk = vec([-value.(e[:,b])' value.(p[:,b+1])' value.(p[:,1])'])
                    res = res + norm(edge.yk-temp);
                end
            end
            agent.q = value.(p[:,1])';
        end
        # Communicate values y to neighborhood
        all_agents = neighborhoodmod.communicate_y_neighborhood(all_agents)
  
        # logtofile.log_est_pos(all_agents, i)

        # Compute z on the edge (each agent)
        for agent in all_agents
            for edge in agent.neighborhood
                z = 0.5*(edge.yk+edge.yj);
                edge.z = z;
                # compute lambda
                edge.lambdak = edge.lambdak +vars.rho * (edge.yk - edge.z)
            end
        end
        if res < 10^-4
            print("location converged at iteration: ", i,'\n')
            break
        end
    end

    # Set the new positions to converged values
    for agent in all_agents
        agent.est_pos = vec(agent.q);
    end
    return all_agents
end

function detect_head_on(all_agents, edge, T)
    deadlock = false;
    for t in 2:T+1
        idxk1, idxk2, idxj1, idxj2 = index_z_mpc(edge.k, edge.j, t, T)
        
        eta = (all_agents[edge.k].nom_traj[:,t]-all_agents[edge.j].nom_traj[:,t])/norm(all_agents[edge.k].nom_traj[:,t]-all_agents[edge.j].nom_traj[:,t])

        deltap = [edge.yk[idxk1], edge.yk[idxk2],0]-[edge.yk[idxk1-1], edge.yk[idxk2-1],0]
        test = cross(vcat(eta,0), deltap)
        if sum(abs.(test))<0.00001
            deadlock = true
        end
    end
    if deadlock == true
        print("deadlock detected\n")
    end
    return deadlock
end

function index_z_mpc(k,j,t,Thor)
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

function check_collision(edge, vars)
    collision = 0
    for t in 1:vars.Thor
        idxk1, idxk2, idxj1, idxj2 = index_z_mpc(edge.k,edge.j,t,vars.Thor)
        if any(norm([edge.yk[idxk1] edge.yk[idxk2]]-[edge.yk[idxj1] edge.yk[idxj2]]) <= vars.r_coll)
            collision = collision+1;
        end
    end
    return collision
end

function collision_constraint(edge, agent,neighbor, vars, z, model, deadlock)
    Rmat = [cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]
    for t in 2:vars.Thor+1
        idxk1, idxk2, idxj1, idxj2 = index_z_mpc(edge.k, edge.j, t, vars.Thor)

        # Define eta based on known nominal trajectories
        delta_bar = agent.nom_traj[:,t]-neighbor.nom_traj[:,t]
        delta = [z[idxk1],z[idxk2]]-[z[idxj1],z[idxj2]]
        eta = delta_bar/norm(delta_bar)

        if deadlock == true
            eta = Rmat*eta
        end

        # Constrain linearized distance between agents to be greater than collision distance
        @constraint(model, norm(delta_bar) + eta'*(delta-delta_bar)>=vars.r_coll)

        # add constraint for non-cooperative agents to follow their nominal trajectory
        # if neighbor.cooperative == 2
        #     @constraint(model, z[idxj1] == neighbor.nom_traj[1,t])
        #     @constraint(model, z[idxj2] == neighbor.nom_traj[2,t])
        # end
    end
    return model
end

function compute_z_mpc(all_agents,vars)
    for agent in all_agents
        for edge in agent.neighborhood
            neighbor = all_agents[edge.j];
            coll = check_collision(edge, vars) 
            coll = 0
            deadlock = detect_head_on(all_agents, edge, vars.Thor)
            if coll != 0 
                model = Model(SCS.Optimizer)
                @variable(model, z[1:size(edge.yk,1)]>=-1000)

                @objective(model, Min, edge.lambdak'*(edge.yk-z)+edge.lambdaj'*(edge.yj-z) +vars.rho/2*((edge.yk - z)'*(edge.yk - z)+(edge.yj - z)'*(edge.yj - z)))

                model = collision_constraint(edge, agent, neighbor, vars, z, model, deadlock)
                
                # # solve for z
                set_silent(model)
                optimize!(model)
                edge.z = vec(value.(z));
            else
                edge.z = 0.5*(edge.yk+edge.yj);
            end
            # compute lambda
            edge.lambdak = edge.lambdak +vars.rho * (edge.yk - edge.z)
            edge.lambdaj = edge.lambdaj +vars.rho * (edge.yj - edge.z)
        end
    end
    return all_agents
end

function control_strat_distributed(all_agents, edges, vars)
    dimz = vars.dim*(vars.Thor+1)*2
    edges = neighborhoodmod.initialize_consensus_variables(edges, dimz)
    all_agents = neighborhoodmod.set_agent_neighborhood(all_agents, edges)
    # set nominal trajectory to be constant velocity
    # for agent in all_agents[1:vars.k_c]
    #     for t in 1:vars.Thor+1
    #         agent.nom_traj[:,t] = agent.est_pos+agent.est_vel*(t-1)
    #     end
    # end

    
    # Iterations
    for i = 1:vars.iter
        res = 0;
        for agent in all_agents
            n_neighbors = size(agent.neighborhood,1)
            model = Model(SCS.Optimizer)
            # model = Model(Mosek.Optimizer)
            @variable(model, p[1:vars.dim,1:n_neighbors+1, 1:vars.Thor+1]>=-1000)
            @variable(model, u[1:vars.dim, 1:vars.Thor]>=-1000)

            @expression(model, J3, sum(1/norm(agent.est_pos-agent.goal)* (p[:,1,t]-agent.goal)'*I(vars.dim)*(p[:,1,t]-agent.goal) for t in 2:vars.Thor+1))
            @expression(model, J4, 0.1*sum((u[:,t]-u[:,t-1])'*I(vars.dim)*(u[:,t]-u[:,t-1]) for t in 2:vars.Thor))
            @expression(model, J5, (u[:,1]-agent.vel)'*I(vars.dim)*(u[:,1]-agent.vel))


            @constraint(model, [t in 1:vars.Thor], p[:,1,t+1] == p[:,1,t] + u[:,t])
            @constraint(model,  p[:,1,1] == agent.est_pos)
            @constraint(model, [t in 1:vars.Thor], u[:,t]'*u[:,t]<=vars.umax^2)
            k=2
            for edge in agent.neighborhood
                neighbor = all_agents[edge.j];
                for t = 1:vars.Thor+1
                    delta_bar = agent.nom_traj[:,t]-neighbor.nom_traj[:,t]
                    delta = p[:,1,t]-p[:,k,t]
                    eta = (delta_bar)/norm(delta_bar)
                    @constraint(model, norm(delta_bar) + eta'*(delta-delta_bar)-vars.r_coll>=0)
                end
                # non-cooperative agents stay on nominal trajectory
                if neighbor.cooperative == 2
                    @constraint(model, [t in 1:vars.Thor+1], p[:,k,t] == neighbor.nom_traj[:,t])
                end
                k=k+1;
            end


            @expression(model, jf1, sum(edge.lambdak'*[vec(p[:,1,:]); vec(p[:,b+1,:])]+vars.rho/2*([vec(p[:,1,:]); vec(p[:,b+1,:])]-edge.z)'*([vec(p[:,1,:]); vec(p[:,b+1,:])]-agent.neighborhood[b].z) for (b,edge) in enumerate(agent.neighborhood) if edge.k < edge.j))
            @expression(model, jf2, sum(edge.lambdak'*[vec(p[:,b+1,:]); vec(p[:,1,:])]+vars.rho/2*([vec(p[:,b+1,:]); vec(p[:,1,:])]-edge.z)'*([vec(p[:,b+1,:]); vec(p[:,1,:])]-edge.z) for (b,edge) in enumerate(agent.neighborhood) if edge.k > edge.j))

            @objective(model, Min, J3+jf1+jf2)

            set_silent(model)

            optimize!(model)

            # # Save for this agent
            j=2
            for edge in agent.neighborhood
                if all_agents[edge.j].cooperative == 1
                    p2 = vec(value.(p[:,j,:])');
                else 
                    p2 = vec(all_agents[edge.j].nom_traj);
                end
                if edge.k < edge.j
                    temp = copy(edge.yk)
                    edge.yk = [vec(value.(p[:,1,:])'); p2];
                    res = res + norm(edge.yk-temp);
                else
                    temp = copy(edge.yk)
                    edge.yk = [p2; vec(value.(p[:,1,:])')];
                    res = res + norm(edge.yk-temp);
                end
            j=j+1
            end
            # Save for this agent
            agent.q = value.(p[:,1,:]);
            agent.est_vel = vec(value.(u[:,1]))
            agent.vel = vec(value.(u[:,1]))
        end
        # logtofile.log_u(all_agents, i)

        # Communicate values y to neighborhood
        all_agents = neighborhoodmod.communicate_y_neighborhood(all_agents)
        all_agents = compute_z_mpc(all_agents,vars)

        # set nominal trajectory
        for agent in all_agents[1:vars.k_c]
            agent.nom_traj = agent.q;
        end
        # print(res,'\n')
        if res < 10^-6
            print("control converged at iteration: ", i,'\n')
            break
        end
    end
    return all_agents
end

function distributed_2_step(all_agents, edges, var_set)
    all_agents_update = estimate_state_distributed(all_agents, edges, var_set)
    all_agents_new = control_strat_distributed(all_agents_update, edges, var_set)
    return all_agents_new
end

end