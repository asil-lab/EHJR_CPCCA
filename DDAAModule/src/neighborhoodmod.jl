module neighborhoodmod

# mutable struct Results
#     positions::Vector{Float64}
#     controls::Vector{Float64}
#     function Results(positions::Matrix{Float64}, controls::Matrix{Float64})
#         return new(positions, controls)
#     end
# end

function set_agent_neighborhood(all_agents,edges)
    for agent in all_agents
        # Find neighbors
        agent.neighborhood = [];
        for edge in edges
            if edge.k == agent.id
                agent.neighborhood = vcat(agent.neighborhood, edge)
            end
        end
    end
    return all_agents
end

function initialize_consensus_variables(edges, dim)
    for edge in edges
        edge.yk = zeros(dim,)
        edge.yj = zeros(dim,)
        edge.z = zeros(dim,)
        edge.lambdak = zeros(dim,)
        edge.lambdaj = zeros(dim,)
    end
    return edges
end

function communicate_y_neighborhood(all_agents)
    for agent in all_agents
        for edge in agent.neighborhood
            if all_agents[edge.j].cooperative ==2
                edge.yj = edge.yk
            else
                for nb in all_agents[edge.j].neighborhood
                    if nb.j == edge.k
                        nb.yj = edge.yk
                    end
                end
            end
        end
    end
    return all_agents
end

function set_converged_values(edges)
    for a in eachindex(all_agents[1:n_cooperative]) 
        all_agents[a].est_pos =  all_agents[a].q 
    end
    # set values for non-cooperative agents
    for edge in edges
        if all_agents[edge.j].cooperative == 2
            idx = dim+dim*Thor+1
            all_agents[edge.j].est_pos = edge.yk[idx:idx+1]
        end
    end
end

end