using Pkg
Pkg.activate("DistributedMPCGP")
# Pkg.instantiate()

using DDAAModule
using Plots
using LinearAlgebra
using Statistics
using Random, Distributions  # hide+
using HDF5
using JLD2
# Pkg.build("DDAAModule")
r=30.0
# Helping functions
function check_collision(all_agents)
    collision = 0
    for k in 1:size(all_agents,1)
        for j in k:size(all_agents,1)
            if k!=j
                for n in 1:size(all_agents,2)
                    if norm(all_agents[k,n].pos-all_agents[j,n].pos)<7
                        # print("collision between ", k, " and ", j, " at time ", n," the distance is ",norm(all_agents[k,n].pos-all_agents[j,n].pos), "\n")
                        collision += 1
                    end
                end
            end
        end
    end
    return collision
end

mutable struct constant_variables
    r::Float64
    r_coll::Float64
    umax::Float64
    Thor::Int
    dim::Int
    Sigma1::Float64
    Distr1::Normal
    Sigma2::Float64
    Distr2::Normal
    k_c::Int
    k_nc::Int
    rho::Float64
    iter::Int
    alpha::Float64
    # noiselevels::Vector{Float64}
end

# set constant_variables
r = 30.0
r_coll = 7.0
umax = 1.5;
Thor = 10;
dim=2;
# Noise levels
Sigma1 = umax*0.2;
Distr1 = Normal(0,Sigma1^2); 
Sigma2 = 0.1;
Distr2 = Normal(0,Sigma2^2);
dimz = dim*2;
rho = 0.1;
iter = 2000;
alpha = 0.95;
k_c = 7;
k_nc = 1;
var_set = constant_variables(r, r_coll, umax, Thor, dim, Sigma1, Distr1, Sigma2, Distr2,k_c,k_nc,rho,iter,alpha);


K = k_c+k_nc


filename = "test_data_8agents_t10.jld2"
N=100

agents, results_pos, results_ctrl = run_simulation(filename,r, N, var_set)




# Save the agents and plots
save_object("agents_test_data_8agents_t10.jld2", agents)



print("the number of collisions is: ", check_collision(agents))


a = Animation()

pos_acc = zeros(K,N,2)
est_pos_acc = zeros(K,N,2)

for t in 1:N
    plt = plot()

    for agent in agents[:,t]
        pos_acc[agent.id,t,:] = agent.pos
        est_pos_acc[agent.id,t,:] = agent.est_pos
        plt = plot!(pos_acc[agent.id,1:t,1],pos_acc[agent.id,1:t,2],marker=:circle)
        plt = plot!(est_pos_acc[agent.id,1:t,1],est_pos_acc[agent.id,1:t,2],marker=:circle)
        plt = plot!(agent.nom_traj[1,:],agent.nom_traj[2,:],linestyle=:dash,marker= :diamond)
    end
    frame(a, plt)
end
	
# # gif(a,fps=3)

function circleShape(x,y,r)
    θ = LinRange(0,2π,100)
    x.+r*cos.(θ), y.+r*sin.(θ)
end

plotcolors = [:red :blue :green :orange :yellow :purple :pink :brown]

times = [1,2,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,40,45,50,90,99]
for t in times
    plt = plot(fmt = :png)
    for (k,agent) in enumerate(agents[:,t])
        plt = plot!(pos_acc[agent.id,1:t,1],pos_acc[agent.id,1:t,2], linecolor = plotcolors[k],linealpha = 0.4)
        plt = plot!(est_pos_acc[agent.id,1:t,1],est_pos_acc[agent.id,1:t,2],linecolor = plotcolors[k])
        plt = plot!(agent.nom_traj[1,:],agent.nom_traj[2,:],linestyle=:dot,linecolor = plotcolors[k],linealpha = 0.5)
        # plt = plot!(circleShape(est_pos_acc[agent.id,t,1],est_pos_acc[agent.id,t,1],50), seriestype = [:shape,], linecolor = plotcolors[k], fillalpha = 0, label = "Communication radius")
        if k == 1
            plt = plot!(circleShape(est_pos_acc[agent.id,t,1],est_pos_acc[agent.id,t,2],r), seriestype = [:shape,], linecolor = plotcolors[k],linealpha=0.1, fillalpha = 0.1, label = "Communication radius")
        end
        # Add start marker
        plt = scatter!([pos_acc[agent.id, 1, 1]], [pos_acc[agent.id, 1, 2]], markercolor = plotcolors[k],marker = :x, markersize = 2, label = "Start")

        # Add goal marker (assuming the goal is the last point in the nominal trajectory)
        plt = scatter!([agent.goal[1]], [agent.goal[2]], markercolor = plotcolors[k],marker = :star5, markersize =2, markerstrokealpha = 0.8, label = "Goal")
    

    # Final plot settings
    end

    plt = plot!(xlabel="x [m]",ylabel="y [m]",legend=false,thickness_scaling = 2, xlims = (0,100), ylims = (0,100),aspect_ratio=1)
    savefig(plt, "T10t$(t).png")
end




# plot()

# for agent in agents[:,1]
#     n = agent.id
#     plot!(pos_acc[agent.id,:,1],pos_acc[agent.id,:,2],marker=:circle, label = ["$n GT" " "])
#     plot!(est_pos_acc[agent.id,:,1],est_pos_acc[agent.id,:,2],marker=:diamond, label = ["$n Est" " "])
# end
# current()


# # Get the number of columns
# _, num_cols = size(results_pos)

# # Initialize a vector to store the averages
# average_pos = Vector{Float64}(undef, num_cols)
# average_ctrl = Vector{Float64}(undef, num_cols)

# # Compute the column-wise average, ignoring NaN values
# for j in 1:num_cols
#     column = results_pos[:, j]
#     average_pos[j] = mean(filter(!isnan,column))
#     column2 = results_ctrl[:, j]
#     average_ctrl[j] = mean(filter(!isnan,column2))
# end

# plt=plot(average_ctrl, label="Control", xlabel="Iteration", ylabel="Average", title="Average Control and Position", linewidth=2)
# plt =plot!(average_pos, label="Position", linewidth=2,thickness_scaling = 1.5)
# savefig(plt, "convergence_03.png")
