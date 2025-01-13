module logtofile
using CSV, DataFrames
using Printf, Dates

function log_pos(agents, iter)
    r=10
    date1 = day(now())
    date2 = month(now())
    filecent = @sprintf("data/NEW_Distributed_POS.csv")
    result = []
    for agent in agents
        result= [result; vec(agent.pos)]
    end
    CSV.write(filecent, Tables.table(result'), append=true, writeheader = false)
end

function log_u(agents, iter)
    r=50
    date1 = day(now())
    date2 = month(now())
    file1 = @sprintf("data/NEW_Distributed_U.csv")
    result = [iter]
    for agent in agents
        result= [result; vec(agent.vel)]
    end
    CSV.write(file1, Tables.table(result'), append=true, writeheader = false)
end

function log_est_pos(agents, iter)
    r=50
    date1 = day(now())
    date2 = month(now())
    file2 = @sprintf("data/NEW_Distributed_ESTPOS.csv")
    result = [iter]
    for agent in agents
        result= [result; vec(agent.q[:,1])]
        # result = [result; vec(agent.est_pos)]
    end
    CSV.write(file2, Tables.table(result'), append=true, writeheader = false)
end

function log_path(agents, iter)
    date1 = day(now())
    date2 = month(now())
    file1 = @sprintf("data/NEW_path_%02d_%02d.csv", date2, date1)
    result = []
    for agent in agents
        result= [result; vec(agent.q)]
    end
    CSV.write(file1, Tables.table(result'), append=true, writeheader = false)
end

end