# Parameter searching api

"""
Parameter searching Api. It can be run in single core or multiple-processing

Performance
------------
```
julia> addprocs(exeflags="--project=$(Base.active_project())");
julia> @time @everywhere using RetroSignalModel
julia> nprocs()
9
julia> @time rs.paramSearching(rs.rtgM4(); num_sim=100, distributed=true, saveall=true)
  2.378919 seconds (78.95 k allocations: 4.964 MiB)
```
"""
function paramSearching(m::RTGmodel;num_sim=10, save_iter=num_sim, distributed=false, saveall=true,batch_size=nprocs(),conditions=getConditions(;df=DataTables.BoolCond), kwags...)
    #create table 
    df = get_search_df(m, conditions)

    ran = 1:save_iter:num_sim
    iter = ran.step
    try
        for _ in 1:length(ran)       
            # parameter searching 
            res =  @showprogress  pmap( x->paramSearching_Once(x[1];conditions=x[2], x[3]...), 
                                    fill([m, conditions, kwags], iter),
                                    distributed=distributed,
                                    batch_size=batch_size)

            # Merge result
            map(r-> r.isValid || saveall ? push_sim!(df, r.test_log, r.p, r.u) : nothing, res)
        end
        return df
    catch e 
        return df
    end

    
end


"""
Sampling one parameter set for validation
"""
function paramSearching_Once(m::RTGmodel;Break=false, idx_hill_coefs=[1], K_dist=K_dist, K_N_dist=K_N_dist,conditions=getConditions(;df=DataTables.BoolCond),show_valid=true,ssmethod=SSMETHOD,kwags...)
    # Sampling parameters 
    p = init_p(m.model;idx_hill_coefs=idx_hill_coefs, K_dist=K_dist, K_N_dist=K_N_dist)
    m.p[:] = p 

    u = getSteadySol(m, m.u; ssmethod=ssmethod).u

    # Test validity
    isValid, test_log = try_conditions(m;conditions=conditions, Break=Break,kwags...)

    # Print valid result
    show_valid && isValid ? @show(p) : nothing

    return (; isValid=isValid, test_log=test_log, u=u, p=p)
end


# Data processing
function get_search_df(m, conditions)
    # Labels
    condNs = ["condition $i" for i in range(1,stop= size(conditions)[1])]
    pNs = string.(keys(init_p(m.model)))
    uNames = string.(keys(init_u(m)))
    # Create empty dataframe
    df_key = [condNs..., pNs...,uNames...]
    trials = DataFrame( [ [] for i in eachindex(df_key) ], df_key)

    return trials
end

function push_sim!(df, conds, ps, us)
    vec = [conds..., ps...,us...]
    push!(df, vec)
end