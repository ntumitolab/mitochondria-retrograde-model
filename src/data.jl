function getExpLevels(;df=DataTables.RNAseq, condition="unstressed", Max_expr=Max_expr,idxName=1)
    prNames = df[!,idxName]
    ex = df[!, condition]
    ex = capVec(ex, Max_expr)
    return NamedTuple{tuple(Symbol.(prNames)...)}(tuple(ex...))
end


function getConditions(;df=DataTables.BoolCond, drop_missing=true) 
    if drop_missing 
        return dropmissing(df)
    else
        return df
    end
end

