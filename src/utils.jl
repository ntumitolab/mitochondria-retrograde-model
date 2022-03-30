using LinearAlgebra

"""Read a CSV file and return a DataFrame"""
read_csv(filename; kwargs...) = DataFrame(CSV.File(filename); kwargs...)

"""Merge csv files under a directory (first layer)."""
function mergeCSV(dirname; format=".csv")

    fns_ = readdir(dirname)

    fns = filter_files(fns_, format)

    length(fns) == 0 ? @warn("No file is found.") : nothing

    fns = joinpath.(dirname, fns)

    df = CSV.read(fns[1], DataFrame)
    for fn in fns[2:end]
        df_ = CSV.read(fn, DataFrame)
        df = vcat(df, df_)
    end

    return df
end


function catalyst_name(m::ReactionSystem; fcall=Catalyst.species, remove_t=true)
    spec_names_t = string.(fcall(m))
    spec_names = remove_t ? replace.(spec_names_t, "(t)" => "") : spec_names_t
    return spec_names
end

function getTables(csv_fns::NamedTuple)
    dfs = []
    for fn in csv_fns
        push!(dfs, read_csv(fn))
    end
    return NamedTuple{keys(csv_fns)}(tuple(dfs...))
end

"""
Get constructor from a data type
"""
get_ctor(datatype) = typeof(datatype).name.wrapper


function find_id(sNames, features; convert_lower=true)
    sNames = convert_lower ? lowercase.(sNames) : sNames
    features = convert_lower ? lowercase.(features) : features
    idxs = Vector{Int}()
    for (ind, i) in enumerate(sNames)
        isEle = true
        for f in features
            if !occursin(f, string(i))
                isEle = false
                break
            end
        end
        if isEle
            push!(idxs, ind)
        end
    end
    return idxs
end

"""
Scale the maximum element of vector `v` to `capVal`.
"""
cap_vec(v::AbstractVector, capVal) = v .* (capVal / norm(v, Inf))
