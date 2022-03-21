using Pkg 
Pkg.activate(@__DIR__)
cd(@__DIR__)
#Pkg.instantiate()

using Weave

function weave_md(filename) 
    @info filename
    Weave.weave(filename; doctype = "github", fig_ext= ".png")
    @info "Done"
end

weave_md("find_valid_solutions.jl")
