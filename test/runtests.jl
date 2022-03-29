using Test

@testset "RetroSignalModel" begin

    @testset "Parameter searching" begin
        include("retrosignalmodel/searching.jl")
    end

    @testset "Data files" begin
        include("retrosignalmodel/load_data.jl")
    end

    @testset "Build models" begin
        include("retrosignalmodel/build_model.jl")
    end

end # @testset "RetroSignalModel"


