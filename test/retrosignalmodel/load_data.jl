using DataFrames
import CSV
import RetroSignalModel as rs

@testset "Check file availability" for fn in rs.DataFiles
    @test isfile(fn)
    df = DataFrame(CSV.File(fn))
end

exp = rs.get_expression_levels()
cond = rs.get_conditions()

prNamesRNA = [:Rtg1, :Rtg2, :Rtg3, :Bmh, :Mks]
prNamesCond = ["Rtg1", "Rtg2", "Rtg3", "Mks"]


@testset "Protein name availability (RNA-Seq)" for pr in prNamesRNA
    @test pr in collect(keys(exp))
end
@testset "Protein name availability (Conditions)" for pr in prNamesCond
    @test pr in names(cond)
end

@test "s" in names(cond)
@test "gfp" in names(cond)
@test "Trans2Nuc" in names(cond)


# Solutions of rtgM4
@testset "Solutions to rtgM4" begin
    df = rs.read_csv(rs.DataFiles.solution_rtgM4)
    num_sols = size(df)[1]
    for i in 1:num_sols
        @test rs.is_valid(rs.rtgM4(i))
    end
end
