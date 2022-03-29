"""
Test multi-processing
"""
# Load library to main worker
@time using Distributed, RetroSignalModel, CSV
# Add workers
addprocs(1, exeflags="--project=$(Base.active_project())");
@show nprocs();

# Load library to workers
@time @everywhere import RetroSignalModel as rs

# Parameter searching
@time rs.paramSearching(rs.rtgM4(); num_sim=3, distributed=true, saveall=true)
