# Test parameter seraching

import RetroSignalModel as rs

# Parameter searching
@time rs.search_params(rs.rtgM4(); num_sim=3, distributed=false, saveall=true)
