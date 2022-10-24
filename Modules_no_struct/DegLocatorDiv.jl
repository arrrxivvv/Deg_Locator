module DegLocatorDiv

using UMat
using Utils
using DivB

include("divBIOFuncs.jl")
export fileNameAttrFunc

include("degLocator_funcs.jl")
export locator_div, locator_div_GUE, locator_div_sin3, locLstDistill, locator_div_GUE_scale, locator_div_GUE_ratio, locator_div_GOE_ratio

include("divBProfile_funcs.jl")
export divB_profile

include("deltaN_funcs.jl")
export deltaN

include("deltaNAvg_funcs.jl")
export deltaN_avg, deltaN_avg_fromFile, deltaN_avg_lst, deltaN_avg_lst_fromFile

include("fFunc_stat.jl")
export fFunc_stat, fFunc_stat_from_file, fFunc_stat_diff_from_file, fFunc_stat_diff_full_from_file, fFunc_stat_diff_full_3d_from_file

end
