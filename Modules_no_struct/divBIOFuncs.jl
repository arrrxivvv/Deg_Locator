
function fileNameAttrFunc( N, param_divide, num_it, seedFedStart; dim = nothing )
	fileNameAttr = string( "N_", N, "_param_divide_", param_divide, "_instanceNum_", num_it, "_seed_", seedFedStart );
	if !isnothing(dim)
		fileNameAttr = string( "dim_", dim, "_", fileNameAttr );
	end
	return fileNameAttr;
end

# function fileNameAttrFunc( N, param_divide, num_it, seedFedStart )
	# fileNameAttr = string( "N_", N, "_param_divide_", param_divide, "_instanceNum_", num_it, "_seed_", seedFedStart );
	# return fileNameAttr;
# end
