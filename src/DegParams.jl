struct DegParams
	divLst::Vector{Int64};
	nDim::Int64;
	N::Int64;
	nonPeriodic::Bool;
	posLst::AbstractArray{CartesianIndex{N}, N} where N;
	
	minLst::Vector{Float64};
	maxLst::Vector{Float64};
	stepLst::Vector{Float64};
	gridLst::Vector{ Vector{Float64} };
	mesh::Array{ Vector{Float64} };
end

function degParamsBase( N, divLst, minLst, maxLst, nDim; isNonPeriodic = false )
	stepLst = ( maxLst .- minLst ) ./ divLst;
	
	if isNonPeriodic
		posLst = CartesianIndices(Tuple(divLst.+1));
		gridLst = [ collect( range( minLst[iDim], maxLst[iDim], length = divLst[iDim]+1 ) ) for iDim = 1:nDim ];
	else
		posLst = CartesianIndices(Tuple(divLst));
		gridLst = [ collect( range( minLst[iDim], maxLst[iDim] - stepLst[iDim], length = divLst[iDim] ) ) for iDim = 1:nDim ];
	end
	
	mesh = [ [ gridLst[j][ind[j]] for j = 1:nDim ] for ind in posLst ];
	
	return DegParams( divLst, nDim, N, isNonPeriodic, posLst, minLst, maxLst, stepLst, gridLst, mesh );
end

function degParamsNonInit( N, divLst, nDim; isNonPeriodic = false )
	minLst = zeros(nDim);
	maxLst = zeros(nDim); 
	
	return degParamsBase( N, divLst, minLst, maxLst, nDim; isNonPeriodic = isNonPeriodic );
end

function degParamsInit( N, divLst, minLst, maxLst, nDim; isNonPeriodic = false )
	if !isa(minLst, Array)
		minLst = fill( minLst, nDim );
	end
	if !isa(maxLst, Array)
		maxLst = fill( maxLst, nDim );
	end
	
	return degParamsBase( N, divLst, minLst, maxLst, nDim; isNonPeriodic = false );
	
	# stepLst = ( maxLst .- minLst ) ./ divLst;
	
	# if isNonPeriodic
		# posLst = CartesianIndices(Tuple(divLst.+1));
		# gridLst = [ collect( range( minLst[iDim], maxLst[iDim], length = divLst[iDim]+1 ) ) for iDim = 1:nDim ];
	# else
		# posLst = CartesianIndices(Tuple(divLst));
		# gridLst = [ collect( range( minLst[iDim], maxLst[iDim] - stepLst[iDim], length = divLst[iDim] ) ) for iDim = 1:nDim ];
	# end
	
	# mesh = [ [ gridLst[j][ind[j]] for j = 1:nDim ] for ind in posLst ];
	
	# return DegParams( divLst, nDim, N, isNonPeriodic, posLst, minLst, maxLst, stepLst, gridLst, mesh );
end

function degParamsPeriodic( N, divLst, minLst, maxLst, nDim )
	isNonPeriodic = false;
	return degParamsInit( N, divLst, minLst, maxLst, nDim; isNonPeriodic = isNonPeriodic );
end

function makeArrOverGrid( type::DataType, params::DegParams )
	# return Array{type,params.nDim}(undef,params.divLst...);
	return similar( params.posLst, type );
end

function updateParamsRange( minLst, maxLst, params::DegParams )
	if !isa( minLst, Array )
		minLst = fill( minLst, params.nDim );
	end
	if !isa( maxLst, Array )
		maxLst = fill( maxLst, params.nDim );
	end
	params.minLst .= minLst;
	params.maxLst .= maxLst;
	params.stepLst .= ( maxLst .- minLst ) ./ params.divLst;
	
	for iDim = 1 : params.nDim
		for ii = 0 : length( params.gridLst[iDim] ) - 1
			params.gridLst[iDim][ii+1] = params.stepLst[iDim] * ii;
		end
	end
	for pos in params.posLst
		for iDim = 1 : params.nDim
			params.mesh[pos][iDim] = params.gridLst[iDim][pos[iDim]];
		end
	end
end

function linIdFromIdVec( idVec::Vector{Int64}, params::DegParams )
	id = idVec[end]-1;
	for iDim = (params.nDim-1):-1:1
		id = id * params.divLst[iDim] + idVec[iDim]-1;
	end
	id = id+1;
end

function wrapIdVec!( idVec::Vector{Int64}, params::DegParams )
	idVec .-= 1;
	idVec .= mod.( idVec, params.divLst );
	idVec .+= 1;
end
