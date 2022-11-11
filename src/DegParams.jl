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
	
	locItThr::ThrArray{Int64,1};
	stepsItThr::ThrArray{Int64,1};
	divsItThr::ThrArray{Int64,1};
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
	
	locItThr = threaded_zeros( Int64, nDim );
	stepsItThr = threaded_zeros( Int64, nDim );
	divsItThr = threaded_zeros( Int64, nDim );
	
	return DegParams( divLst, nDim, N, isNonPeriodic, posLst, minLst, maxLst, stepLst, gridLst, mesh, locItThr, stepsItThr, divsItThr );
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

function shLinId!( idLin, iSh, params::DegParams )
	iSh -= 1;
	for iDim = 1 : params.nDim
		idVec[iDim] = (iSh & 1);
		iSh = iSh >> 1;
	end
end

function shIdVec!( params::DegParams, loc, iSh; shD = 1 )
	broadcastAssign!( params.locItThr, loc );
	
	shLocIt!( params, iSh; shD=shD );
end

function shLocIt!( params::DegParams, iSh; shD = 1 )
	iSh -= 1;	
	for iDim = 1:params.nDim
		params.locItThr[iDim] += (iSh & 1) * 	( isa(shD,Array) ? shD[iDim] : shD );
		iSh = iSh >> 1;
	end
	
	wrapIdVec!( getThrInst( params.locItThr ), params )
end

function locItCorner!( params::DegParams, iSh )
	broadcastAssign!( params.locItThr, 1 );
	
	shLocIt!( params, iSh; shD = params.divLst );
end

function linItCurrent( params::DegParams )
	return linIdFromIdVec( getThrInst( params.itLocThr ), params );
end
