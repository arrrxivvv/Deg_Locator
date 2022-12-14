using ShiftedArrays
using ThreadedArrays

struct DegParams
	divLst::Vector{Int64};
	nDim::Int64;
	N::Int64;
	nonPeriodicLst::Vector{Bool};
	posLst::AbstractArray{CartesianIndex{N}, N} where N;
	
	minLst::Vector{Float64};
	maxLst::Vector{Float64};
	stepLst::Vector{Float64};
	gridLst::Vector{ Vector{Float64} };
	mesh::Array{ Vector{Float64} };
	
	locItThr::ThrArray{Int64,1};
	linIdThr::ThrStruct{Int64};
	stepsItThr::ThrArray{Int64,1};
	divsItThr::ThrArray{Int64,1};
	
	posLstSh::Array{CircShiftedArray};
end

function degParamsBase( N, divLst, minLst, maxLst, nDim; isNonPeriodic = false )
	stepLst = ( maxLst .- minLst ) ./ divLst;
	
	if isa(isNonPeriodic, Bool)
		isNonPeriodic = fill(isNonPeriodic,nDim);
	end
	
	szLst = copy(divLst);
	maxGridLst = copy(maxLst);
	szLst .+= isNonPeriodic;
	maxGridLst .-= .!isNonPeriodic .* stepLst;
	
	posLst = CartesianIndices( Tuple(szLst) );
	gridLst = [ collect( range( minLst[iDim], maxGridLst[iDim], length = szLst[iDim] ) ) for iDim = 1:nDim ];
	
	mesh = [ [ gridLst[j][ind[j]] for j = 1:nDim ] for ind in posLst ];
	
	locItThr = threaded_zeros( Int64, nDim );
	linIdThr = thrStructFill( Int64, 0 );
	stepsItThr = threaded_zeros( Int64, nDim );
	divsItThr = threaded_zeros( Int64, nDim );
	
	posLstSh = [ 
		ShiftedArrays.circshift( 
			posLst, 
			ntuple(( i -> i==iDim ? (-1)^iSgn : 0 ), nDim ) )
		for iDim = 1 : nDim, iSgn = 1 : 2
		 ];
	
	return DegParams( deepcopy(divLst), nDim, N, isNonPeriodic, posLst, deepcopy(minLst), deepcopy(maxLst), stepLst, gridLst, mesh,  locItThr, linIdThr, stepsItThr, divsItThr, posLstSh );
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
	
	return degParamsBase( N, divLst, minLst, maxLst, nDim; isNonPeriodic = isNonPeriodic );
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
		for ii = 1 : length( params.gridLst[iDim] )
			params.gridLst[iDim][ii] = params.stepLst[iDim] * (ii-1) + minLst[iDim];
		end
	end
	for pos in params.posLst
		for iDim = 1 : params.nDim
			params.mesh[pos][iDim] = params.gridLst[iDim][pos[iDim]];
		end
	end
end

function updateParamsRangeSteps( minLst, stepLst, params::DegParams )
	if !isa( minLst, Array )
		minLst = fill( minLst, params.nDim );
	end
	params.minLst .= minLst
	params.stepLst .= stepLst;
	params.maxLst .= minLst .+ params.divLst .* stepLst;
	
	refreshMeshGrid( params );
end

function updateParamsMin( minLst, params::DegParams )
	params.minLst .= minLst;
	params.maxLst .= minLst .+ params.divLst .* stepLst;
	refreshMeshGrid( params );
end

function refreshMeshGrid( params::DegParams )
	for iDim = 1 : params.nDim
		for ii = 1 : length( params.gridLst[iDim] )
			params.gridLst[iDim][ii] = params.stepLst[iDim] * (ii-1) + params.minLst[iDim];
		end
	end
	for pos in params.posLst
		for iDim = 1 : params.nDim
			params.mesh[pos][iDim] = params.gridLst[iDim][pos[iDim]];
		end
	end
end

function linIdFromIdVecArr( idVec, arr::Array )
	id = idVec[ndims(arr)]-1;
	for iDim = (ndims(arr)-1):-1:1
		id = id * size(arr,iDim) + idVec[iDim]-1;
	end
	id = id+1;
	return id;
end

function linIdFromIdVec( idVec, params::DegParams )
	id = idVec[end]-1;
	for iDim = (params.nDim-1):-1:1
		id = id * size(params.posLst,iDim) + idVec[iDim]-1;
	end
	# iDim = params.nDim-1;
	# while iDim > 0
		# id = id * size(params.posLst,iDim) + idVec[iDim]-1;
		# iDim -= 1;
	# end
	id = id+1;
	return id;
end

function getCurrentLinId( params::DegParams )
	params.linIdThr.data[Threads.threadid()] = params.locItThr[end];
	params.linIdThr.data[Threads.threadid()] -= 1;
	# for iDim = (params.nDim-1):-1:1
		# id = id * size(params.posLst,iDim) + idVec[iDim]-1;
	# end
	iDim = params.nDim-1;
	while iDim > 0
		params.linIdThr.data[Threads.threadid()] *= size(params.posLst,iDim);
		params.linIdThr.data[Threads.threadid()] += params.locItThr[iDim]-1;
		iDim -= 1;
	end
	params.linIdThr.data[Threads.threadid()] += 1;
end

function wrapIdVecArr!( idVec::Vector{Int64}, arr::Array )
	idVec .-= 1;
	idVec .= mod.( idVec, size(arr) );
	idVec .+= 1;
end

function wrapIdVec!( idVec, params::DegParams )
	for iDim = 1 : params.nDim
		if !params.nonPeriodicLst[iDim]
			idVec[iDim] -= 1;
			idVec[iDim] = mod( idVec[iDim], params.divLst[iDim] );
			idVec[iDim] += 1;
		end
	end
end

function wrapCurrentIdVec( params::DegParams )
	wrapIdVec!( getThrInst( params.locItThr ), params );
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
	return linIdFromIdVec( getThrInst( params.locItThr ), params );
end

# function setCurrentLoc( params::DegParams, loc )
	# if isa(loc, CartesianIndex)
		# for iDim = 1 : params.nDim
			# params.locItThr[iDim] = loc[iDim];
		# end
	# elseif isa( loc, Vector ) || isa( loc, SubArray )
		# broadcastAssign!( params.locItThr, loc );
	# end
# end

function setCurrentLoc( params::DegParams, loc )
	for iDim = 1 : params.nDim
		params.locItThr[iDim] = loc[iDim];
	end
	
end

function setDoubleLoc( params::DegParams )
	for iDim = 1 : params.nDim
		params.locItThr[iDim] = 2 *params.locItThr[iDim] - 1;
	end
end

function isPosSurface( pos, params::DegParams )
	isSurface = false;
	for iDim = 1 : params.nDim
		isSurface = ( isSurface 
			|| ( pos[iDim] == 1 )
			|| ( pos[iDim] == size(params.posLst,iDim) ) );
		if isSurface
			break;
		end
	end
	return isSurface;
end

function getArrShifted( arr::Array, pos, iDim, iSh, params::DegParams )
	if isa( pos, CartesianIndex )
		idSh = iSh > 0 ? 1 : 2;
		return arr[ params.posLstSh[iDim, idSh][pos] ];
	elseif isa( pos, Vector )
		setCurrentLoc( params, pos );
		params.locItThr[iDim] += iSh;
		wrapIdVecArr( getThrInst(params.locItThr), arr );
		linId = linIdFromIdVecArr( getThrInst(params.locItThr), arr );
		params.locItThr[iDim] -= iSh;
		wrapIdVecArr( getThrInst(params.locItThr), arr );
		return arr[linId];
	end
end

function getCurrentArrSh( arr::Array, params::DegParams; dimSh = 0, iSh = 0 )
	if dimSh != 0
		params.locItThr[dimSh] += iSh;
		wrapIdVecArr!( getThrInst(params.locItThr), arr );
		linId = linIdFromIdVecArr( getThrInst(params.locItThr), arr );
		params.locItThr[dimSh] -= iSh;
		wrapIdVecArr!( getThrInst(params.locItThr), arr );
	else
		@time linId = linIdFromIdVecArr( 
		getThrInst( params.locItThr ), arr );
		# @time linId = linIdFromIdVec( 
		# getThrInst( params.locItThr ), params );
		# getCurrentLinId( params );
	end
	# return arr[getThrInst( params.linIdThr )];
	return arr[linId];
end

function getCurrentArr( arr::Array, params::DegParams )
	return getCurrentArrSh( arr, params );
end

function needInitArr( arr::Array, pos )
	if !isa(pos, Int64)
		linId  = linIdFromIdVecArr( pos, arr );
	end
	return isassigned( arr, linId ), linId;
end
