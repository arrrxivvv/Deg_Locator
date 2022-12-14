function locateDiv( degBerrys::DegBerrys, non0Arr; HmatFun )
	posNLst, negNLst, posLocLst, negLocLst, BfieldLst, divBLst = locateDiv_detailedOutput( degBerrys, non0Arr; HmatFun = HmatFun );
	return posNLst, negNLst, posLocLst, negLocLst;
end

function locateDiv_detailedOutput( degBerrys::DegBerrys, non0Arr; HmatFun )
	thresNon0 = 1e-6;
	
	if degBerrys.enumSaveMem >= memEig
		@info("Eigen and Link layered:")
		Utils.@timeInfo linksCalcAllLayered( degBerrys, HmatFun );
	else
		@info("Eigen:")
		startNextEigen( degBerrys.degMats );
		Utils.@timeInfo eigenAll( degBerrys.degMats; HmatFun = HmatFun );

		@info("Link:")
		Utils.@timeInfo linksCalcAll( degBerrys );
	end
	
	@info("Bfield:")
	Utils.@timeInfo BfieldCalcAll( degBerrys );
	
	if degBerrys.params.nDim >= 3
		@info("DivB:")
		Utils.@timeInfo divBCalcAll( degBerrys );
		non0ArrCmplx = degBerrys.divBLst;
	elseif degBerrys.params.nDim == 2
		non0ArrCmplx = degBerrys.BfieldLst[1];
	end
	Threads.@threads for pos in degBerrys.params.posLst
		for n = 1:degBerrys.params.N
			non0Arr[pos,n] = real( non0ArrCmplx[pos][n] );
		end
	end
	
	@info("find Non0")
	Utils.@timeInfo NLstPol, locLstPol = findNon0Locs( non0Arr, degBerrys, thresNon0 );
	@info("GC")
	Utils.@timeInfo GC.gc();
	
	if degBerrys.params.nDim == 2
		NLstPol[1] .+= NLstPol[2];
		locLstPol[1] .= vcat.( locLstPol[1], locLstPol[2] );
	end
	
	return NLstPol[1], NLstPol[2], locLstPol[1], locLstPol[2], degBerrys.BfieldLst, degBerrys.divBLst; 
end

# function locateDiv( degBerrys::DegBerrys, non0Arr; HmatFun )
	# thresNon0 = 1e-6;
	
	# if degBerrys.enumSaveMem >= memEig
		# @info("Eigen and Link layered:")
		# Utils.@timeInfo linksCalcAllLayered( degBerrys, HmatFun );
	# else
		# @info("Eigen:")
		# startNextEigen( degBerrys.degMats );
		# Utils.@timeInfo eigenAll( degBerrys.degMats; HmatFun = HmatFun );

		# @info("Link:")
		# Utils.@timeInfo linksCalcAll( degBerrys );
	# end
	
	# @info("Bfield:")
	# Utils.@timeInfo BfieldCalcAll( degBerrys );
	
	# if degBerrys.params.nDim >= 3
		# @info("DivB:")
		# Utils.@timeInfo divBCalcAll( degBerrys );
		# non0ArrCmplx = degBerrys.divBLst;
	# elseif degBerrys.params.nDim == 2
		# non0ArrCmplx = degBerrys.BfieldLst[1];
	# end
	# Threads.@threads for pos in degBerrys.params.posLst
		# for n = 1:degBerrys.params.N
			# non0Arr[pos,n] = real( non0ArrCmplx[pos][n] );
		# end
	# end
	
	# @info("find Non0")
	# Utils.@timeInfo NLstPol, locLstPol = findNon0Locs( non0Arr, degBerrys, thresNon0 );
	# @info("GC")
	# Utils.@timeInfo GC.gc();
	
	# if degBerrys.params.nDim == 2
		# NLstPol[1] .+= NLstPol[2];
		# locLstPol[1] .= vcat.( locLstPol[1], locLstPol[2] );
	# end
	
	# return NLstPol[1], NLstPol[2], locLstPol[1], locLstPol[2]; 
# end

function findNon0Locs( non0Arr, degBerrys::DegBerrys, thres )
	isDegArr = [[zeros(Bool,degBerrys.params.divLst...)
		for n = 1:degBerrys.params.N]
		for iPol = 1:2];
	NLstPol = [ zeros(Int64, degBerrys.params.N) for iPol = 1:2 ];
	
	locLstPol = [
		Vector{Array{Int64,2}}(undef,degBerrys.params.N)
		for iPol = 1:2];
	
	for iPol = 1:2, n = 1:degBerrys.params.N
		locLstPol[iPol][n] = 
			cartIndLstToArr( 
				findall( (x->(-1)^(iPol-1)*x>thres), 
					selectdim(non0Arr,degBerrys.params.nDim+1,n)), 
				degBerrys.params.nDim );
		NLstPol[iPol][n] = 
			size( locLstPol[iPol][n],1 );
	end
	
	return NLstPol, locLstPol;
end

function non0ArrInit( params::DegParams )
	return zeros(params.divLst..., params.N);
end

function degTmpArrs( params::DegParams, enumSaveMem::EnumSaveMem )
	typeHElem = ComplexF64;
	# if params.nDim == 2
		# typeHElem = Float64;
	# end
	if enumSaveMem == memNone
		matsFull = matsGridHThreaded( params, threaded_zeros(typeHElem,params.N,params.N); typeElm = typeHElem );
		matsGridInitAll( matsFull );
		degBerrysFull = degBerrysInit( params, matsFull; isFullInit = true );
	elseif enumSaveMem >= memEig
		degBerrysFull = degBerrysEigLayered( params; typeElm = typeHElem );
	end
	non0Arr = non0ArrInit( params );
	return degBerrysFull, non0Arr;
end

function locateDivCell( degBerrys::DegBerrys, non0Arr, gapLn; HmatFun )
	posNLst, negNLst, posLocLst, negLocLst, BfieldLst, divBLst = locateDivCell_detailedOutput( degBerrys, non0Arr, gapLn; HmatFun );
	
	return posNLst, negNLst, posLocLst, negLocLst;
end

function locateDivCell_detailedOutput( degBerrys::DegBerrys, non0Arr, gapLn; HmatFun )
	@info("divB celled")
	Utils.@timeInfo divBCellLayered( degBerrys, HmatFun, gapLn );
	Threads.@threads for pos in degBerrys.params.posLst
		for n = 1 : degBerrys.params.N
			non0Arr[pos,n] = real( degBerrys.divBLst[pos][n] );
		end
	end
	
	@info("find Non0")
	thresNon0 = 1e-6;
	Utils.@timeInfo NLstPol, locLstPol = findNon0Locs( non0Arr, degBerrys, thresNon0 );
	@info("GC")
	Utils.@timeInfo GC.gc();
	return NLstPol[1], NLstPol[2], locLstPol[1], locLstPol[2], degBerrys.BfieldLstSurface, degBerrys.divBLst; 
end

# function locateDivCell( degBerrys::DegBerrys, non0Arr, gapLn; HmatFun )
	# @info("divB celled")
	# Utils.@timeInfo divBCellLayered( degBerrys, HmatFun, gapLn );
	# Threads.@threads for pos in degBerrys.params.posLst
		# for n = 1 : degBerrys.params.N
			# non0Arr[pos,n] = real( degBerrys.divBLst[pos][n] );
		# end
	# end
	
	# @info("find Non0")
	# thresNon0 = 1e-6;
	# Utils.@timeInfo NLstPol, locLstPol = findNon0Locs( non0Arr, degBerrys, thresNon0 );
	# @info("GC")
	# Utils.@timeInfo GC.gc();
	# return NLstPol[1], NLstPol[2], locLstPol[1], locLstPol[2]; 
# end

function degTmpArrsCell( params::DegParams, gapLn::Int64 )
	degBerrysCell = degBerrysFineSurfaceInit( params, gapLn );
	non0Arr = non0ArrInit( params );
	
	return degBerrysCell, non0Arr, gapLn;
end
