function locateDiv( degBerrys::DegBerrys, non0Arr; HmatFun )
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
	# ((x,y)->(x.=real.(y))).(non0Arr, non0ArrCmplx );
	
	@info("find Non0")
	Utils.@timeInfo NLstPol, locLstPol = findNon0Locs( non0Arr, degBerrys, thresNon0 );
	@info("GC")
	Utils.@timeInfo GC.gc();
	
	# @infiltrate
	
	return NLstPol[1], NLstPol[2], locLstPol[1], locLstPol[2]; 
	# , Hlst;
end

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

function degTmpArrs( params::DegParams, enumSaveMem::EnumSaveMem )
	if enumSaveMem == memNone
		matsFull = matsGridHThreaded( params, threaded_zeros(ComplexF64,params.N,params.N) );
		matsGridInitAll( matsFull );
		degBerrysFull = degBerrysInit( params, matsFull; isFullInit = true );
	elseif enumSaveMem >= memEig
		degBerrysFull = degBerrysEigLayered( params );
	end
	non0Arr = zeros(params.divLst...,params.N);
	return degBerrysFull, non0Arr;
end
