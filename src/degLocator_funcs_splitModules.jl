# using Infiltrator
function locateDiv( non0Arr, degBerrys::DegBerrys )
	Hlst = DegLocatorDiv.HlstFunc(H_GUE,degBerrys.params.nDim,degBerrys.params.N);
	HmatFun = (H,xLst) -> Hmat_3comb_ratio!( H, xLst, Hlst );
	
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
	
	# thres = 1e-6;
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
	
	return NLstPol[1], NLstPol[2], locLstPol[1], locLstPol[2], Hlst;
end

function findNon0Locs( non0Arr, degBerrys::DegBerrys, thres )
	# posNLst = zeros(Int64, degBerrys.params.N);
	# negNLst = zeros(Int64, degBerrys.params.N);
	isDegArr = [[zeros(Bool,degBerrys.params.divLst...)
		for n = 1:degBerrys.params.N]
		for iPol = 1:2];
	NLstPol = [ zeros(Int64, degBerrys.params.N) for iPol = 1:2 ];
	
	# isNon0 = false;
	# thresN = -thres;
	# non0Tmp = threaded_zeros(degBerrys.params.N);
	# NLstPolThr = threaded_zeros( Int64, degBerrys.params.N, 2 );
	# @time begin
	# Threads.@threads for pos in degBerrys.params.posLst
		# getThrInst( non0Tmp ) .= non0Arr[pos];
		# for iPol = 1:2, n = 1:degBerrys.params.N
			# isNon0 = (-1)^(iPol-1) * 
				# non0Arr[pos][n] > thres;
			# isNon0 = iPol == 1 ? 
				# non0Arr[pos,n] > thres : 
				# non0Arr[pos,n] < thresN;
			# if isNon0
				# NLstPolThr[n,iPol] += 1;
			# end
			# end
		# end
	# end
	# end
	# @infiltrate
	# for iPol = 1:2, n = 1:degBerrys.params.N
		# NLstPol[iPol][n] = 0;
		# for iTh = 1 : Threads.nthreads()
			# NLstPol[iPol][n] += NLstPolThr.data[iTh][n,iPol];
		# end
	# end
	
	# locLstPol = [[
		# ones(Int64, NLstPol[iPol][n], degBerrys.params.nDim)
		# for n = 1:degBerrys.params.N]
		# for iPol = 1:2];
	locLstPol = [
		Vector{Array{Int64,2}}(undef,degBerrys.params.N)
		for iPol = 1:2];
	
	# @time begin
	for iPol = 1:2, n = 1:degBerrys.params.N
		locLstPol[iPol][n] = 
			cartIndLstToArr( 
				findall( (x->(-1)^(iPol-1)*x>thres), 
					selectdim(non0Arr,degBerrys.params.nDim+1,n)), 
				degBerrys.params.nDim );
		NLstPol[iPol][n] = 
			size( locLstPol[iPol][n],1 );
	end
	# end
	# @time begin
	# for iPol = 1:2, n = 1:degBerrys.params.N
		# cnt = 1;
		# for pos in degBerrys.params.posLst
			# if (-1)^(iPol-1) * 
				# non0Arr[pos,n] > thres
				# for iDim = 1 : degBerrys.params.nDim
					# locLstPol[iPol][n][cnt,iDim] = pos[iDim];
				# end
				# cnt += 1;
			# end
		# end
	# end
	# end
	# @infiltrate
	
	return NLstPol, locLstPol;
end
