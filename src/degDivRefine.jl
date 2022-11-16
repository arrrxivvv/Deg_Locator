using Infiltrator

function degDivRefineFromFile( mSz, divLst, itNum, seed; fMod = "", fOutMod = "", dim = 3, fExt = jld2Type, numRes = 3, thres = 0.1 )
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim=dim );
	fMain = fDeg;
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	posLocLst = load( fName, varNamePosLoc );
	negLocLst = load( fName, varNameNegLoc );
	HLstLst = load( fName, varNameHlst );
	
	locLstPol = [posLocLst,negLocLst];
	
	param_min_num = 0;
	param_max_num = 2*pi;
	params = degParamsPeriodic( mSz, divLst, param_min_num, param_max_num, dim );
	
	HmatThr = threaded_zeros( ComplexF64, params.N, params.N );
	degMatsGrid = matsGridHThreaded( params, HmatThr );
	
	# numRes = 3;
	
	divBLst = [ [ 
		[ zeros( size(locLstPol[iPol][it][lev],1), numRes ) for lev = 1:mSz ] 
			for it = 1:itNum ] 
				for iPol = 1:2 ];
	BfieldLst = [
		[ [ zeros( 2, 3, size( locLstPol[iPol][it][n], 1 ), numRes )
			for n = 1 : mSz ]
				for it = 1 : itNum ] 
					for iPol = 1 : 2];
	
	
	divLstCube = fill(2,params.nDim);
	
	matsGridCubeLst = Vector{DegMatsOnGrid}(undef, numRes);
	degBerrysLst = Vector{DegBerrys}(undef,numRes);
	
	divLstCube = ones(Int64,params.nDim);
	for iRes = 1 : numRes
		divLstCube .= 2^(iRes-1);
		matsGridCubeLst[iRes] = 
			matsGridHThreaded( 
				degParamsNonInit( params.N, divLstCube, params.nDim; isNonPeriodic = true ), 
				HmatThr
			);
		degBerrysLst[iRes] = 
			degBerrysInit( 
				matsGridCubeLst[iRes].params, 
				matsGridCubeLst[iRes] );
	end
	# @infiltrate
	
	numSh = 2^params.nDim;
	maxLstTmp = zeros(params.nDim);
	iShTmp = zeros(Int64,params.nDim);
	for it = 1 : itNum
		# println( "it = ", it, "/", itNum );
		HLst = HLstLst[it];
		HmatFun = ( (H,xLst) -> Hmat_3comb_ratio!(H,xLst,HLst) );
		# HmatGridFun reevaluate
		HmatGridFun = (xLst) -> HmatFun( getHLoc( iLst, degMatsGrid ), 
			degMatsGrid.mesh[
				linIdFromIdVec( iLst, params )] 
		);
		startNextEigen( degMatsGrid );
		for n = 1 : mSz, iPol = 1 : 2
			print( "it = ", it, " / ", itNum, ", ", "n = ", n, " / ", mSz, "              \r" );
			# println( "n = ", n, );
			locLst = locLstPol[iPol][it][n];
			# Threads.@threads 
			@info( "re eigen corners" )
			tFull = @timed begin
			for iLoc = 1 : size(locLst,1)
				loc = @view(locLst[iLoc,:]);
				
				for iSh = 1 : numSh
					setCurrentLoc( degMatsGrid.params, loc );
					shLocIt!( degMatsGrid.params, iSh );
					
					eigenAtCurrentLoc( degMatsGrid; HmatFun = HmatFun );
				end
			end
			end
			@info( timeMemStr( tFull.time, tFull.bytes ) )
			@info( "finer eigen each cell" )
			tFull = @timed begin
			# Threads.@threads 
			for iLoc = 1 : size(locLst, 1)
				loc = @view(locLst[iLoc,:]);
				linId = linIdFromIdVec( loc, params );
				maxLstTmp .= params.mesh[linId] .+ params.stepLst;
				for iRes = 1 : numRes
					matsCube = matsGridCubeLst[iRes];
					updateParamsRange( params.mesh[linId], maxLstTmp, matsCube.params );
					
					startNextEigen( matsCube );
					if iRes == 1
						matsGridTransfer!( matsCube, degMatsGrid; locS = loc );
					else
						matsGridTransferSurfaceDouble!(matsCube, matsGridCubeLst[iRes-1]);
					end
					
					divBSurfaceOutput( degBerrysLst[iRes], HmatFun; transferredResults = true );
					divBLst[iPol][it][n][iLoc,iRes] = 
					real(degBerrysLst[iRes].divBSurface[n]);
					for iDim = 1 : params.nDim, iBnd = 1 : 2
						BfieldLst[iPol][it][n][iBnd,iDim,iLoc,iRes] = real( degBerrysLst[iRes].BfieldLstSurface[iBnd,iDim,n] );
					end
				end
				# @infiltrate
			end
			end
			@info( timeMemStr( tFull.time, tFull.bytes ) )
			# @infiltrate
		end
	end
	# @infiltrate
	 
	attrLstRes = deepcopy( attrLst );
	valLstRes = deepcopy( valLst );
	push!( attrLstRes, "nRes" );
	push!( valLstRes, numRes );
	fResMain = "divBrefined";
	fOutModFull = fMod;
	if fOutMod != ""
		fOutModFull = fMod * fOutMod;
	end
	fResName = fNameFunc( fResMain, attrLstRes, valLstRes, fExt; fMod = fOutModFull );
	save( fResName, "divBLst", divBLst, "BfieldLst", BfieldLst );
	
	# return divBLst;
end

function divBRefinedStats( mSz, divLst, itNum, seed; fMod = "", dim = 3, numRes = 3, fModOut = "", fExt = jld2Type, thres=0.1 )
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim = dim, attrMoreLst = ["nRes"], valMoreLst = [numRes] );
	fMain = "divBrefined";
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	divBLst = load( fName, "divBLst" );
	BfieldLst = load( fName, "BfieldLst" );
	
	divBLstVcat = [ 
		vcat( ((x->vcat(x...)).( divBLst[iPol]) )... )
		for iPol = 1 : 2 ];
		
	BMaxLst = [ [ [ 
		maxDropDims( BfieldLst[iPol][it][n]; dims = [1,2] )
		for n = 1 : mSz ]
		for it = 1 : itNum ]
		for iPol = 1 : 2 ];
	BRemainLst = [ [ [ 
		zeros( 2*dim-1, size( divBLst[iPol][it][n], 1 ), numRes )
		for n = 1 : mSz ]
		for it = 1 : itNum ]
		for iPol = 1 : 2 ];
	BDiffLst = [[[
		BfieldLst[iPol][it][n][:,:,:,end] - BfieldLst[iPol][it][n][:,:,:,1]
		for n = 1 : mSz]
		for it = 1 : itNum]
		for iPol = 1 : 2];
	BDiffLstCat = [ 
		cat( (x->cat(x...; dims=3)).(BDiffLst[iPol]) ...; dims = 3 )
		for iPol = 1 : 2 ];
	BDiffLstVcat = [ reshape(BDiffLstCat[iPol],length(BDiffLstCat[iPol])) for iPol = 1 : 2 ]
	
	for iPol = 1:2, it = 1:itNum, n = 1:mSz
		for iRes = 1 : numRes, iLoc = 1 : size( divBLst[iPol][it][n], 1 )
			isFnd = false;
			valMax = BMaxLst;
			id = 1;
			for iBnd = 1:2, iDim = 1 : dim
				if isFnd || BfieldLst[iPol][it][n][iBnd,iDim,iLoc,iRes] != BMaxLst[iPol][it][n][iLoc,iRes]
					BRemainLst[iPol][it][n][id,iLoc,iRes] = BfieldLst[iPol][it][n][iBnd,iDim,iLoc,iRes];
					id += 1;
				else
					isFnd = true;
				end
			end
		end
	end
	
	BMaxLstVcat = [ 
		vcat( (x->vcat(x...)).( BMaxLst[iPol] )... )
		for iPol = 1:2 ];
	BRemainLstCat2 = [
		cat( (x->cat(x...; dims=2)).( BRemainLst[iPol] ) ...; dims=2 )
		for iPol = 1:2 ];
	BRemainLstVcat = [
		reshape( 
		cat( (x->cat(x...; dims=2)).( BRemainLst[iPol] ) ...; dims=2 ), :, numRes )
		for iPol = 1:2 ];
		
	BMaxLstErrVcat = [ 
		BMaxLstVcat[iPol][
			abs.(@view(divBLstVcat[iPol][:,3])).<(1-thres)
			,:] 
		for iPol = 1:2 ];
	@infiltrate
	BMaxLstNonErrVcat = [ 
		BMaxLstVcat[iPol][
			abs.(@view(divBLstVcat[iPol][:,3])).>(1-thres)
			,:] 
		for iPol = 1:2 ];
	BRemainLstErrCat2 = [ 
		BRemainLstCat2[iPol][:,
			abs.(@view(divBLstVcat[iPol][:,3])).<(1-thres)
			,:] 
		for iPol = 1:2 ];
	BRemainLstNonErrCat2 = [ 
		BRemainLstCat2[iPol][:,
			abs.(@view(divBLstVcat[iPol][:,3])).>(1-thres)
			,:] 
		for iPol = 1:2 ];
	BRemainLstErrVcat = [
		reshape( BRemainLstErrCat2[iPol], :, 3 )
		for iPol = 1:2];
	BRemainLstNonErrVcat = [
		reshape( BRemainLstNonErrCat2[iPol], :, 3 )
		for iPol = 1:2];
	BDiffLstErrCat = [ 
		BDiffLstCat[iPol][:,:,abs.(@view(divBLstVcat[iPol][:,3])).<(1-thres)]
		for iPol = 1 : 2 ];
	BDiffLstNonErrCat = [ 
		BDiffLstCat[iPol][:,:,abs.(@view(divBLstVcat[iPol][:,3])).>(1-thres)]
		for iPol = 1 : 2 ];
	BMaxDiffLstErrCat = [ 
		(-1)^(iPol-1) *minDropDims( (-1)^(iPol-1) * BDiffLstErrCat[iPol]; dims = [1,2] )
		for iPol = 1 : 2 ];
	
	BRemainDiffLstErrCat = [
		zeros(2*dim-1,length(BMaxDiffLstErrCat[iPol]))
		for iPol = 1 : 2];
	
	@infiltrate
	for iPol = 1 : 2, ii = 1 : length(BMaxDiffLstErrCat[iPol])
		isFnd = false;
		id = 1;
		for iBnd = 1:2, iDim = 1 : dim
			if isFnd || BDiffLstErrCat[iPol][iBnd,iDim,ii] != BMaxDiffLstErrCat[iPol][ii]
				BRemainDiffLstErrCat[iPol][id,ii] = BDiffLstErrCat[iPol][iBnd,iDim,ii];
				id += 1;
			else
				isFnd = true;
			end
		end
	end
	
	BRemainDiffLstErrVcat = [
		reshape( BRemainDiffLstErrCat[iPol], length(BRemainDiffLstErrCat[iPol]) )
		for iPol = 1 : 2];
	
	fOutMain = "divBRefinedVcat";
	
	fOutName = fNameFunc( fOutMain, attrLst, valLst, fExt; fMod = [fMod, fModOut] );
	
	save( fOutName, "divBLstVcat", divBLstVcat , "BMaxLst", BMaxLst, "BRemainLst", BRemainLst, "BMaxLstVcat", BMaxLstVcat, "BRemainLstVcat", BRemainLstVcat, "BMaxLstErrVcat", BMaxLstErrVcat, "BMaxLstNonErrVcat", BMaxLstNonErrVcat, "BRemainLstErrVcat", BRemainLstErrVcat, "BRemainLstNonErrVcat", BRemainLstNonErrVcat, "BDiffLstVcat", BDiffLstVcat, "BMaxDiffLstErrCat", BMaxDiffLstErrCat, "BRemainDiffLstErrCat", BRemainDiffLstErrCat, "BRemainDiffLstErrVcat", BRemainDiffLstErrVcat);
end
