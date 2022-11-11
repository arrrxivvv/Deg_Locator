# using Infiltrator

function degDivRefineFromFile( mSz, divLst, itNum, seed; fMod = "", dim = 3, fExt = jld2Type )
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
	
	numRes = 3;
	
	divBlst = [ [ 
		[ zeros( size(locLstPol[iPol][it][lev],1), numRes ) for lev = 1:mSz ] 
			for it = 1:itNum ] 
				for iPol = 1:2 ];
	
	divLstCube = fill(2,params.nDim);
	# paramsCubeThr = thrStructCopy( degParamsNonInit( params.N, divLst, params.nDim; isNonPeriodic = true ) );
	# matsGridCubeThr = thrStructCopy( 	 
		# matsGridHThreaded( 
			# degParamsNonInit( params.N, divLst, params.nDim; isNonPeriodic = true ), 
			# HmatThr ) 
	# );
	
	matsGridCubeLst = Vector{DegMatsOnGrid}(undef, numRes);
	
	divLstCube = ones(params.nDim);
	for iRes = 1 : numRes
		divLstCube .= 2^iRes;
		matsGridCubeLst[iRes] = 
			matsGridHThreaded( 
				degParamsNonInit( params.N, divLstCube, params.nDim; isNonPeriodic = true ), 
				HmatThr
			);
	end
	
	numSh = 2^params.nDim;
	maxLstTmp = zeros(params.nDim);
	iShTmp = zeros(Int64,params.nDim);
	for it = 1 : itNum
		HLst = HLstLst[it];
		HmatFun = ( (H,xLst) -> Hmat_3comb_ratio!(H,xLst,HLst) );
		# HmatGridFun reevaluate
		HmatGridFun = (xLst) -> HmatFun( getHLoc( iLst, degMatsGrid ), 
			degMatsGrid.mesh[
				linIdFromIdVec( iLst, params )] 
		);
		startNextEigen( degMatsGrid );
		for n = 1 : mSz, iPol = 1 : 2
			locLst = locLstPol[iPol][it][n];
			# Threads.@threads 
			for iLoc = 1 : size(locLst,1)
				loc = @view(locLst[iLoc,:]);
				locSh = ones( Int64, params.nDim ) ;
				dLoc = ones( Int64, params.nDim );
				
				for iSh = 0 : numSh-1
					locSh .= loc;
					iShTmp .= iSh;
					for iDimSh = 1 : params.nDim
						if iShTmp & 1 != 0
							locSh[iDimSh] += 1;
						end
						iShTmp = iShTmp >> 1;
					end
					wrapIdVec!(locSh, params);
					
					HmatGridFun( locSh );
					eigenAtLoc( locSh, degMatsGrid );
				end
			end
			for iLoc = 1 : size(locLst, 1)
				loc = @view(locLst[iLoc,:]);
				linId = linIdFromIdVec( loc, params );
				maxLstTmp .= params.mesh[linId] .+ params.stepLst;
				for iRes = 1 : numRes
					matsCube = matsGridCubeLst[iRes];
					updateParamsRange( params.mesh[linId], maxLstTmp, matsCube.params );
					
					startNextEigen( matsCube );
					matsGridTransfer( matsCube, degMatsGrid, loc );
					
					
				end
			end
			# @infiltrate
		end
	end
	# @infiltrate
end
