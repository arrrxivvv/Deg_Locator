#Loops_MC module


function loops_MC_methods( divNum = 64, itNum = 10000; updaterType::(Type{T} where T <: LoopsUpdater), fMod = "", cArea = 1, cPerim = 1, cFerro = 0, isBeta = false, beta = nothing, itNumSample = 100, itStartSample = 50, isInit0 = false )
	cFerroSigned = cFerro;
	nDim = 3;
	params = ParamsLoops( divNum, nDim );
	nDimLayer = nDim-1;
	divLst = params.divLst;
	divTup = Tuple(divLst);
	posLst = params.posLst;
	posLstShLst = params.posLstShLst;
	linkDimLst = params.linkDimLst;
	linkDimShLst = params.linkDimShLst;
	
	itStep = max( Int64( floor(itNum / itNumSample) ), 1 );
	lnSample = Int64( floor( itNum / itStep ) );
	itStartSample = min( itStartSample, itNum );
	
	zakLstLst = zeros( Bool, divNum, divNum, nDim, itNum );

	BfieldLst = [ zeros( Bool, divTup ) for dim = 1 : nDim ];
	linkLst = [ zeros( Bool, divTup ) for dim = 1 : nDim ];
	linkFerroLst = [ zeros( Bool, divTup ) for lnkDim = 1 : nDim-1, BDim = 1 : nDim ];
	
	numBfieldLst = zeros(Int64, itNum, nDim);
	numLinkLst = zeros(Int64, itNum, nDim);
	
	probInit = exp( -cArea ) / ( 1 + exp( -cArea ) );
	distInit = Binomial( 1, probInit );
	if !isInit0
		for dim = 1 : nDim
			rand!( distInit, BfieldLst[dim] );
		end
	end
	updateLinkFrom0ByB( BfieldLst, linkLst, linkFerroLst, params );
	
	BfieldSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	linkSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	BfieldStartSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];
	linkStartSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];
	
	pFlipLst = genPFlipLst( cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );

	updater = updaterType( params; cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
	
	itSample = 1;
	for it = 1 : itNum
		print( "it = ", it, "         \r" )
		for dim = 1 : nDim
			@view(zakLstLst[:,:,dim,it]) .= dropdims( reduce( xor, BfieldLst[dim]; dims = dim, init = false ); dims = dim );
		end
		
		if mod( it, itStep ) == 0
			Threads.@threads for dim = 1 : nDim
				linkSampleLst[itSample][dim] .= linkLst[dim];
				BfieldSampleLst[itSample][dim] .= BfieldLst[dim];
			end
			itSample += 1;
		end
		
		if it <= itStartSample
			Threads.@threads for dim = 1 : nDim
				linkStartSampleLst[it][dim] .= linkLst[dim];
				BfieldStartSampleLst[it][dim] .= BfieldLst[dim];
			end
		end
		
		updateLoops( updater, BfieldLst, linkLst, linkFerroLst, params );
		
		for dim = 1 : nDim
			numBfieldLst[it,dim] = sum(BfieldLst[dim]);
			numLinkLst[it,dim] = sum( linkLst[dim] );
		end
	end
	
	xyDims = (1,2);
	zakMeanLst = dropdims( mean( zakLstLst; dims = xyDims ); dims = xyDims );
	
	zakLstSampleLst = zakLstLst[:,:,:,[itStep:itStep:itNum;]];
	
	fModOut = getFModLoopsMC( fMod, updaterType, isInit0 );
	
	fMain = fMainLoopsMC;
	valLst = Any[divNum, itNum, cArea, cPerim];
	attrLst, valLst = getAttrValLstLoopsMC( divNum, itNum, cArea, cPerim; beta = beta, cFerro = cFerro );
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
	
	save( fName, "divNum", divNum, "itNum", itNum, "cArea", cArea, "cPerim", cPerim, "beta", beta );
	
	fMainZakLst = fMainLoopsMC * "_" * "zakLstAllMean";
	fMainZakSample = fMainLoopsMC * "_" * "zakLstSampleMean";
	
	fNameZakLst = fNameFunc( fMainZakLst, attrLst, valLst, jld2Type; fMod = fModOut );
	fNameZakSample = fNameFunc( fMainZakSample, attrLst, valLst, jld2Type; fMod = fModOut );
	
	save( fNameZakLst, "zakLstLst", zakLstLst, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	save( fNameZakSample, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	
	oFNameLoops = fNameFunc( oFNameLoopsMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoops, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "linkSampleLst", linkSampleLst, "BfieldSampleLst", BfieldSampleLst );
	
	oFNameLoopsNum = fNameFunc( oFNameLoopsNumMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsNum, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst );
	
	oFNameLoopsStart = fNameFunc( oFNameLoopsStartMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsStart, "linkStartSampleLst", linkStartSampleLst, "BfieldStartSampleLst", BfieldStartSampleLst );
	
	return fName;
end
