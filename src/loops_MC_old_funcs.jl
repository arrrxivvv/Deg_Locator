#Loops_MC funcs
# obsolete

function loops_MC_smart_old( divNum = 64, itNum = 10000; fMod = "", cArea = 1, cPerim = 1, cFerro = 0, beta = 1, itNumSample = 100, itStartSample = 50, isInit0 = false )
	cFerroSigned = cFerro;
	nDim = 3;
	divLst = fill(div, nDim);
	posLst = CartesianIndices(ntuple(x->divNum,nDim));
	posLstShLst = [ ShiftedArrays.circshift( posLst, ntuple( ( i -> ( i == dim ? (-1)^iPol : 0 ) ), nDim ) ) for dim = 1 : nDim, iPol = 1 : 2 ];
	
	probInit = exp( -cArea ) / ( 1 + exp( -cArea ) );
	
	distInit = Binomial( 1, probInit );
	
	itStep = Int64( floor(itNum / itNumSample) );
	lnSample = Int64( floor( itNum / itStep ) );
	
	linkDimLst = [ zeros(Int64, nDim-1) for dim = 1 : nDim ];
	dim = 1;
	for dim = 1:nDim
		iDLink = 1;
		for dimLink = 1 : nDim
			if dimLink != dim
				linkDimLst[dim][iDLink] = dimLink;
				iDLink += 1;
			end
		end
	end
	
	numFlipTypes = 2;
	numPerimTypes = 2*numFlipTypes+1;
	flipStep = 2;
	pFlipLst = zeros( numPerimTypes, numPerimTypes, numFlipTypes );
	Bval = 1;
	for iB = 1 : numFlipTypes
		Eval = 0;
		Eval -= cArea * Bval;
		lnkVal = numFlipTypes * flipStep;
		for iL = 1 : numPerimTypes
			Eval -= cPerim * lnkVal;
			lnkFerroVal = numFlipTypes * flipStep;
			for iLFerro = 1 : numPerimTypes
				Eval -= cFerro * lnkFerroVal;
				dE = -2*Eval;
				pFlipLst[iLFerro, iL, iB] = 1 / ( 1 + exp( dE ) );
				lnkFerroVal -= flipStep;
			end
			lnkVal -= flipStep;
		end
		Bval *= -1;
	end

	posLstDimLst = [ selectdim( posLst, dim, it ) for dim = 1 : nDim, it = 1 : divNum ];
	posLayerLst = CartesianIndices( ntuple(x->divNum,nDim-1) );
	
	posABLst = [ Vector{CartesianIndex{nDim}}(undef,Int64(divNum^(nDim)/2)) for iAB = 1 : 2, iDim = 1 : nDim ];
	
	for iDim = 1 : nDim
		iA = 1; iB = 1;
		for pos in posLst
			sumId = 0;
			for dim = 1 : nDim-1
				sumId += pos[linkDimLst[iDim][dim]];
			end
			if sumId % 2 == 0
				posABLst[1,iDim][iA] = pos;
				iA += 1;
			else
				posABLst[2,iDim][iB] = pos;
				iB += 1;
			end
		end
	end

	lnLayer = divNum^2;

	zakLstLst = zeros( Bool, divNum, divNum, nDim, itNum );

	BfieldLst = [ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ];
	linkLst = [ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ];
	linkFerroLst = [ zeros( Bool, divNum, divNum, divNum ) for lnkDim = 1 : nDim-1, BDim = 1 : nDim ];
	
	if !isInit0
		for dim = 1 : nDim
			rand!( distInit, BfieldLst[dim] );
		end
	end
	
	BfieldSampleLst = [[ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	linkSampleLst = [[ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	BfieldStartSampleLst = [[ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];
	linkStartSampleLst = [[ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];
	# @infiltrate
	for pos in posLst, dim = 1 : nDim
		if BfieldLst[dim][pos]
			for dimLink in linkDimLst[dim]
				linkLst[dimLink][pos] = !linkLst[dimLink][pos];
				for dimLinkSh in linkDimLst[dim]
					if dimLinkSh != dimLink
						linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]] = !linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]];
					end
				end
			end
			for lnkDim = 1 : nDim-1
				linkFerroLst[lnkDim,dim][pos] = !linkFerroLst[lnkDim,dim][pos];
				dimLink = linkDimLst[dim][lnkDim];
				for dimLinkSh in linkDimLst[dim]
					if dimLinkSh != dimLink
						linkFerroLst[lnkDim,dim][posLstShLst[dimLinkSh,1][pos]] = linkFerroLst[lnkDim,dim][posLstShLst[dimLinkSh,1][pos]];
					end
				end
			end
		end
	end
	
	numBfieldLst = zeros(Int64, itNum, nDim);
	numLinkLst = zeros(Int64, itNum, nDim);

	idLayerLst = zeros( UInt, divNum );
	dELst = zeros(divNum);
	
	rejLst = zeros(divNum);
	itSample = 1;
	@time for it = 1 : itNum
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
		
		for dim = 1 : nDim
			rand!( idLayerLst );
			idLayerLst .= mod.( idLayerLst, lnLayer ) .+ 1;
			dELst .= 0;
			rand!(rejLst);
			# Threads.@threads for iLayer = 1 : divNum
			for iAB = 1 : 2
				Threads.@threads for pos in posABLst[iAB,dim]
					dE = 0;
					dE -= cArea * boolToOnePN( BfieldLst[dim][pos] );
					for dimLink in linkDimLst[dim]
						dE -= cPerim * boolToOnePN( linkLst[dimLink][pos] );
						for dimLinkSh in linkDimLst[dim]
							if dimLinkSh != dimLink
								dE -= cPerim * boolToOnePN( linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]] );
							end
						end
					end
					for lnkDim = 1 : nDim-1
						dimLink = linkDimLst[dim][lnkDim];
						dE -= cFerroSigned * boolToOnePN( linkFerroLst[lnkDim,dim][pos] );
						for dimLinkSh in linkDimLst[dim]
							if dimLinkSh != dimLink
								dE -= cFerroSigned * boolToOnePN( linkFerroLst[lnkDim,dim][posLstShLst[dimLinkSh,1][pos]] );
							end
						end
					end
					pSwitch = rand();
					if pSwitch < exp( - beta * dE ) / ( 1 + exp( - beta * dE ) )
						BfieldLst[dim][pos] = !BfieldLst[dim][pos];
						for dimLink in linkDimLst[dim]
							linkLst[dimLink][pos] = !linkLst[dimLink][pos];
							for dimLinkSh in linkDimLst[dim]
								if dimLinkSh != dimLink
									linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]] = !linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]];
								end
							end
						end
						for lnkDim = 1 : nDim-1
							linkFerroLst[lnkDim,dim][pos] = !linkFerroLst[lnkDim,dim][pos];
							dimLink = linkDimLst[dim][lnkDim];
							for dimLinkSh in linkDimLst[dim]
								if dimLinkSh != dimLink
									linkFerroLst[lnkDim,dim][posLstShLst[dimLinkSh,1][pos]] = !linkFerroLst[lnkDim,dim][posLstShLst[dimLinkSh,1][pos]];
								end
							end
						end
					end
				end
			end
			# end
		end
		for dim = 1 : nDim
			numBfieldLst[it,dim] = sum(BfieldLst[dim]);
			numLinkLst[it,dim] = sum( linkLst[dim] );
		end
	end
	# @infiltrate
	
	xyDims = (1,2);
	zakMeanLst = dropdims( mean( zakLstLst; dims = xyDims ); dims = xyDims );
	
	zakLstSampleLst = zakLstLst[:,:,:,[itStep:itStep:itNum;]];
	
	fModOut = fMod;
	if isInit0
		if !isempty(fMod)
			fModOut *= "_";
		end
		fModOut *= "isInit0";
	end
	
	fMain = fMainLoopsMC;
	attrLst = cFerro == 0 ? attrLstLoops : attrLstLoopsFerro;
	# ["divNum","itNum","cArea","cPerim","beta"];
	valLst = Any[divNum, itNum, cArea, cPerim, beta];
	if cFerro != 0
		push!( valLst, cFerro );
	end
	# @infiltrate
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
	
	# save( fName, "zakLstLst", zakLstLst, "divNum", divNum, "itNum", itNum, "cArea", cArea, "cPerim", cPerim, "beta", beta, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst, "linkSampleLst", linkSampleLst, "BfieldSampleLst", BfieldSampleLst );
	save( fName, "divNum", divNum, "itNum", itNum, "cArea", cArea, "cPerim", cPerim, "beta", beta );
	
	fMainZakLst = fMainLoopsMC * "_" * "zakLstAllMean";
	fMainZakSample = fMainLoopsMC * "_" * "zakLstSampleMean";
	
	fNameZakLst = fNameFunc( fMainZakLst, attrLst, valLst, jld2Type; fMod = fModOut );
	fNameZakSample = fNameFunc( fMainZakSample, attrLst, valLst, jld2Type; fMod = fModOut );
	
	save( fNameZakLst, "zakLstLst", zakLstLst, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	save( fNameZakSample, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	
	# oFNameLoopsMain = "loopsSample";
	oFNameLoops = fNameFunc( oFNameLoopsMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoops, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "linkSampleLst", linkSampleLst, "BfieldSampleLst", BfieldSampleLst );
	
	# oFNameLoopsStartMain = "loopsStartSample";
	oFNameLoopsStart = fNameFunc( oFNameLoopsStartMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsStart, "linkStartSampleLst", linkStartSampleLst, "BfieldStartSampleLst", BfieldStartSampleLst );
	
	return fName;
end

function loops_MC_staggeredCube_old( divNum = 64, itNum = 10000; fMod = "", cArea = 1, cPerim = 1, cFerro = 0, beta = 1, itNumSample = 100, itStartSample = 50, isInit0 = false )
	if isempty( fMod )
		fMod *= "_";
	end
	fMod *= "staggeredCube";
	
	cFerroSigned = cFerro;
	nDim = 3;
	iDimLst = 1:nDim;
	iIsShLst = 1:2;
	divLst = fill(divNum, nDim);
	posLst = CartesianIndices(ntuple(x->divNum,nDim));
	posLstAdvanced = [ ShiftedArrays.circshift( posLst, ntuple( ( i -> ( i == dim ? -1 : 0 ) ), nDim ) ) for dim = 1 : nDim ];
	posLstAdvOrNot = pushfirst!( posLstAdvanced, ShiftedArrays.circshift( posLst, ntuple(x->0,nDim) ) );
	posShOrNotLst = [ [ posLst, posLstAdvanced[dim] ] for dim = 1 : nDim ];
	posLstShLst = [ ShiftedArrays.circshift( posLst, ntuple( ( i -> ( i == dim ? (-1)^iPol : 0 ) ), nDim ) ) for dim = 1 : nDim, iPol = 1 : 2 ];
	
	probInit = exp( -cArea ) / ( 1 + exp( -cArea ) );
	
	distInit = Binomial( 1, probInit );
	
	itStep = Int64( floor(itNum / itNumSample) );
	lnSample = Int64( floor( itNum / itStep ) );
	
	linkDimLst = [ zeros(Int64, nDim-1) for dim = 1 : nDim ];
	dim = 1;
	for dim = 1:nDim
		iDLink = 1;
		for dimLink = 1 : nDim
			if dimLink != dim
				linkDimLst[dim][iDLink] = dimLink;
				iDLink += 1;
			end
		end
	end
	
	numFlipTypes = 2;
	numPerimTypes = 2*numFlipTypes+1;
	flipStep = 2;
	pFlipLst = zeros( numPerimTypes, numPerimTypes, numFlipTypes );
	Bval = 1;
	for iB = 1 : numFlipTypes
		lnkVal = numFlipTypes * flipStep;
		for iL = 1 : numPerimTypes
			lnkFerroVal = numFlipTypes * flipStep;
			for iLFerro = 1 : numPerimTypes
				Eval = - ( cArea * Bval + cPerim * lnkVal + cFerroSigned * lnkFerroVal );
				dE = -2*Eval;
				pFlipLst[iLFerro, iL, iB] = 1 / ( 1 + exp( dE ) );
				lnkFerroVal -= flipStep;
			end
			lnkVal -= flipStep;
		end
		Bval *= -1;
	end

	posLstDimLst = [ selectdim( posLst, dim, it ) for dim = 1 : nDim, it = 1 : divNum ];
	
	coordsStagA = ntuple( iDim -> 1:2:divLst[iDim], nDim );
	coordsStagB = ntuple( iDim -> 2:2:divLst[iDim], nDim );
	posStagCubeLstA = @view( posLst[coordsStagA...] );
	posStagCubeLstB = @view( posLst[coordsStagB...] );
	posStagCubeLst = [ cat( @view( posLstAdvOrNot[iAdv][coordsStagA...]), @view(posLstAdvOrNot[iAdv][coordsStagB...]); dims = nDim+1 ) for iAdv = 1 : nDim+1 ];
	# cat( posStagCubeLstA, posStagCubeLstB, nDim+1 );

	zakLstLst = zeros( Bool, divNum, divNum, nDim, itNum );

	BfieldLst = [ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ];
	linkLst = [ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ];
	linkFerroLst = [ zeros( Bool, divNum, divNum, divNum ) for lnkDim = 1 : nDim-1, BDim = 1 : nDim ];
	
	if !isInit0
		for dim = 1 : nDim
			rand!( distInit, BfieldLst[dim] );
		end
	end
	
	BfieldSampleLst = [[ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	linkSampleLst = [[ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	BfieldStartSampleLst = [[ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];
	linkStartSampleLst = [[ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];
	# @infiltrate
	for pos in posLst, dim = 1 : nDim
		if BfieldLst[dim][pos]
			for dimLink in linkDimLst[dim]
				linkLst[dimLink][pos] = !linkLst[dimLink][pos];
				for dimLinkSh in linkDimLst[dim]
					if dimLinkSh != dimLink
						linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]] = !linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]];
					end
				end
			end
			for lnkDim = 1 : nDim-1
				linkFerroLst[lnkDim,dim][pos] = !linkFerroLst[lnkDim,dim][pos];
				dimLink = linkDimLst[dim][lnkDim];
				for dimLinkSh in linkDimLst[dim]
					if dimLinkSh != dimLink
						linkFerroLst[lnkDim,dim][posLstShLst[dimLinkSh,1][pos]] = linkFerroLst[lnkDim,dim][posLstShLst[dimLinkSh,1][pos]];
					end
				end
			end
		end
	end
	
	numBfieldLst = zeros(Int64, itNum, nDim);
	numLinkLst = zeros(Int64, itNum, nDim);

	idLayerLst = zeros( UInt, divNum );
	dELst = zeros(divNum);
	
	rejLst = zeros(divNum);
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
		
		for iAdv = 1 : nDim+1
			Threads.@threads for posCube in posStagCubeLst[iAdv]
				dim = rand(iDimLst);
				isShId = rand(iIsShLst);
				pos = posShOrNotLst[dim][isShId][posCube];
				# @infiltrate
				dE = 0;
				dE -= cArea * boolToOnePN( BfieldLst[dim][pos] );
				for dimLink in linkDimLst[dim]
					dE -= cPerim * boolToOnePN( linkLst[dimLink][pos] );
					for dimLinkSh in linkDimLst[dim]
						if dimLinkSh != dimLink
							dE -= cPerim * boolToOnePN( linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]] );
						end
					end
				end
				for lnkDim = 1 : nDim-1
					dimLink = linkDimLst[dim][lnkDim];
					dE -= cFerroSigned * boolToOnePN( linkFerroLst[lnkDim,dim][pos] );
					for dimLinkSh in linkDimLst[dim]
						if dimLinkSh != dimLink
							dE -= cFerroSigned * boolToOnePN( linkFerroLst[lnkDim,dim][posLstShLst[dimLinkSh,1][pos]] );
						end
					end
				end
				pSwitch = rand();
				if pSwitch < exp( - beta * dE ) / ( 1 + exp( - beta * dE ) )
					BfieldLst[dim][pos] = !BfieldLst[dim][pos];
					for dimLink in linkDimLst[dim]
						linkLst[dimLink][pos] = !linkLst[dimLink][pos];
						for dimLinkSh in linkDimLst[dim]
							if dimLinkSh != dimLink
								linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]] = !linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]];
							end
						end
					end
					for lnkDim = 1 : nDim-1
						linkFerroLst[lnkDim,dim][pos] = !linkFerroLst[lnkDim,dim][pos];
						dimLink = linkDimLst[dim][lnkDim];
						for dimLinkSh in linkDimLst[dim]
							if dimLinkSh != dimLink
								linkFerroLst[lnkDim,dim][posLstShLst[dimLinkSh,1][pos]] = !linkFerroLst[lnkDim,dim][posLstShLst[dimLinkSh,1][pos]];
							end
						end
					end
				end
			end
		end
		
		for dim = 1 : nDim
			numBfieldLst[it,dim] = sum(BfieldLst[dim]);
			numLinkLst[it,dim] = sum( linkLst[dim] );
		end
	end
	
	xyDims = (1,2);
	zakMeanLst = dropdims( mean( zakLstLst; dims = xyDims ); dims = xyDims );
	
	zakLstSampleLst = zakLstLst[:,:,:,[itStep:itStep:itNum;]];
	
	fModOut = fMod;
	if isInit0
		if !isempty(fMod)
			fModOut *= "_";
		end
		fModOut *= "isInit0";
	end
	
	fMain = fMainLoopsMC;
	attrLst = cFerro == 0 ? attrLstLoops : attrLstLoopsFerro;
	# ["divNum","itNum","cArea","cPerim","beta"];
	valLst = Any[divNum,itNum, cArea, cPerim,beta];
	if cFerro != 0
		push!( valLst, cFerro );
	end
	# @infiltrate
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
	
	save( fName, "divNum", divNum, "itNum", itNum, "cArea", cArea, "cPerim", cPerim, "beta", beta );
	
	fMainZakLst = fMainLoopsMC * "_" * "zakLstAllMean";
	fMainZakSample = fMainLoopsMC * "_" * "zakLstSampleMean";
	
	fNameZakLst = fNameFunc( fMainZakLst, attrLst, valLst, jld2Type; fMod = fModOut );
	fNameZakSample = fNameFunc( fMainZakSample, attrLst, valLst, jld2Type; fMod = fModOut );
	
	save( fNameZakLst, "zakLstLst", zakLstLst, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	save( fNameZakSample, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	
	oFNameLoops = "loopsSample";
	
	oFNameLoops = fNameFunc( oFNameLoops, attrLst, valLst, jld2Type; fMod = fModOut );
	
	save( oFNameLoops, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "linkSampleLst", linkSampleLst, "BfieldSampleLst", BfieldSampleLst );
	
	oFNameLoopsStart = "loopsStartSample";
	oFNameLoopsStart = fNameFunc( oFNameLoopsStart, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsStart, "linkStartSampleLst", linkStartSampleLst, "BfieldStartSampleLst", BfieldStartSampleLst  )
	
	return fName;
end

