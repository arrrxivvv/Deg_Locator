module Loops_MC

using ShiftedArrays
using FilenameManip
using Random
using JLD2
using Statistics
using Distributions
 
# using Infiltrator

export loops_MC, loops_MC_smart, loops_MC_staggeredCube

fMainLoopsMC = "loops_MC";
attrLstLoops = ["divNum","itNum","cArea","cPerim","beta"];
attrLstLoopsFerro = deepcopy(attrLstLoops);
attrLstLoopsFerro = push!( attrLstLoopsFerro, "cFerro" );

struct ParamsLoops
	nDim::Int64;
	divNum::Int64;
	divLst::Vector{Int64};
	posLst::AbstractArray{CartesianIndex{N}, N} where N;
	posLstShLst::Vector{AbstractArray{CartesianIndex{N}, N}} where N;
	
end

function loops_MC( divNum = 64, itNum = 10000; fMod = "", cArea = 1, cPerim = 1, beta = 1, itNumSample = 100, isInit0 = false )
	nDim = 3;
	divLst = fill(div, nDim);
	posLst = CartesianIndices((divNum,divNum,divNum));
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
	# @infiltrate

	posLstDimLst = [ selectdim( posLst, dim, it ) for dim = 1 : nDim, it = 1 : divNum ];

	lnLayer = divNum^2;

	zakLstLst = zeros( Bool, divNum, divNum, nDim, itNum );

	BfieldLst = [ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ];
	linkLst = [ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ];
	
	if !isInit0
		for dim = 1 : nDim
			rand!( distInit, BfieldLst[dim] );
		end
	end
	
	BfieldSampleLst = [[ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	linkSampleLst = [[ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	# @infiltrate
	for pos in posLst, dim = 1 : nDim
		if BfieldLst[dim][pos]
			for dimLink in linkDimLst[dim]
				linkLst[dimLink][pos] = !linkLst[dimLink][pos];
				for dimLinkSh in linkDimLst[dim]
					if dimLinkSh != dimLink
						linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]] = !linkLst[dimLink][posLstShLst[dimLinkSh,1][pos]];
					end
					# @infiltrate
				end
			end
			# @infiltrate
		end
	end
	# @infiltrate
	
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
			
			# @infiltrate
		end
		
		if mod( it, itStep ) == 0
			Threads.@threads for dim = 1 : nDim
				linkSampleLst[itSample][dim] .= linkLst[dim];
				BfieldSampleLst[itSample][dim] .= BfieldLst[dim];
			end
			itSample += 1;
		end
		
		for dim = 1 : nDim
			rand!( idLayerLst );
			idLayerLst .= mod.( idLayerLst, lnLayer ) .+ 1;
			dELst .= 0;
			rand!(rejLst);
			Threads.@threads for iLayer = 1 : divNum
				posLayer = posLstDimLst[dim,iLayer][idLayerLst[iLayer]];
				dELst[iLayer] -= cArea * boolToOnePN( BfieldLst[dim][posLayer] );
				for dimLink in linkDimLst[dim]
					dELst[iLayer] -= cPerim * boolToOnePN( linkLst[dimLink][posLayer] );
					for dimLinkSh in linkDimLst[dim]
						if dimLinkSh != dimLink
							dELst[iLayer] -= cPerim * boolToOnePN( linkLst[dimLink][posLstShLst[dimLinkSh,1][posLayer]] );
						end
					end
				end
				# @infiltrate
				if dELst[iLayer] < 0 || rejLst[iLayer] < exp( - beta * dELst[iLayer] )
					BfieldLst[dim][posLayer] = !BfieldLst[dim][posLayer];
					for dimLink in linkDimLst[dim]
						linkLst[dimLink][posLayer] = !linkLst[dimLink][posLayer];
						for dimLinkSh in linkDimLst[dim]
							if dimLinkSh != dimLink
								linkLst[dimLink][posLstShLst[dimLinkSh,1][posLayer]] = !linkLst[dimLink][posLstShLst[dimLinkSh,1][posLayer]];
							end
						end
					end
					# @infiltrate
				end
			end
			# @infiltrate it >= 200
		end
		for dim = 1 : nDim
			numBfieldLst[it,dim] = sum(BfieldLst[dim]);
			numLinkLst[it,dim] = sum( linkLst[dim] );
		end
	end
	
	xyDims = (1,2);
	zakMeanLst = dropdims( mean( zakLstLst; dims = xyDims ); dims = xyDims );
	
	zakLstSampleLst = zakLstLst[:,:,:,[itStep:itStep:itNum;]];
	
	fMain = fMainLoopsMC;
	attrLst = attrLstLoops;
	# ["divNum","itNum","cArea","cPerim","beta"];
	valLst = Any[divNum,itNum, cArea, cPerim,beta];
	fModOut = fMod;
	if isInit0
		if !isempty(fMod)
			fModOut *= "_";
		end
		fModOut *= "isInit0";
	end
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
	
	save( fName, "zakLstLst", zakLstLst, "divNum", divNum, "itNum", itNum, "cArea", cArea, "cPerim", cPerim, "beta", beta, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst, "linkSampleLst", linkSampleLst, "BfieldSampleLst", BfieldSampleLst );
	
	return fName;
end

function loops_MC_smart( divNum = 64, itNum = 10000; fMod = "", cArea = 1, cPerim = 1, cFerro = 0, beta = 1, itNumSample = 100, itStartSample = 50, isInit0 = false )
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

	posLstDimLst = [ selectdim( posLst, dim, it ) for dim = 1 : nDim, it = 1 : divNum ];
	posLayerLst = CartesianIndices( ntuple(x->divNum,nDim-1) );
	# posALst = Vector{CartesianIndex{nDim}}(undef,Int64(divNum^(nDim)/2));
	# posBLst = Vector{CartesianIndex{nDim}}(undef,Int64(divNum^(nDim)/2));
	
	# posABLst = [posALst, posBLst];
	posABLst = [ Vector{CartesianIndex{nDim}}(undef,Int64(divNum^(nDim)/2)); for iAB = 1 : 2, iDim = 1 : nDim ];

	# for iDim = 1 : nDim
		# iA = 1; iB = 1;
		# for pos in posLst
			# sumId = 0;
			# for dim = 1 : nDim-1
				# sumId += pos[dim];
			# end
			# if sumId % 2 == 0
				# posALst[iA,iDim] = pos;
				# iA += 1;
			# else
				# posBLst[iB,iDim] = pos;
				# iB += 1;
			# end
		# end
	# end
	
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
	
	# save( fName, "zakLstLst", zakLstLst, "divNum", divNum, "itNum", itNum, "cArea", cArea, "cPerim", cPerim, "beta", beta, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst, "linkSampleLst", linkSampleLst, "BfieldSampleLst", BfieldSampleLst );
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

function loops_MC_staggeredCube( divNum = 64, itNum = 10000; fMod = "", cArea = 1, cPerim = 1, cFerro = 0, beta = 1, itNumSample = 100, itStartSample = 50, isInit0 = false )
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

function loopsSamplesResaveFromFile( divNum, itNum; fMod = "", cArea = 1, cPerim = 1, beta = 1, isValLstFloat = false )
	fMain = fMainLoopsMC;
	attrLst = attrLstLoops;
	valLstAny = Any[ divNum, itNum, cArea, cPerim, beta ];
	valLstFloat = [ divNum, itNum, cArea, cPerim, beta ];
	
	valLst = valLstAny;
	if isValLstFloat
		valLst = valLstFloat;
	end
	
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	numBfieldLst = load( fName, "numBfieldLst" );
	numLinkLst = load( fName, "numLinkLst" );
	linkSampleLst = load( fName, "linkSampleLst" );
	BfieldSampleLst = load( fName, "BfieldSampleLst" );
	
	oFName = "loopsSample";
	
	oFName = fNameFunc( oFName, attrLst, valLst, jld2Type; fMod = fMod );
	
	save( oFName, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "linkSampleLst", linkSampleLst, "BfieldSampleLst", BfieldSampleLst );
end

function zakResaveFromFile( divNum, itNum; fMod = "", cArea = 1, cPerim = 1, beta = 1, isValLstFloat = false )
	fMain = fMainLoopsMC;
	attrLst = attrLstLoops;
	valLstAny = Any[ divNum, itNum, cArea, cPerim, beta ];
	valLstFloat = [ divNum, itNum, cArea, cPerim, beta ];
	
	valLst = valLstAny;
	if isValLstFloat
		valLst = valLstFloat;
	end
	
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	zakLstLst = load( fName, "zakLstLst" );
	zakMeanLst = load( fName, "zakMeanLst" );
	zakLstSampleLst = load( fName, "zakLstSampleLst" );
	
	fMainOut = fMainLoopsMC * "_" * "zakLstAllMean";
	fMainOutSample = fMainLoopsMC * "_" * "zakLstSampleMean";
	
	fNameOut = fNameFunc( fMainOut, attrLstLoops, valLst, jld2Type; fMod = fMod );
	fNameOutSample = fNameFunc( fMainOutSample, attrLstLoops, valLst, jld2Type; fMod = fMod );
	
	save( fNameOut, "zakLstLst", zakLstLst, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	save( fNameOutSample, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
end

function linkBfieldLstGenerateFromFile( divNum, itNum; cArea = 1, cPerim = 1, cFerro = 0, beta = 1, fMod = "", isValLstFloat = false, noSave = false, itLst = nothing )
	fMain = "loopsSample";
	attrLst = cFerro == 0 ? attrLstLoops : attrLstLoopsFerro;;
	valLstAny = Any[ divNum, itNum, cArea, cPerim, beta ];
	valLstFloat = [ divNum, itNum, cArea, cPerim, beta ];
	
	valLst = valLstAny;
	if isValLstFloat
		valLst = valLstFloat;
	end
	if cFerro != 0
		push!( valLst, cFerro );
	end
	
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	linkSampleLst = load( fName, "linkSampleLst" );
	BfieldSampleLst = load( fName, "BfieldSampleLst" );
	
	lnSample = length( linkSampleLst );
	nDim = 3;
	iDimLst = 1:nDim;
	
	linkNumLst = [ sum( linkSampleLst[it][dim] ) for it = 1 : lnSample, dim = 1 : nDim ];
	linkPlotLst = [ zeros(nDim,2,linkNumLst[it,dim]) for it = 1 : lnSample, dim = 1 : nDim ];
	
	# @infiltrate
	
	posLst = CartesianIndices( linkSampleLst[1][1] );
	
	divStep = 1 / divNum;
	
	coordDimLst = fill( divNum, nDim + 1 );
	coordDimLst[1] = nDim;
	coordFloatLst = zeros(coordDimLst...);
	
	dimLst = fill(divNum,nDim);
	divStepLst = 1 ./ dimLst;
	
	for iDim = iDimLst
		coordFloatLstDim = selectdim( coordFloatLst, 1, iDim );
		for iD = 1 : dimLst[iDim]
			selectdim( coordFloatLstDim, iDim, iD ) .= (iD-1) * divStepLst[iDim];
		end
	end
	
	coordFloatShLst = [ ShiftedArrays.circshift( coordFloatLst, ntuple( ( iD -> iD == iDim+1 ? -1 : 0 ), nDim ) )  for iDim = 1 : nDim ];
	
	# @time for it = 1 : lnSample
		# print( " it = ", it, "             \r" )
		# for dim in iDimLst
			# iLnk = 0;
			# for pos in posLst
				# if linkSampleLst[it][dim][pos]
					# iLnk += 1;
					# for iD = iDimLst
						# linkPlotLst[it,dim][iD,1,iLnk] = pos[iD];
						# linkPlotLst[it,dim][iD,2,iLnk] = pos[iD];
					# end
					# linkPlotLst[it,dim][dim,2,iLnk] += 1;
				# end
			# end
			# linkPlotLst[it,dim] .-= 1;
			# linkPlotLst[it,dim] .*= divStep;
		# end
	# end
	
	# dim = 1; iD = 1; iLnk = 0;
	# coordFloatShLstDim = coordFloatShLst[1];
	# linkPlotLstItDim = linkPlotLst[1,1];
	# linkSampleLstThis = linkSampleLst[1][1];
	# pos = posLst[1];
	if isnothing( itLst )
		itLst = 1 : lnSample
	end
	@time for it in itLst
		print( " it = ", it, "             \r" )
		for dim in iDimLst
			coordFloatShLstDim = coordFloatShLst[dim];
			linkPlotLstItDim = linkPlotLst[it,dim];
			linkSampleLstThis = linkSampleLst[it][dim];
			iLnk = 0;
			for pos in posLst
				if (linkSampleLstThis[pos])
					iLnk += 1;
					for iD in iDimLst
						linkPlotLstItDim[iD,1,iLnk] = coordFloatLst[iD,pos];
						# linkPlotLst[it,dim][iD,2,iLnk] = coordFloatLst[iD,pos];
						linkPlotLstItDim[iD,2,iLnk] = coordFloatShLstDim[iD,pos];
					end
				end
			end		
			# @infiltrate	
			# @view( linkPlotLst[it,dim][dim,2,:] ) .+= divStepLst[dim];
		end
	end
	
	# @infiltrate
	
	oFMain = "linkLstForMathematica";
	oFName = fNameFunc( oFMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	if !noSave
		save( oFName, "linkPlotLst", linkPlotLst );
	end
	
	GC.gc();
end

function linkBfieldPlotLstSliceResave( divNum, itNum; itLst, cArea = 1, cPerim = 1, cFerro = 0, beta = 1, isValLstFloat = false, fMod = "" )
	fMain = "linkLstForMathematica";
	attrLst = attrLstLoops;
	valLstAny = Any[ divNum, itNum, cArea, cPerim, beta ];
	valLstFloat = [ divNum, itNum, cArea, cPerim, beta ];
	
	valLst = valLstAny;
	if isValLstFloat
		valLst = valLstFloat;
	end
	
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	linkPlotLst = load( fName, "linkPlotLst" );
	
	# @infiltrate
	
	linkPlotLstSlice = linkPlotLst[itLst, :];
	
	oFMain = fMain * "_" * "sliced";
	oFName = fNameFunc( oFMain, attrLst, valLst, jld2Type; fMod = fMod ); 
	
	save( oFName, "linkPlotLstSlice", linkPlotLstSlice );
end

function zakAvgFromFile( divNum, itNum; fMod = "", cArea = 1, cPerim = 1, beta = 1, isValLstFloat = false )
	fMain = fMainLoopsMC;
	attrLst = attrLstLoops;
	valLstAny = Any[ divNum, itNum, cArea, cPerim, beta ];
	valLstFloat = [ divNum, itNum, cArea, cPerim, beta ];
	
	valLst = valLstAny;
	if isValLstFloat
		valLst = valLstFloat;
	end
	
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	zakLstLst = load(fName, "zakLstLst");
	
	xyDims = (1,2);
	zakMeanLst = dropdims( mean( zakLstLst; dims = xyDims ); dims = xyDims );
	
	fMainOut = fMainLoopsMC * "_zakMeanLst";
	fNameOut = fNameFunc( fMainOut, attrLst, valLstAny, jld2Type; fMod = fMod );
	
	save( fNameOut, "zakMeanLst", zakMeanLst );
end

function zakCorrFromFile( divNum, itNum; fMod = "", cArea = 1, cPerim = 1, beta = 1, isValLstFloat = false )
	fMain = fMainLoopsMC;
	attrLst = attrLstLoops;
	valLstAny = Any[ divNum, itNum, cArea, cPerim, beta ];
	valLstFloat = [ divNum, itNum, cArea, cPerim, beta ];
	
	valLst = valLstAny;
	if isValLstFloat
		valLst = valLstFloat;
	end
	
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	zakLstLstBool = load(fName, "zakLstLst");;
	
	shArr = [ (-x+1,-y+1) for x=1:divNum, y=1:divNum ];
	dArea = (1/divNum)^2;
	
	zakCorrLst = zeros( size(zakLstLstBool) );
	zakLstLst = zeros( Int64, size(zakLstLstBool) )
	zakLstLst .= 2 .* zakLstLstBool .- 1;
	
	zakLstLstShLst = [ ShiftedArrays.circshift( zakLstLst, shArr[x,y] ) for x = 1:divNum, y = 1:divNum ];
	
	zakXYColsLst = [ @view( zakLstLst[x:x,y:y,:,:] ) for x = 1 : divNum, y = 1 : divNum ];
	
	# @infiltrate
	
	zakLstCopy = similar(zakCorrLst);
	
	itDim = 4;
	
	for x = 1 : divNum, y = 1 : divNum
		print( "x = ", x, ", y = ", y, ",      \r" );
		# zakLstLstSh = ShiftedArrays.circshift( zakLstLst, shArr[x,y] );
		
		# @time zakLstCopy .= zakLstLstShLst[x,y];
		# @time zakLstCopy .*= zakXYColsLst[x,y];
		# @time zakLstCopy .*= dArea;
		# @time zakCorrLst .+= zakLstCopy;
		# Threads.@threads 
		for it = 1 : itNum
			@view(zakCorrLst[:,:,:,it]) .+= @view(zakLstLstShLst[x,y][:,:,:,it]) .* @view(zakLstLst[x:x,y:y,:,it]) .* dArea;
			# @infiltrate
		end
		# zakArrSh2 = ShiftedArrays.circshift( zakLstLst, shArr[x,y] );
			# Threads.@threads for it = 1 : itNum #Less
				# selectdim( zakCorrLst, itDim, it ) .= selectdim( zakCorrLst, itDim, it ) .+ @view(zakLstLst[x:x,y:y,:,it]) .* selectdim(zakArrSh2,itDim,it) .* dArea;
			# end
	end
	
	dimAvg = 4;
	zakCorrAvgLst = dropdims( mean( zakCorrLst; dims = dimAvg ); dims = dimAvg );
	zakCorrStdLst = dropdims( std( zakCorrLst; dims = dimAvg ); dims = dimAvg );
	
	# @infiltrate
	
	oFmain = "loops_MC_zakCorr";
	oFName = fNameFunc( oFmain, attrLst, valLst, jld2Type; fMod );
	save( oFName, "zakCorrLst", zakCorrLst, "zakCorrAvgLst", zakCorrAvgLst, "zakCorrStdLst", zakCorrStdLst );
end

function zakCorrResaveFromFile( divNum, itNum; fMod = "", cArea = 1, cPerim = 1, beta = 1, isValLstFloat = false, isAvgStdSaved = true )
	fMain = "loops_MC_zakCorr";
	
	attrLst = attrLstLoops;
	valLstAny = Any[ divNum, itNum, cArea, cPerim, beta ];
	valLstFloat = [ divNum, itNum, cArea, cPerim, beta ]; 
	
	valLst = valLstAny;
	if isValLstFloat
		valLst = valLstFloat;
	end
	
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	if isAvgStdSaved
		zakCorrAvgLst = load( fName, "zakCorrAvgLst" );
		zakCorrStdLst = load( fName, "zakCorrStdLst" );
	else
		zakCorrLst = load( fName, "zakCorrLst" );
		dimAvg = 4;
		zakCorrAvgLst = dropdims( mean( zakCorrLst; dims = dimAvg ); dims = dimAvg );
		zakCorrStdLst = dropdims( std( zakCorrLst; dims = dimAvg ); dims = dimAvg );
	end
	
	oFMain = "loops_MC_zakCorrStats";
	
	oFName = fNameFunc( oFMain, attrLst, valLst, jld2Type; fMod = fMod );
	save( oFName, "zakCorrAvgLst", zakCorrAvgLst, "zakCorrStdLst", zakCorrStdLst );
end

function zakCorrSampleFromFile( divNum, itNum; itSampleStep = 100, fMod = "", cArea = 1, cPerim = 1, beta = 1, isValLstFloat = false )
	fMain = "loops_MC_zakCorr";
	
	attrLst = attrLstLoops;
	valLstAny = Any[ divNum, itNum, cArea, cPerim, beta ];
	valLstFloat = [ divNum, itNum, cArea, cPerim, beta ]; 
	
	valLst = valLstAny;
	if isValLstFloat
		valLst = valLstFloat;
	end
	
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	zakCorrLst = load( fName, "zakCorrLst" );
	
	zakCorrSampleLst = zakCorrLst[:,:,:,itSampleStep:itSampleStep:itNum];
	
	oFMain = "loops_MC_zakCorrSample";
	oFName = fNameFunc( oFMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	save( oFName, "zakCorrSampleLst", zakCorrSampleLst );
end

function zakSampleLstFromFile( divNum, itNum; itNumSample = 100, fMod = "", cArea = 1, cPerim = 1, beta = 1, isValLstFloat = false )
	fMain = fMainLoopsMC;
	attrLst = attrLstLoops;
	valLstAny = Any[ divNum, itNum, cArea, cPerim, beta ];
	valLstFloat = [ divNum, itNum, cArea, cPerim, beta ];
	
	valLst = valLstAny;
	if isValLstFloat
		valLst = valLstFloat;
	end
	
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	zakLstLst = load(fName, "zakLstLst");
	
	itStep = Int64( floor(itNum / itNumSample) );
	
	zakLstSampleLst = zakLstLst[:,:,:,[itStep:itStep:itNum;]];
	
	fMainOut = fMain * "_zakSamples";
	fNameOut = fNameFunc( fMainOut, attrLst, valLstAny, jld2Type; fMod = fMod );
	
	save( fNameOut, "zakLstSampleLst", zakLstSampleLst );
end

function MCLinkUpdate( BfieldLst, linkLst,  )
	for dim = 1 : nDim
	end
end

function boolToOnePN( varBool::Bool )
	return -(-1).^varBool;
end

end #endmodule
