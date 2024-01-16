module Loops_MC

using ShiftedArrays
using FilenameManip
using Random
using JLD2
using Statistics
using Distributions
using Utils
 
using Infiltrator

export loops_MC_methods, loops_MC, loops_MC_smart, loops_MC_staggeredCube

oFNameLoopsMain = "loopsSample";
oFNameLoopsStartMain = "loopsStartSample";
oFNameLoopsNumMain = "loopsNum";

fMainLoopsMC = "loops_MC";
attrLstLoops = ["divNum","itNum","cArea","cPerim","beta"];
attrLstLoopsFerro = deepcopy(attrLstLoops);
attrLstLoopsFerro = push!( attrLstLoopsFerro, "cFerro" );

struct ParamsLoops{N}
	nDim::Int64;
	nDimLayer::Int64;
	divNum::Int64;
	divLst::Vector{Int64};
	posLst::CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}};
	posLstShLst:: Matrix{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}}}; 
	linkDimLst::Vector{Vector{Int64}};
	linkDimShLst::Vector{Vector{Int64}};
	
	function ParamsLoops( divNum::Int64, nDim::Int64 )
		nDimLayer = nDim - 1;
		divLst = fill(divNum, nDim);
		posLst = CartesianIndices(ntuple(x->divNum,nDim));
		posLstShLst = Utils.arrShAdvRetLstFunc( posLst, nDim );
		linkDimLst = [ zeros(Int64, nDim-1) for dim = 1 : nDim ];
		for dim = 1:nDim
			iDLink = 1;
			for dimLink = 1 : nDim
				if dimLink != dim
					linkDimLst[dim][iDLink] = dimLink;
					iDLink += 1;
				end
			end
		end
		linkDimShLst = [ circshift( linkDimLst[dim], 1 ) for dim = 1 : nDim ];
		
		new{nDim}( nDim, nDimLayer, divNum, divLst, posLst, posLstShLst, linkDimLst, linkDimShLst );
	end
end

abstract type LoopsUpdater end

struct ABUpdater{N,N_1} <: LoopsUpdater
	posLstDimLst;
	posLayerLst::CartesianIndices{N_1,Tuple{Vararg{Base.OneTo{Int64},N_1}}};
	posABLst::Matrix{Vector{CartesianIndex{N}}};
	
	pFlipLst::Array{Float64};
	
	function ABUpdater( params::ParamsLoops; cArea, cPerim, cFerroSigned )
		posLstDimLst = [ selectdim( params.posLst, dim, it ) for dim = 1 : params.nDim, it = 1 : params.divNum ];
		posLayerLst = CartesianIndices( ntuple(x->params.divNum,params.nDim-1) );
		
		posABLst = [ Vector{CartesianIndex{params.nDim}}(undef,Int64(params.divNum^(params.nDim)/2)); for iAB = 1 : 2, iDim = 1 : params.nDim ];
		
		pFlipLst = genPFlipLst( cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
		
		for iDim = 1 : params.nDim
			iA = 1; iB = 1;
			for pos in params.posLst
				sumId = 0;
				for dim = 1 : params.nDim-1
					sumId += pos[params.linkDimLst[iDim][dim]];
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
		
		new{params.nDim,params.nDim-1}( posLstDimLst, posLayerLst, posABLst, pFlipLst );
	end
end

function updateLoops( BfieldLst, linkLst, params::ParamsLoops; updater::LoopsUpdater )
	error( "Loops_MC: updater not defined yet" );
end

function updateLoops( BfieldLst, linkLst, linkFerroLst, params::ParamsLoops; updater::ABUpdater )
	for dim = 1 : params.nDim
		for iAB = 1 : 2
			Threads.@threads for pos in updater.posABLst[iAB,dim]
				iArea = BfieldLst[dim][pos] + 1;
				iL = 1;
				for iLnkDim in 1 : params.nDimLayer
					dimLink = params.linkDimLst[dim][iLnkDim];
					dimLinkSh = params.linkDimShLst[dim][iLnkDim];
					iL += linkLst[dimLink][pos];
					iL += linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
				end
				iLFerro = 1;
				for iLnkDim = 1 : params.nDimLayer
					dimLink = params.linkDimLst[dim][iLnkDim];
					dimLinkSh = params.linkDimShLst[dim][iLnkDim];
					iLFerro += linkFerroLst[iLnkDim,dim][pos];
					iLFerro += linkFerroLst[iLnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]];
				end
				pSwitchRand = rand();
				pFlip = updater.pFlipLst[iLFerro, iL, iArea];
				if pSwitchRand < pFlip
					BfieldLst[dim][pos] = !BfieldLst[dim][pos];
					for iLnkDim in 1 : params.nDimLayer
						dimLink = params.linkDimLst[dim][iLnkDim];
						dimLinkSh = params.linkDimShLst[dim][iLnkDim];
						linkLst[dimLink][pos] = !linkLst[dimLink][pos];
						linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]] = !linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
					end
					for iLnkDim = 1 : params.nDimLayer
						linkFerroLst[iLnkDim,dim][pos] = !linkFerroLst[iLnkDim,dim][pos];
						dimLink = params.linkDimLst[dim][iLnkDim];
						dimLinkSh = params.linkDimShLst[dim][iLnkDim];
						linkFerroLst[iLnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]] = !linkFerroLst[iLnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]];
					end
				end
			end
		end
	end
end

function loops_MC_methods( divNum = 64, itNum = 10000; fMod = "", cArea = 1, cPerim = 1, cFerro = 0, beta = 1, itNumSample = 100, itStartSample = 50, isInit0 = false )
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
	
	itStep = Int64( floor(itNum / itNumSample) );
	lnSample = Int64( floor( itNum / itStep ) );
	
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

	abUpdater = ABUpdater( params; cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
	posLstDimLst = abUpdater.posLstDimLst;
	posLayerLst = abUpdater.posLayerLst;
	posABLst = abUpdater.posABLst;

	lnLayer = divNum^2;

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
		
		updateLoops( BfieldLst, linkLst, linkFerroLst, params; updater = abUpdater );
		
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
	valLst = Any[divNum, itNum, cArea, cPerim, beta];
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
	
	# oFNameLoopsMain = "loopsSample";
	oFNameLoops = fNameFunc( oFNameLoopsMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoops, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "linkSampleLst", linkSampleLst, "BfieldSampleLst", BfieldSampleLst );
	
	oFNameLoopsNum = fNameFunc( oFNameLoopsNumMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsNum, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst );
	
	# oFNameLoopsStartMain = "loopsStartSample";
	oFNameLoopsStart = fNameFunc( oFNameLoopsStartMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsStart, "linkStartSampleLst", linkStartSampleLst, "BfieldStartSampleLst", BfieldStartSampleLst );
	
	return fName;
end

function updateLinkFrom0ByB( BfieldLst, linkLst, linkFerroLst, params::ParamsLoops )
	for pos in params.posLst, dim = 1 : params.nDim
		if BfieldLst[dim][pos]
			for dimLink in params.linkDimLst[dim]
				linkLst[dimLink][pos] = !linkLst[dimLink][pos];
				for dimLinkSh in params.linkDimLst[dim]
					if dimLinkSh != dimLink
						linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]] = !linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
					end
				end
			end
			for lnkDim = 1 : params.nDim-1
				linkFerroLst[lnkDim,dim][pos] = !linkFerroLst[lnkDim,dim][pos];
				dimLink = params.linkDimLst[dim][lnkDim];
				for dimLinkSh in params.linkDimLst[dim]
					if dimLinkSh != dimLink
						linkFerroLst[lnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]] = linkFerroLst[lnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]];
					end
				end
			end
		end
	end
end

function genPFlipLst( ; cArea, cPerim, cFerroSigned = 0 )
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
	
	return pFlipLst;
end

function boolToOnePN( varBool::Bool )
	# return -(-1).^varBool;
	return varBool ? 1 : -1;
end

include("loops_MC_methods_seperate_funcs.jl")
export loops_MC, loops_MC_smart, loops_MC_staggeredCube

include("loops_MC_resave_funcs.jl")

include("loops_MC_old_funcs.jl")

end #endmodule
