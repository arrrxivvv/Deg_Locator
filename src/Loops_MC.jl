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

dirLog = "./log/";

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

function updateLoops( updater::LoopsUpdater, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	error( "Loops_MC: updater not defined yet" );
end

function getUpdaterFMod( updaterType::(Type{T} where T <: LoopsUpdater) )
	error( "Loops_MC: updater not defined yet" );
end

struct ABUpdater{N,N_1} <: LoopsUpdater
	posLstDimLst;
	posLayerLst::CartesianIndices{N_1,Tuple{Vararg{Base.OneTo{Int64},N_1}}};
	posABLst::Matrix{Vector{CartesianIndex{N}}};
	
	pFlipLst::Array{Float64};
	
	function ABUpdater( params::ParamsLoops; cArea::Float64, cPerim::Float64, cFerroSigned::Float64 )
		posLstDimLst = [ selectdim( params.posLst, dim, it ) for dim = 1 : params.nDim, it = 1 : params.divNum ];
		posLayerLst = CartesianIndices( ntuple(x->params.divNum,params.nDim-1) );
		
		posABLst = [ Vector{CartesianIndex{params.nDim}}(undef,Int64(params.divNum^(params.nDim)/2)) for iAB = 1 : 2, iDim = 1 : params.nDim ];
		
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

function updateLoops( updater::ABUpdater, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	for dim = 1 : params.nDim
		for iAB = 1 : 2
			Threads.@threads for pos in updater.posABLst[iAB,dim]
				# @time begin
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
				# end
				# @time begin
				pSwitchRand = rand();
				if pSwitchRand < updater.pFlipLst[iLFerro, iL, iArea];
					flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = dim );
				end
				# end
				# @infiltrate
			end
		end
	end
end

function getUpdaterFMod( updaterType::Type{ABUpdater} )
	return "upAB";
end

struct StaggeredCubeUpdater{N,Nplus1} <: LoopsUpdater
	posLstSh0::CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}};
	posLstAdvOrNot::Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}} }};
	posShOrNotLst::Vector{Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}}}};
	posStagCubeLst::Vector{Array{CartesianIndex{N},Nplus1}};
	idStagLst::CartesianIndices{Nplus1,NTuple{Nplus1,Base.OneTo{Int64}}};
	
	pFlipLst::Array{Float64};
	
	iDimLst::UnitRange{Int64};
	iIsShLst::UnitRange{Int64};
	randIDimLst::Array{Int64,Nplus1};
	randIShLst::Array{Int64,Nplus1};
	
	function StaggeredCubeUpdater( params::ParamsLoops; cArea::Float64, cPerim::Float64, cFerroSigned::Float64 )
		posLstSh0 = ShiftedArrays.circshift( params.posLst, ntuple(x->0,params.nDim) );
		posLstAdvanced = [ ShiftedArrays.circshift( params.posLst, ntuple( ( i -> ( i == dim ? -1 : 0 ) ), params.nDim ) ) for dim = 1 : params.nDim ];
		posLstAdvOrNot = pushfirst!( posLstAdvanced, ShiftedArrays.circshift( params.posLst, ntuple(x->0,params.nDim) ) );
		posShOrNotLst = [ [ posLstSh0, posLstAdvanced[dim] ] for dim = 1 : params.nDim ];
		posShOrNotArr = [ (iSh == 1 ? posLstSh0 : posLstAdvanced[dim]) for dim = 1 : params.nDim, iSh = 1:2 ];
		
		coordsStagA = ntuple( iDim -> 1:2:params.divLst[iDim], params.nDim );
		coordsStagB = ntuple( iDim -> 2:2:params.divLst[iDim], params.nDim );
		posStagCubeLstA = @view( params.posLst[coordsStagA...] );
		posStagCubeLstB = @view( params.posLst[coordsStagB...] );
		posStagCubeLst = [ cat( @view( posLstAdvOrNot[iAdv][coordsStagA...]), @view(posLstAdvOrNot[iAdv][coordsStagB...]); dims = params.nDim+1 ) for iAdv = 1 : params.nDim+1 ];
		idStagLst = CartesianIndices( posStagCubeLst[1] );
		
		pFlipLst = genPFlipLst( cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
		
		iDimLst = 1:params.nDim;
		iIsShLst = 1:2;
		randIDimLst = similar( posStagCubeLst[1], Int64 );
		randIShLst = similar(randIDimLst);
		
		new{params.nDim,params.nDim+1}( posLstSh0, posLstAdvOrNot, posShOrNotLst, posStagCubeLst, idStagLst, pFlipLst, iDimLst, iIsShLst, randIDimLst, randIShLst );
	end
end

function updateLoops( updater::StaggeredCubeUpdater, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	for iAdv = 1 : params.nDim+1
		rand!( updater.randIDimLst, updater.iDimLst );
		rand!( updater.randIShLst, updater.iIsShLst );
		Threads.@threads for idStag in updater.idStagLst
			posCube = updater.posStagCubeLst[iAdv][idStag];
			pos = updater.posShOrNotLst[updater.randIDimLst[idStag]][updater.randIShLst[idStag]][posCube];
			iArea = BfieldLst[updater.randIDimLst[idStag]][pos] + 1;
			iL = 1;
			for iLnkDim in 1 : params.nDimLayer
				dimLink = params.linkDimLst[updater.randIDimLst[idStag]][iLnkDim];
				dimLinkSh = params.linkDimShLst[updater.randIDimLst[idStag]][iLnkDim];
				iL += linkLst[dimLink][pos];
				iL += linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
			end
			iLFerro = 1;
			for iLnkDim = 1 : params.nDimLayer
				dimLink = params.linkDimLst[updater.randIDimLst[idStag]][iLnkDim];
				dimLinkSh = params.linkDimShLst[updater.randIDimLst[idStag]][iLnkDim];
				iLFerro += linkFerroLst[iLnkDim,updater.randIDimLst[idStag]][pos];
				iLFerro += linkFerroLst[iLnkDim,updater.randIDimLst[idStag]][params.posLstShLst[dimLinkSh,1][pos]];
			end
			pSwitch = rand();
			pFlip = updater.pFlipLst[iLFerro, iL, iArea];
			if pSwitch < pFlip
				flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = updater.randIDimLst[idStag] );
			end
			# end
			# @infiltrate
		end
		# @infiltrate
	end
end

function getUpdaterFMod( updaterType::Type{StaggeredCubeUpdater} )
	return "upStagCube";
end

function loops_MC_methods( divNum = 64, itNum = 10000; updaterType::(Type{T} where T <: LoopsUpdater), fMod = "", cArea = 1, cPerim = 1, cFerro = 0, beta = 1, itNumSample = 100, itStartSample = 50, isInit0 = false )
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
	
	fModOut = fMod;
	fModOut = Utils.strAppendWith_( fModOut, getUpdaterFMod(updaterType) );
	if isInit0
		# if !isempty(fMod)
			# fModOut *= "_";
		# end
		# fModOut *= "isInit0";
		fModOut = Utils.strAppendWith_( fModOut, "isInit0" );
	end
	
	fMain = fMainLoopsMC;
	attrLst = cFerro == 0 ? attrLstLoops : attrLstLoopsFerro;
	# ["divNum","itNum","cArea","cPerim","beta"];
	valLst = Any[divNum, itNum, cArea, cPerim, beta];
	if cFerro != 0
		push!( valLst, cFerro );
	end
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
	for dim = 1 : params.nDim, pos in params.posLst
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
						linkFerroLst[lnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]] = !linkFerroLst[lnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]];
					end
				end
			end
			@infiltrate
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

function flipBLinkAtPos( params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}; pos::CartesianIndex{D}, dim::Int64 ) where {D}
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

function boolToOnePN( varBool::Bool )
	# return -(-1).^varBool;
	return varBool ? 1 : -1;
end

include("loops_MC_run_funcs.jl")
export runLoopMC_withParams

include("loops_MC_methods_seperate_funcs.jl")
export loops_MC, loops_MC_smart, loops_MC_staggeredCube

include("loops_MC_resave_funcs.jl")

include("loops_MC_old_funcs.jl")



end #endmodule
