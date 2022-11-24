module NelderMeadTest

using Utils
using Statistics
# using Infiltrator

export nmOpt!, cosXY

function cosXY( xyLst )
	return cos(xyLst[1]) + cos(xyLst[2]);
end

export nmArrsConstruct, NelderMeadArrs

struct NelderMeadArrs
	nDim::Int64;
	nPts::Int64;
	
	ptLstOrg::Array{Float64,2};
	ptLstTmp::Array{Float64,2};;
	ptLst::Array{Float64,2};
	valLstTmp::Vector{Float64};
	valLst::Vector{Float64};
	ixLst::Vector{Int64};
	
	volMat::Matrix{Float64};
	
	oPt::Vector{Float64};
	rPt::Vector{Float64};
	ePt::Vector{Float64};
	
	rVal::Base.RefValue{Float64};
	eVal::Base.RefValue{Float64};
	
	auxLst;
	auxLstTmp;
end

function nmArrsConstruct( nDim; auxLst = nothing, auxLstTmp = nothing )
	nPts = nDim+1;
	valLst = zeros(nPts);
	valLstTmp = zeros(nPts);
	ptLst = zeros(nDim, nPts);
	ptLstTmp = copy(ptLst);
	ptLstOrg = copy(ptLst);
	ixLst = zeros(Int64, nPts);
	
	volMat = zeros(nDim,nDim);
	
	oPt = zeros(nDim);
	rPt = zeros(nDim);
	ePt = zeros(nDim);
	
	rVal = Float64(0);
	eVal = Float64(0);
	
	NelderMeadArrs( nDim, nPts, ptLstOrg, ptLstTmp, ptLst, valLstTmp, valLst, ixLst, volMat, oPt, rPt, ePt, Ref(rVal), Ref(eVal), auxLst, auxLstTmp );
end

include("NelderMead_funcs_obj.jl")

function nmOpt!( funOrg, dim, ptLstOrg, valLst, ixLst, ptLst, ptLstTmp, valLstTmp, threshold = 1e-10, cntCut = 1000; alpha = 1, gamma = 2, rho = 1/2, sigma = 1/2, funInit! = nothing, isRecVals = true, yesAux = false, auxLst = nothing, auxLstTmp = nothing )
	if yesAux
		fun! = funOrg;
	else
		fun! = (pt,aux,auxOut)->( funOrg(pt) );
		funInit! = (pt,auxOut)->( funOrg(pt) );
		auxLst = similar( valLst, Nothing );
		auxLstTmp = similar( auxLst );
	end
	numPt = dim+1;
	ptLst .= ptLstOrg;
	valLstInit!( valLst, ptLst, funInit!, auxLst );
	rVal = 0;
	eVal = 0;
	oPt = zeros(dim);
	rPt = zeros(dim);
	ePt = zeros(dim);
	rAux = deepcopy( auxLst[1] );
	eAux = deepcopy( auxLst[1] );
	oAux = deepcopy( auxLst[1] );
	
	valVar = std( valLst );
	cnt = 0;
	if isRecVals
		valRecLst = [];
	end
	
	while valVar > threshold && cnt <= cntCut
		cnt += 1;
		println(cnt);
		# print( "\e[s", "count: ", cnt, ", val: ", valLst[1], "\e[u" );
		# print( "\e[s", "count: ", cnt, ", val: ", valLst[1], "\e[u" );
		sortperm!( ixLst, valLst );
		permute1d!( valLst, valLstTmp, ixLst );
		permute1d!( auxLst, auxLstTmp, ixLst );
		permuteCol2d!( ptLst, ptLstTmp, ixLst );
		
		for id = 1 : dim
			oPt[id] = sum( @view(ptLst[id,1:numPt-1]) ) / (numPt-1);
		end
		rPt .= oPt .+ alpha .* ( oPt .- @view(ptLst[:,end]) )
		rVal = fun!( rPt, auxLst[end], rAux );
		if rVal < valLst[1]
			ePt .= oPt .+ gamma .* ( rPt .- oPt );
			eVal = fun!( ePt, auxLst[end], eAux );
			if eVal < rVal
				updatePtVal!( ptLst, valLst, auxLst, ePt, eVal, eAux );
			else
				updatePtVal!( ptLst, valLst, auxLst, rPt, rVal, rAux );
			end
		elseif rVal >= valLst[numPt-1]
			if rVal <= valLst[end]
				ePt .= oPt .+ rho .* (rPt .- oPt);
				eVal = fun!(ePt, auxLst[end], eAux);
				if eVal < rVal
					updatePtVal!( ptLst, valLst, auxLst, ePt, eVal, eAux );
				else
					shrinkPtLst!( ptLst, valLst, sigma, fun!, auxLst, auxLstTmp );
				end
			else
				ePt .= oPt .+ rho .* ( @view(ptLst[:,end]) .- oPt );
				eVal = fun!( ePt, auxLst[end], eAux );
				if eVal < valLst[end]
					updatePtVal!( ptLst, valLst, auxLst, ePt, eVal, eAux );
				else
					shrinkPtLst!( ptLst, valLst, sigma, fun!, auxLst, auxLstTmp );
				end
			end
		else
			updatePtVal!( ptLst, valLst, auxLst, rPt, rVal, rAux );
		end
		valVar = std( valLst );
		if isRecVals
			push!( valRecLst, valLst[end] );
		end
	end
	return cnt;
	# @infiltrate
end

function updatePtVal!( ptLst, valLst, auxLst, pt, val, aux )
	@view(ptLst[:,end]) .= pt;
	valLst[end] = val;
	structAssign!( auxLst[end], aux );
	# auxLst[end] = aux;
end

function shrinkPtLst!( ptLst, valLst, sigma, fun!, auxLst, auxLstTmp )
	for i2 = 2 : size( ptLst, 2 )
		for i1 = 1 : size( ptLst, 1 )
			ptLst[i1,i2] = ptLst[i1,1] + sigma * ( ptLst[i1,i2] - ptLst[i1,1] );
		end
	end
	valLstRefresh!( valLst, ptLst, fun!, auxLst; auxLstTmp = auxLstTmp );
end

function valLstInit!( valLst, ptLst, fun!, auxLst )
	valLstRefresh!( valLst, ptLst, fun!, auxLst; isInit = true );
end

function valLstRefresh!( valLst, ptLst, fun!, auxLst; auxLstTmp = nothing, isInit = false )
	for iVal = 1 : length(valLst)
		if isInit
			valLst[iVal] = fun!( ptLst[:,iVal], auxLst[iVal] );
		else
			structAssign!( auxLstTmp[iVal], auxLst[iVal] );
			valLst[iVal] = fun!( ptLst[:,iVal], auxLstTmp[iVal], auxLst[iVal] );
		end
		# valLst[iVal] = ( isInit ? fun!( ptLst[:,iVal], auxLst[iVal] ) : fun!( ptLst[:,iVal], auxLst[iVal] ) );
	end
end

end
