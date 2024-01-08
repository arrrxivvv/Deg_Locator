module NelderMeadTest

using Utils
using Statistics
# using Infiltrator

export nmOpt!, cosXY

function cosXY( xyLst )
	return cos(xyLst[1]) + cos(xyLst[2]);
end

function nmOpt!( fun, dim, ptLstOrg, valLst, ixLst, ptLst, ptLstTmp, valLstTmp, threshold = 1e-10, cntCut = 1000; alpha = 1, gamma = 2, rho = 1/2, sigma = 1/2, isRecVals = true )
	numPt = dim+1;
	ptLst .= ptLstOrg;
	valLstRefresh!( valLst, ptLst, fun )
	rVal = 0;
	eVal = 0;
	oPt = zeros(dim);
	rPt = zeros(dim);
	ePt = zeros(dim);
	valVar = std( valLst );
	cnt = 0;
	if isRecVals
		valRecLst = [];
	end
	
	while valVar > threshold && cnt <= cntCut
		cnt += 1;
		# print( "\e[s", "count: ", cnt, ", val: ", valLst[1], "\e[u" );
		# print( "\e[s", "count: ", cnt, ", val: ", valLst[1], "\e[u" );
		sortperm!( ixLst, valLst );
		permute1d!( valLst, valLstTmp, ixLst );
		permuteCol2d!( ptLst, ptLstTmp, ixLst );
		# @infiltrate
		
		for id = 1 : dim
			oPt[id] = sum( @view(ptLst[id,1:numPt-1]) ) / (numPt-1);
		end
		rPt .= oPt .+ alpha .* ( oPt .- @view(ptLst[:,end]) )
		rVal = fun( rPt );
		if rVal < valLst[1]
			ePt .= oPt .+ gamma .* ( rPt .- oPt );
			eVal = fun( ePt );
			if eVal < rVal
				updatePtVal!( ptLst, valLst, ePt, eVal );
			else
				updatePtVal!( ptLst, valLst, rPt, rVal );
			end
		elseif rVal >= valLst[numPt-1]
			if rVal <= valLst[end]
				ePt .= oPt .+ rho .* (rPt .- oPt);
				eVal = fun(ePt);
				if eVal < rVal
					updatePtVal!( ptLst, valLst, ePt, eVal );
				else
					shrinkPtLst!( ptLst, valLst, sigma, fun );
				end
			else
				ePt .= oPt .+ rho .* ( @view(ptLst[:,end]) .- oPt );
				eVal = fun( ePt );
				if eVal < valLst[end]
					updatePtVal!( ptLst, valLst, ePt, eVal );
				else
					shrinkPtLst!( ptLst, valLst, sigma, fun );
				end
			end
		else
			updatePtVal!( ptLst, valLst, rPt, rVal );
			# @view(ptLst[:,end]) .= rPt;
			# valLst[end] = rVal;
		end
		valVar = std( valLst );
		# @infiltrate cond = cnt >= 500;
		if isRecVals
			push!( valRecLst, valLst[end] );
		end
	end
	# @infiltrate # cond = valLst[1] <= 1e3 *threshold
	return cnt;
end

function updatePtVal!( ptLst, valLst, pt, val )
	@view(ptLst[:,end]) .= pt;
	valLst[end] = val;
end

function shrinkPtLst!( ptLst, valLst, sigma, fun )
	for i2 = 2 : size( ptLst, 2 )
		for i1 = 1 : size( ptLst, 1 )
			ptLst[i1,i2] = ptLst[i1,1] + sigma * ( ptLst[i1,i2] - ptLst[i1,1] );
		end
	end
	valLstRefresh!( valLst, ptLst, fun )
end

function valLstRefresh!( valLst, ptLst, fun )
	for iVal = 1 : length(valLst)
		valLst[iVal] = fun( ptLst[:,iVal] );
	end
end

end
