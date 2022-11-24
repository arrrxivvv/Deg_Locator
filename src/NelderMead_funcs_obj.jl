using Statistics
using LinearAlgebra

function nmOpt!( nmArrs::NelderMeadArrs, funOrg; ptLstStart = nothing, thresVal, thresSz, cntCut = 1000, alpha = 1, gamma = 2, rho = 1/2, sigma = 1/2, funInit! = nothing, isRecVals = true )
	if !isnothing(ptLstStart)
		nmArrs.ptLstOrg .= ptLstStart; 
	end
	nmArrs.ptLst .= nmArrs.ptLstOrg;
	
	valLstRefresh!( funOrg, nmArrs );
	
	valVar = std(nmArrs.valLst);
	volume = volumeNmArrs!( nmArrs );
	volHedronThres = 0.1 * thresSz^nmArrs.nDim;
	
	cnt = 0;
	while ( valVar > thresVal || volume > volHedronThres ) && cnt < cntCut
		cnt += 1;
		print("\r","cnt = $cnt             ")
		
		sortperm!( nmArrs.ixLst, nmArrs.valLst );
		permute1d!( nmArrs.valLst, nmArrs.valLstTmp, nmArrs.ixLst );
		permuteCol2d!( nmArrs.ptLst, nmArrs.ptLstTmp, nmArrs.ixLst );
		
		for iDim = 1 : nmArrs.nDim
			nmArrs.oPt[iDim] = sum( @view( nmArrs.ptLst[iDim,1:nmArrs.nPts-1] ) ) / (nmArrs.nPts-1);
		end
		nmArrs.rPt .= nmArrs.oPt .+ alpha .* ( nmArrs.oPt .- @view(nmArrs.ptLst[:,end]) );
		nmArrs.rVal[] = funOrg( nmArrs.rPt );
		
		if nmArrs.rVal[] < nmArrs.valLst[1]
			nmArrs.ePt .= nmArrs.oPt .+ gamma .* (nmArrs.rPt - nmArrs.oPt); 
			nmArrs.eVal[] = funOrg( nmArrs.ePt );
			if nmArrs.eVal[] < nmArrs.rVal[]
				updateNmArrsWorst!( nmArrs, nmArrs.ePt, nmArrs.eVal[] );
			else
				updateNmArrsWorst!( nmArrs, nmArrs.rPt, nmArrs.rVal[] );
			end
		elseif nmArrs.rVal[] >= nmArrs.valLst[nmArrs.nPts-1]
			if nmArrs.rVal[] <= nmArrs.valLst[end]
				nmArrs.ePt .= nmArrs.oPt .+ rho .* ( nmArrs.rPt - nmArrs.oPt );
				nmArrs.eVal[] = funOrg( nmArrs.ePt );
				if nmArrs.eVal[] < nmArrs.rVal[]
					updateNmArrsWorst!( nmArrs, nmArrs.ePt, nmArrs.eVal[] );
				else
					shrinkNmArr!(nmArrs, sigma, funOrg);
				end
			else
				nmArrs.ePt .= nmArrs.oPt .+ rho .* ( @views(nmArrs.ptLst[:,end] .- nmArrs.oPt) );
				nmArrs.eVal[] = funOrg( nmArrs.ePt );
				if nmArrs.eVal[] < nmArrs.valLst[end]
					updateNmArrsWorst!( nmArrs, nmArrs.ePt, nmArrs.eVal[] );
				else
					shrinkNmArr!( nmArrs, sigma, funOrg );
				end
			end
		else
			updateNmArrsWorst!( nmArrs, nmArrs.rPt, nmArrs.rVal[] );
		end
		valVar = std(nmArrs.valLst);
		volume = volumeNmArrs!(nmArrs);
	end
	return cnt;
end

function valLstRefresh!( fun, nmArrs::NelderMeadArrs )
	for iPt = 1 : nmArrs.nPts
		nmArrs.valLst[iPt] = fun( @view( nmArrs.ptLst[:,iPt] ) );
	end
end

function shrinkNmArr!( nmArrs::NelderMeadArrs, sigma, fun )
	for iPt = 2 : nmArrs.nPts, iDim = 1 : nmArrs.nDim
		nmArrs.ptLst[iDim,iPt] = 
			sigma * (nmArrs.ptLst[iDim,iPt] - nmArrs.ptLst[iDim,1]) +
			nmArrs.ptLst[iDim,1];
	end
	valLstRefresh!( fun, nmArrs );
end

function updateNmArrsWorst!( nmArrs::NelderMeadArrs, pt, val )
	@view(nmArrs.ptLst[:,end]) .= pt;
	nmArrs.valLst[end] = val;
end

function volumeNmArrs!( nmArrs::NelderMeadArrs )
	for iPt = 2 : nmArrs.nPts
		@views nmArrs.volMat[:,iPt-1] .= nmArrs.ptLst[:,iPt] .- nmArrs.ptLst[:,1];
	end
	return det( nmArrs.volMat );
end
