module Loops_MC

using ShiftedArrays
using FilenameManip
using Random
using JLD2
using Statistics

# using Infiltrator

export loops_MC

fMainLoopsMC = "loops_MC";
attrLstLoops = ["divNum","itNum","cArea","cPerim","beta"];

function loops_MC( divNum = 64, itNum = 10000; fMod = "", cArea = 1, cPerim = 1, beta = 1, itNumSample = 100 )
	nDim = 3;
	divLst = fill(div, nDim);
	posLst = CartesianIndices((divNum,divNum,divNum));
	posLstShLst = [ ShiftedArrays.circshift( posLst, ntuple( ( i -> ( i == dim ? (-1)^iPol : 0 ) ), nDim ) ) for dim = 1 : nDim, iPol = 1 : 2 ];
	
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

	BfieldLst = [ rand( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ];
	linkLst = [ zeros( Bool, divNum, divNum, divNum ) for dim = 1 : nDim ];
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
	
	numBfieldLst = zeros(Int64, itNum, nDim);
	numLinkLst = zeros(Int64, itNum, nDim);

	idLayerLst = zeros( UInt, divNum );
	dELst = zeros(divNum);
	
	rejLst = zeros(divNum);
	for it = 1 : itNum
		print( "it = ", it, "         \r" )
		for dim = 1 : nDim
			@view(zakLstLst[:,:,dim,it]) .= dropdims( reduce( xor, BfieldLst[dim]; dims = dim, init = false ); dims = dim );
			
			# @infiltrate
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
	
	itStep = Int64( floor(itNum / itNumSample) );
	zakLstSampleLst = zakLstLst[:,:,:,[itStep:itStep:itNum;]];
	
	fMain = fMainLoopsMC;
	attrLst = attrLstLoops;
	# ["divNum","itNum","cArea","cPerim","beta"];
	valLst = Any[divNum,itNum, cArea, cPerim,beta];
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	save( fName, "zakLstLst", zakLstLst, "divNum", divNum, "itNum", itNum, "cArea", cArea, "cPerim", cPerim, "beta", beta, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	
	return fName;
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
