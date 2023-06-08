module Loops_MC

using ShiftedArrays
using FilenameManip
using Random
using JLD2

# using Infiltrator

export loops_MC

function loops_MC( divNum = 64, itNum = 10000; fMod = "", cArea = 1, cPerim = 1, beta = 1 )
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

	idLayerLst = zeros( UInt, divNum );
	dELst = zeros(divNum);

	# cArea = 1;
	# cPerim = 1;
	# beta = 1;
	
	rejLst = zeros(divNum);
	for it = 1 : itNum
		print( "it = ", it, "         \r" )
		for dim = 1 : nDim
			@view(zakLstLst[:,:,dim,it]) .= dropdims( reduce( xor, BfieldLst[dim]; dims = dim, init = false ); dims = dim );
			
			for pos in posLst
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
		end
		for dim = 1 : nDim
			rand!( idLayerLst );
			idLayerLst .= mod.( idLayerLst, lnLayer ) .+ 1;
			dELst .= 0;
			rand!(rejLst);
			# Threads.@threads 
			for iLayer = 1 : divNum
				dELst[iLayer] -= cArea * boolToOnePN( BfieldLst[dim][posLstDimLst[dim][idLayerLst[iLayer]]] );
				for dimLink in linkDimLst[dim]
					dELst[iLayer] -= cPerim * boolToOnePN( linkLst[dimLink][posLstDimLst[dim][idLayerLst[iLayer]]] );
					for dimLinkSh in linkDimLst[dim]
						if dimLinkSh != dimLink
							dELst[iLayer] -= cPerim * boolToOnePN( linkLst[dimLink][posLstShLst[dimLinkSh,1][posLstDimLst[dim][idLayerLst[iLayer]]]] );
						end
					end
				end
				# @infiltrate
				if dELst[iLayer] < 0 || rejLst[iLayer] < exp( - beta * dELst[iLayer] )
					BfieldLst[dim][posLstDimLst[dim][idLayerLst[iLayer]]] = !BfieldLst[dim][posLstDimLst[dim][idLayerLst[iLayer]]];
					for dimLink in linkDimLst[dim]
						linkLst[dimLink][posLstDimLst[dim][idLayerLst[iLayer]]] = !linkLst[dimLink][posLstDimLst[dim][idLayerLst[iLayer]]];
						for dimLinkSh in linkDimLst[dim]
							if dimLinkSh != dimLink
								linkLst[dimLink][posLstShLst[dimLinkSh,1][posLstDimLst[dim][idLayerLst[iLayer]]]] = !linkLst[dimLink][posLstShLst[dimLinkSh,1][posLstDimLst[dim][idLayerLst[iLayer]]]];
							end
						end
					end
					# @infiltrate
				end
				# @infiltrate
			end
		end
	end
	
	fMain = "loops_MC";
	attrLst = ["divNum","itNum","cArea","cPerim","beta"];
	valLst = [Int64(divNum),Int64(itNum), cArea, cPerim,beta];
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	save( fName, "zakLstLst", zakLstLst, "divNum", divNum, "itNum", itNum, "cArea", cArea, "cPerim", cPerim, "beta", beta );
end

function MCLinkUpdate( BfieldLst, linkLst,  )
	for dim = 1 : nDim
	end
end

function boolToOnePN( varBool::Bool )
	return -(-1).^varBool;
end

end #endmodule
