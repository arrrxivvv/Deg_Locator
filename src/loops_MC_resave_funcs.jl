#Loops_MC funcs
#resave from files


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