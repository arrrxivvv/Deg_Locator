# Module Loops_MC

function getAttrValLstLttcLst( divNumLst::Vector{Int64}, itNumLst::Vector{Int64}, nDim::Int64 )
	attrLstFLst = getAttrLstLttcBase();
	valLstFLst = Any[];
	append!( valLstFLst, summarizeArrAttr.( [divNumLst, itNumLst] ), [nDim] );
	
	return attrLstFLst, valLstFLst;
end

function getAttrValLttcGroup( divNumLst::Vector{Int64}, itNumLst::Vector{Int64}, nDim::Int64,  paramsGrp::ParamsGroup; rndDigits = 3 )
	attrLstFLst, valLstFLst = getAttrValLstLttcLst( divNumLst, itNumLst, nDim );
	append!( attrLstFLst, getGroupAttr( paramsGrp ) );
	append!( valLstFLst, getGroupValLst( paramsGrp; digs = rndDigits ) );
	
	return attrLstFLst, valLstFLst;
end

function getAttrValGroupsSummarized( divNumLst::Vector{Int64}, itNumLst::Vector{Int64}, nDim::Int64, paramsGrpLst::Vector{ParamsGroup}; rndDigits = 3 )
	attrLstFLst, valLstFLst = getAttrValLstLttcLst( divNumLst, itNumLst, nDim );
	for iGrp = 1 : length(paramsGrpLst)
		append!( attrLstFLst, getGroupAttr( paramsGrpLst[iGrp] ) );
		append!( valLstFLst, getGroupValLst( paramsGrpLst[iGrp]; digs = rndDigits ) );
	end
	
	return attrLstFLst, valLstFLst;
end

function getAttrValGroupsItNumLstSummarized( divNumLst::Vector{Int64}, itNumLstLst::Vector{Vector{Int64}}, nDim::Int64, paramsGrpLstLst::Vector{Vector{ParamsGroup}}; isAbbrev = false, rndDigits = 3 )
	itNumLst = deepcopy(itNumLstLst);
	itNumLst = append!(itNumLst...);
	paramsGrpLst = deepcopy(paramsGrpLstLst);
	paramsGrpLst = append!( paramsGrpLst... );
	attrLstFLst, valLstFLst = getAttrValLstLttcLst( divNumLst, itNumLst, nDim );
	for iGrp = 1 : length(paramsGrpLst)
		attrLstTmp = getGroupAttr( paramsGrpLst[iGrp] );
		valLstTmp = getGroupValLst( paramsGrpLst[iGrp]; digs = rndDigits );
		if isAbbrev
			attrLstTmp = attrLstTmp[1:1];
			valLstTmp = valLstTmp[1:1];
		end
		append!( attrLstFLst, attrLstTmp );
		append!( valLstFLst, valLstTmp );
	end
	
	return attrLstFLst, valLstFLst;
end

function runLoopMC_withParamsGroup( updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, itNumLstLst::Vector{Vector{Int64}}, divNumLst::Vector{Int64}, paramsGrpLstLst::Vector{Vector{ParamsGroup}}; fMod = "", nDim::Int64 = 3 )
	if length(itNumLstLst) != length(paramsGrpLstLst)
		throw(ArgumentError("itNumLstLst and paramsGrpLstLst lengths different"));
	end
	
	for iItLst = 1 : length(itNumLstLst)
		runLoopMC_withParamsGroup( updaterType, initializerType, itNumLstLst[iItLst], divNumLst,  paramsGrpLstLst[iItLst] );
	end
end

function runLoopMC_withParamsGroup( updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, paramsGrpLst::Vector{ParamsGroup}; fMod = "", nDim::Int64 = 3 )
	for iGrp = 1 : length(paramsGrpLst)
		cAreaLst, cPerimLst, cFerroLst = paramsLstFromGroup( paramsGrpLst[iGrp] );
		runLoopMC_withParamsBase( updaterType, initializerType, itNumLst, divNumLst, cAreaLst, cPerimLst, cFerroLst; fMod = fMod, nDim );
	end
end

function runLoopMC_withParamsBase( updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, cAreaLst::Array{Float64}, cPerimLst::Array{Float64}, cFerroLst::Array{Float64}; fMod = "", nDim::Int64 = 3 )
	@time for itNum in itNumLst, divNum in divNumLst
		for idParam in CartesianIndices(cAreaLst)
			cArea = cAreaLst[idParam];
			cPerim = cPerimLst[idParam];
			cFerro = cFerroLst[idParam];
			
			@time fNameSmart = loops_MC_methods_cALF( divNum, itNum; nDim = nDim, cPerim = cPerim, cArea = cArea, cFerro = cFerro, fMod = fMod, updaterType = updaterType, initializerType = initializerType );
			GC.gc()
		end
	end
end

function genFNameToLogs( fNameLst::Vector{String}, fLstName::String, fArrLstName::String )
	writedlm( fLstName, fNameLst );
	
	open(dirLog * fNameFileLstLst, "w") do io
		println(io, fLstName);
	end
	open( dirLog * fNameFileLstJld2Lst, "w" ) do io
		println(io, fArrLstName);
	end
end

function genFNameLstLoopMC( updaterType::Type{<:LoopsUpdater}, itNumLstLst::Vector{Vector{Int64}}, divNumLst::Vector{Int64}, paramsGrpLstLst::Vector{Vector{ParamsGroup}}; fMain::String, fMod = "", rndDigits::Int64 = 3, isAbbrev = true, nDim::Int64 = 3 )
	lnItNums = length(itNumLstLst);
	if lnItNums != length(paramsGrpLstLst)
		throw(ArgumentError("itNumLstLst and paramsGrpLstLst lengths different"));
	end	
	
	fNameLst, fItGrpLstName = genFNameItGrpLstArrSaved( updaterType, itNumLstLst, divNumLst, paramsGrpLstLst; fMain = fMain, fMod = fMod, rndDigits = rndDigits, isAbbrev = isAbbrev, nDim = nDim );
	
	fLstMain = fMain * "_fLst";
	attrLstFLst, valLstFLst = getAttrValGroupsItNumLstSummarized( divNumLst, itNumLstLst, nDim, paramsGrpLstLst; isAbbrev = isAbbrev );
	fModFLst = getFModLoopsMC( fMod, updaterType );
	fLstName = fNameFunc( fLstMain, attrLstFLst, valLstFLst, ".txt"; fMod = fModFLst );
	genFNameToLogs( fNameLst, fLstName, fItGrpLstName );
	
	return fLstName;
end

function genFNameLstLoopMC( updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, paramsGrpLst::Vector{ParamsGroup}; fMain::String, fMod = "", rndDigits::Int64 = 3, nDim::Int64 = 3 )
	fNameLst, fGrpLstName = genFNameGrpLstArrSaved( updaterType, initializerType, itNumLst, divNumLst, paramsGrpLst; fMain = fMain, fMod = fMod, rndDigits = rndDigits, nDim = nDim );
	
	fLstMain = fMain * "_fLst";
	attrLstFLst, valLstFLst = getAttrValGroupsSummarized( divNumLst, itNumLst, nDim, paramsGrpLst; rndDigits = rndDigits );
	fModFLst = getFModLoopsMC( fMod, updaterType );
	fLstName = fNameFunc( fLstMain, attrLstFLst, valLstFLst, ".txt"; fMod = fModFLst );
	genFNameToLogs( fNameLst, fLstName, fGrpLstName );
	
	return fLstName;
end

function genFNameItGrpLstArrSaved( updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, itNumLstLst::Vector{Vector{Int64}}, divNumLst::Vector{Int64}, paramsGrpLstLst::Vector{Vector{ParamsGroup}}; fMain::String, fMod = "", rndDigits::Int64 = 3, isAbbrev = true, nDim::Int64 = 3 )
	fNameLstLst = Vector{Vector{String}}(undef, length(itNumLstLst));
	fNameJld2Lst = Vector{String}(undef,length(itNumLstLst));
	for iIt = 1 : length(itNumLstLst)
		fNameLstLst[iIt], fNameJld2Lst[iIt] = genFNameGrpLstArrSaved( updaterType, initializerType, itNumLstLst[iIt], divNumLst, paramsGrpLstLst[iIt]; fMain = fMain, fMod = fMod, rndDigits = rndDigits, nDim = nDim );
	end
	fNameLst = append!(fNameLstLst...);
	
	fLstMain = fMain * "_fLst";
	fItGrpLstMain = fMain * "_fArr" * "ItLstGrp";
	attrFLstItLstGrp, valFLstItLstGrp = getAttrValGroupsItNumLstSummarized( divNumLst, itNumLstLst, nDim, paramsGrpLstLst; isAbbrev = isAbbrev, rndDigits = rndDigits );
	fModFull = getFModLoopsMC( fMod, updaterType );
	fItGrpLstName = fNameFunc( fItGrpLstMain, attrFLstItLstGrp, valFLstItLstGrp, jld2Type; fMod = fModFull );
	
	save( fItGrpLstName, "fLstJld2MasterNameLst", fNameJld2Lst );
	
	return fNameLst, fItGrpLstName;
end

function genFNameGrpLstArrSaved( updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, paramsGrpLst::Vector{ParamsGroup}; fMain::String, fMod = "", rndDigits::Int64 = 3, nDim::Int64 = 3 )
	fNameLstLst = Vector{Vector{String}}(undef, length(paramsGrpLst));
	fNameJld2Lst = Vector{String}(undef,length(paramsGrpLst));
	for iGrp = 1 : length(paramsGrpLst)
		fNameLstLst[iGrp], fNameJld2Lst[iGrp] = genFNameLstArrSaved( updaterType, initializerType, itNumLst, divNumLst, paramsGrpLst[iGrp]; fMain = fMain, fMod = fMod, rndDigits = rndDigits, nDim = nDim );
	end
	fNameLst = append!(fNameLstLst...);
	
	fLstMain = fMain * "_fLst";
	fGrpLstMain = fMain * "_fArr" * "Grp";
	attrFLstGrp, valFLstGrp = getAttrValGroupsSummarized( divNumLst, itNumLst, nDim, paramsGrpLst; rndDigits = rndDigits );
	fModFull = getFModLoopsMC( fMod, updaterType );
	fGrpLstName = fNameFunc( fGrpLstMain, attrFLstGrp, valFLstGrp, jld2Type; fMod = fModFull );
	
	save( fGrpLstName, "fLstNameLst", fNameJld2Lst );
	
	return fNameLst, fGrpLstName;
end

function genFNameLstArrSaved( updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, paramsGrp::ParamsGroup; fMain::String, fMod = "", rndDigits::Int64 = 3, nDim::Int64 = 3 )
	cAreaLst, cPerimLst, cFerroLst = paramsGrp.paramsLst;
	fNameLst, fNameArr = genFNameLstInJulia( updaterType, initializerType, itNumLst, divNumLst, cAreaLst, cPerimLst, cFerroLst; fMain = fMain, fMod = fMod, nDim = nDim );
	
	fLstMain = fMain * "_fLst";
	fLstJld2Main = fMain * "_fArr" * "Jld2";
	attrLstFLst, valLstFLst = getAttrValLttcGroup( divNumLst, itNumLst, nDim, paramsGrp; rndDigits = rndDigits );
	fModFull = getFModLoopsMC( fMod, updaterType );
	fLstJld2Name = fNameFunc( fLstJld2Main, attrLstFLst, valLstFLst, jld2Type; fMod = fModFull );
	save(fLstJld2Name, "fNameArr", fNameArr);
	
	return fNameLst, fLstJld2Name;
end

function genFNameLstInJulia( updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, cAreaLst::Array{Float64}, cPerimLst::Array{Float64}, cFerroLst::Array{Float64}; fMain::String, fMod = "", nDim::Int64 = 3 )
	isFileNameOnly = true;
	fNameLst = Vector{String}(undef,0);
	fNameArr = Array{String}(undef, length(itNumLst), length(divNumLst), size(cAreaLst)...);
	@time for iIt = 1 : length(itNumLst), iDiv = 1 : length(divNumLst), idParam in CartesianIndices(cAreaLst)
		cArea = cAreaLst[idParam];
		cPerim = cPerimLst[idParam];
		cFerro = cFerroLst[idParam];
		divNum = divNumLst[iDiv];
		itNum = itNumLst[iIt];
		
		fName = loops_MC_methods_cALF( divNum, itNum; nDim = nDim, cPerim = cPerim, cArea = cArea, cFerro = cFerro, fMod = fMod, updaterType = updaterType, initializerType = initializerType, isFileNameOnly = isFileNameOnly, fMainOutside = fMain );
		GC.gc()
		push!(fNameLst,fName);
		fNameArr[iIt, iDiv, idParam] = fName;
	end
	
	return fNameLst, fNameArr;
end

function saveParamsLoopMC( updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, paramsGrpLst::Vector{ParamsGroup}; fMod = "", rndDigits = 3, nDim::Int64 = 3 )
	cAreaLst, cPerimLst, cFerroLst = getParamLstFromGroupLst( paramsGrpLst );
	cAreaLstStr = (a->string.(a)).(cAreaLst);
	cPerimLstStr = (a->string.(a)).(cPerimLst);
	cFerroLstStr = (a->string.(a)).(cFerroLst);
	
	fParamMain = "loopsMC_params";
	attrLstParamLst, valLstParamLst = getAttrValGroupsSummarized( divNumLst, itNumLst, nDim, paramsGrpLst; rndDigits = rndDigits );
	fModFull = getFModLoopsMC( fMod, updaterType );
	
	grpNameLst = getGroupName.( paramsGrpLst );
	grpParamsNameLst = getGroupAttr.( paramsGrpLst );
	grpParamsLst = getGroupParamsLst.( paramsGrpLst );
	
	fParamName = fNameFunc( fParamMain, attrLstParamLst, valLstParamLst, jld2Type; fMod = fModFull );
	save( fParamName, "itNumLst", itNumLst, "divNumLst", divNumLst, "nDim", nDim, "grpNameLst", grpNameLst, "grpParamsNameLst", grpParamsNameLst, "grpParamsLst", grpParamsLst, "cAreaLst", cAreaLst, "cAreaLstStr", cAreaLstStr, "cPerimLst", cPerimLst, "cPerimLstStr", cPerimLstStr, "cFerroLst", cFerroLst, "cFerroLstStr", cFerroLstStr );
	
	open( dirLog * fNameSaveParamsLst, "w" ) do io
		println( io, fParamName );
	end
	
	return fParamName;
end

function saveParamsLoopMC( updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, itNumLstLst::Vector{Vector{Int64}}, divNumLst::Vector{Int64}, paramsGrpLstLst::Vector{Vector{ParamsGroup}}; fMod = "", rndDigits = 3, isAbbrev = true, nDim::Int64 = 3 )
	cParamsLstLst = getParamLstFromGroupLst.( paramsGrpLstLst );
	cAreaLst = [ cParamsLstLst[iIt][1] for iIt = 1 : length(itNumLstLst) ];
	cPerimLst = [ cParamsLstLst[iIt][2] for iIt = 1 : length(itNumLstLst) ];
	cFerroLst = [ cParamsLstLst[iIt][3] for iIt = 1 : length(itNumLstLst) ];
	toString2Depth = b -> (a->string.(a)).(b);
	cAreaLstStr = toString2Depth.(cAreaLst);
	cPerimLstStr = toString2Depth.(cPerimLst);
	cFerroLstStr = toString2Depth.(cFerroLst);
	
	fParamMain = "loopsMC_params";
	attrLstParamLst, valLstParamLst = getAttrValGroupsItNumLstSummarized( divNumLst, itNumLstLst, nDim, paramsGrpLstLst; isAbbrev = isAbbrev, rndDigits = rndDigits );
	fModFull = getFModLoopsMC( fMod, updaterType );
	
	grpNameLst = ( (a)->getGroupName.(a) ).( paramsGrpLstLst );
	grpParamsNameLst = ( (a)->getGroupAttr.(a) ).( paramsGrpLstLst );
	grpParamsLst = (  (a)->getGroupParamsLst.(a) ).( paramsGrpLstLst );
	
	fParamName = fNameFunc( fParamMain, attrLstParamLst, valLstParamLst, jld2Type; fMod = fModFull );
	save( fParamName, "itNumLst", itNumLstLst, "divNumLst", divNumLst, "nDim", nDim, "grpNameLst", grpNameLst, "grpParamsNameLst", grpParamsNameLst, "grpParamsLst", grpParamsLst, "cAreaLst", cAreaLst, "cAreaLstStr", cAreaLstStr, "cPerimLst", cPerimLst, "cPerimLstStr", cPerimLstStr, "cFerroLst", cFerroLst, "cFerroLstStr", cFerroLstStr );
	
	open( dirLog * fNameSaveParamsLst, "w" ) do io
		println( io, fParamName );
	end
	
	return fParamName;
end
