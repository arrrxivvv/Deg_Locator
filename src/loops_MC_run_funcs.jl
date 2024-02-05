#Loops_MC funcs
#run with different params
using DelimitedFiles

fNameFileLstLst = "fNameFileLstLst.txt";
fNameFileLstJld2Lst = "fNameFileLstJld2Lst.txt";
fNameSaveParamsLst = "fNameSaveParamsLst.txt";

function runLoopMC_withParams( updaterType::(Type{T} where T <: LoopsUpdater), itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, betaLst::Vector{Float64}, cRatioLst::Vector{Float64}, cFerroRatioLst::Vector, sgnAreaLst::Vector{Int64}, sgnPerimLst::Vector{Int64}, isInit0Lst::Vector{Bool}; fMod = "", betaBase = 1 )
	@time for itNum in itNumLst, divNum in divNumLst, beta in betaLst, cRatio in cRatioLst, cFerroRatio in cFerroRatioLst, sgnArea in sgnAreaLst, sgnPerim in sgnPerimLst, isInit0 in isInit0Lst
		cRatioSq = sqrt(cRatio);
		cArea = sgnArea * beta * cRatioSq;
		cPerimAbs = beta / cRatioSq;
		cPerim = sgnPerim * cPerimAbs;
		cFerro = cPerimAbs * cFerroRatio;
		# @time fNameSmart = methodLoops( divNum, itNum; cPerim = cPerim, cArea = cArea, cFerro = cFerro, beta = betaBase, fMod = fModSmart, isInit0 = isInit0 );
		@time fNameSmart = loops_MC_methods( divNum, itNum; cPerim = cPerim, cArea = cArea, cFerro = cFerro, beta = betaBase, fMod = fMod, isInit0 = isInit0, updaterType = updaterType );
		GC.gc()
	end
end

function genFNameLstLoopMC( updaterType::(Type{T} where T <: LoopsUpdater), itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, betaLst::Vector{Float64}, cRatioLst::Vector{Float64}, cFerroRatioLst::Vector, sgnAreaLst::Vector{Int64}, sgnPerimLst::Vector{Int64}, isInit0Lst::Vector{Bool}; fMain::String, fMod = "", betaBase = 1, isFModMethod = true, rndDigits = 3 )
	fModWithMethod = fMod;
	if isFModMethod
		fModWithMethod = Utils.strAppendWith_( fMod, getUpdaterFMod( updaterType ) )
	end
	fNameLst = Vector{String}(undef,0);
	fNameArr = Array{String}(undef, length(itNumLst), length(divNumLst), length(betaLst), length(cRatioLst), length(cFerroRatioLst), length(sgnAreaLst), length(sgnPerimLst), length(isInit0Lst));
	for iIt = 1 : length(itNumLst), iDiv = 1 : length(divNumLst), iBeta = 1 : length(betaLst), iRatio = 1 : length(cRatioLst), iFerroRatio = 1 : length(cFerroRatioLst), iSgnArea = 1 : length(sgnAreaLst), iSgnPerim = 1 : length(sgnPerimLst), iInit0 = 1 : length(isInit0Lst)
		itNum = itNumLst[iIt];
		divNum = divNumLst[iDiv];
		cRatio = cRatioLst[iRatio];
		beta = betaLst[iBeta];
		cFerroRatio = cFerroRatioLst[iFerroRatio];
		sgnArea = sgnAreaLst[iSgnArea];
		sgnPerim = sgnPerimLst[iSgnPerim];
		isInit0 = isInit0Lst[iInit0];
		cRatioSq = sqrt(cRatio);
		cArea = sgnArea * beta * cRatioSq;
		cPerimAbs = beta / cRatioSq;
		cPerim = sgnPerim * cPerimAbs;
		cFerro = cPerimAbs * cFerroRatio;
		
		attrLst = cFerro == 0 ? attrLstLoops : attrLstLoopsFerro;
		valLst = Any[divNum, itNum, cArea, cPerim, betaBase];
		if cFerro != 0
			push!( valLst, cFerro );
		end
		fModFull = fModWithMethod;
		if isInit0
			fModFull = Utils.strAppendWith_( fModWithMethod, "isInit0" );
		end
		fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModFull );
		
		push!( fNameLst, fName );
		fNameArr[iIt,iDiv,iBeta,iRatio,iFerroRatio,iSgnArea,iSgnPerim,iInit0] = fName;
		GC.gc()
	end
	
	fLstMain = fMain * "_fileLst";
	fLstJld2Main = fLstMain * "Jld2";
	attrLstFLst = ["div", "it", "beta", "cRatio", "cFerroRatio", "isInit0", "sgnArea", "sgnPerim"];
	valLstFLst = Any[divNumLst[[1,end]],itNumLst[[1,end]],betaLst[[1,end]],round.(cRatioLst[[1,end]]; digits=rndDigits),cFerroRatioLst[[1,end]], isInit0Lst[[1,end]], sgnAreaLst[[1,end]], sgnPerimLst[[1,end]]];
	fLstName = fNameFunc( fLstMain, attrLstFLst, valLstFLst, ".txt"; fMod = fModWithMethod );
	fLstJld2Name = fNameFunc( fLstJld2Main, attrLstFLst, valLstFLst, jld2Type; fMod = fModWithMethod );
	writedlm( fLstName, fNameLst );
	save( fLstJld2Name, "fNameArr", fNameArr );
	
	open(dirLog * fNameFileLstLst, "w") do io
		println(io, fLstName);
	end
	open( dirLog * fNameFileLstJld2Lst, "w" ) do io
		println(io, fLstJld2Name);
	end
	# println( open(fNameFileLstLst), fLstName );	
	
	return fLstName;
end

function saveParamsLoopMC( updaterType::(Type{T} where T <: LoopsUpdater), itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, betaLst::Vector{Float64}, cRatioLst::Vector{Float64}, cFerroRatioLst::Vector, sgnAreaLst::Vector{Int64}, sgnPerimLst::Vector{Int64}, isInit0Lst::Vector{Bool}; fMod = "", rndDigits = 3, isFModMethod = true )
	fModWithMethod = fMod;
	if isFModMethod
		fModWithMethod = Utils.strAppendWith_( fMod, getUpdaterFMod( updaterType ) )
	end
	cAreaLst = zeros( length(betaLst), length(cRatioLst), length(cFerroRatioLst), length(sgnAreaLst), length(sgnPerimLst), length(isInit0Lst) );
	cPerimLst = similar(cAreaLst);
	cFerroLst = similar(cAreaLst);
	for iIt = 1 : length(itNumLst), iDiv = 1 : length(divNumLst), iBeta = 1 : length(betaLst), iRatio = 1 : length(cRatioLst), iFerroRatio = 1 : length(cFerroRatioLst), iSgnArea = 1 : length(sgnAreaLst), iSgnPerim = 1 : length(sgnPerimLst), iInit0 = 1 : length(isInit0Lst)
		itNum = itNumLst[iIt];
		cRatio = cRatioLst[iRatio];
		beta = betaLst[iBeta];
		cFerroRatio = cFerroRatioLst[iFerroRatio];
		sgnArea = sgnAreaLst[iSgnArea];
		sgnPerim = sgnPerimLst[iSgnPerim];
		cRatioSq = sqrt(cRatio);
		cArea = sgnArea * beta * cRatioSq;
		cPerimAbs = beta / cRatioSq;
		cPerim = sgnPerim * cPerimAbs;
		cFerro = cPerimAbs * cFerroRatio;
		
		cAreaLst[iBeta, iRatio, iFerroRatio, iInit0, iSgnArea, iSgnPerim] = cArea;
		cPerimLst[iBeta, iRatio, iFerroRatio, iInit0, iSgnArea, iSgnPerim] = cPerim;
		cFerroLst[iBeta, iRatio, iFerroRatio, iInit0, iSgnArea, iSgnPerim] = cFerro;
	end
	
	cAreaLstStr = string.(cAreaLst);
	cPerimLstStr = string.(cPerimLst);
	cFerroLstStr = string.(cFerroLst);
	
	fParamMain = "loopsMC_params";
	attrLstParamLst = ["beta", "cRatio", "cFerroRatio", "isInit0", "sgnArea", "sgnPerim"];
	valLstParamLst = Any[betaLst[[1,end]],round.(cRatioLst[[1,end]]; digits=rndDigits),cFerroRatioLst[[1,end]], isInit0Lst[[1,end]], sgnAreaLst[[1,end]], sgnPerimLst[[1,end]]];
	fParamName = fNameFunc( fParamMain, attrLstParamLst, valLstParamLst, jld2Type; fMod = fModWithMethod );
	save( fParamName, "divNumLst", divNumLst, "itNumLst", itNumLst, "betaLst", betaLst, "cRatioLst", cRatioLst, "cFerroRatioLst", cFerroRatioLst, "sgnAreaLst", sgnAreaLst, "sgnPerimLst", sgnPerimLst, "isInit0Lst", isInit0Lst, "cAreaLst", cAreaLst, "cAreaLstStr", cAreaLstStr, "cPerimLst", cPerimLst, "cPerimLstStr", cPerimLstStr, "cFerroLst", cFerroLst, "cFerroLstStr", cFerroLstStr );
	
	open( dirLog * fNameSaveParamsLst, "w" ) do io
		println( io, fParamName );
	end
	
	return fParamName;
end
