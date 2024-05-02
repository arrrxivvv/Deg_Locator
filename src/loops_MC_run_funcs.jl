#Loops_MC funcs
#run with different params
using DelimitedFiles

fNameFileLstLst = "fNameFileLstLst.txt";
fNameFileLstJld2Lst = "fNameFileLstJld2Lst.txt";
fNameSaveParamsLst = "fNameSaveParamsLst.txt";

abstract type ParamsGroup end

errMsgParamsGrpUndefined = "params group not implemented";

function getGroupName( typeGrp::Type{<:ParamsGroup} )
	error( errMsgParamsGrpUndefined );
end

getGroupName( grp::ParamsGroup ) = getGroupName( typeof(grp) );
# function getGroupName( grp::ParamsGroup )
	# return getGroupName( type(grp) );
# end

function getGroupAttr( typeGrp::Type{<:ParamsGroup} )
	error(errMsgParamsGrpUndefined);
end

getGroupAttr( grp::ParamsGroup ) = getGroupAttr( typeof(grp) );

# function getGroupAttr( grp::ParamsGroup )
	# return getGroupAttr( type(grp) );
# end

function getGroupValLst( grp::ParamsGroup )
	error(errMsgParamsGrpUndefined);
end

function getGroupParamsLst( grp::ParamsGroup )
	return grp.groupParamsLst;
	# error(errMsgParamsGrpUndefined);
end

function paramsLstFromGroup( grp::ParamsGroup )
	# error(errMsgParamsGrpUndefined);
	cAreaLst = grp.paramsLst[1];
	cPerimLst = grp.paramsLst[2];
	cFerroLst = grp.paramsLst[3];
	
	return cAreaLst, cPerimLst, cFerroLst;
end

struct BetaParamsGroup <: ParamsGroup
	betaLst::Vector{Float64};
	cRatioLst::Vector{Float64};
	cFerroRatioLst::Vector{Float64};
	sgnAreaLst::Vector{Int64};
	sgnPerimLst::Vector{Int64};
	
	groupParamsLst::Vector{<:Vector};
	paramsLst::Vector{<:Array{Float64}};
	
	function BetaParamsGroup( betaLst, cRatioLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst )
		 groupParamsLst = [betaLst, cRatioLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst];
		 cAreaLst, cPerimLst, cFerroLst = paramsLstFromBetaRatio( betaLst, cRatioLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst );
		 paramsLst = [cAreaLst, cPerimLst, cFerroLst];
		 new( betaLst, cRatioLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst, groupParamsLst, paramsLst );
	end
end

function getGroupName( grpType::Type{BetaParamsGroup} )
	return "BetaCRatioGroup";
end

function getGroupAttr( groupType::Type{BetaParamsGroup} )
	attrLst = ["beta", "cRatio", "cFerroRatio", "sgnArea", "sgnPerim"];
	return attrLst;
end

function getGroupValLst( grp::BetaParamsGroup; digs = 3 )
	valLst = Any[];
	append!( valLst, summarizeArrAttrRound.( [grp.betaLst, grp.cRatioLst, grp.cFerroRatioLst]; digs = digs ), summarizeArrBoolAttr.([grp.sgnAreaLst, grp.sgnPerimLst]) );
	return valLst;
end

# function paramsLstFromGroup( grp::BetaParamsGroup )
	# betaLst = grp.betaLst;
	# cRatioLst = grp.cRatioLst;
	# cFerroRatioLst = grp.cFerroRatioLst;
	# sgnAreaLst = grp.sgnAreaLst;
	# sgnPerimLst = grp.sgnPerimLst;
	# cAreaLst = zeros( length(betaLst), length(cRatioLst), length(cFerroRatioLst), length(sgnAreaLst), length(sgnPerimLst) );
	# cPerimLst = similar(cAreaLst);
	# cFerroLst = similar(cAreaLst);
	
	# for iBeta = 1 : length(betaLst), iRatio = 1 : length(cRatioLst), iFerroRatio = 1 : length(cFerroRatioLst), iSgnArea = 1 : length(sgnAreaLst), iSgnPerim = 1 : length(sgnPerimLst), iInit0 = 1 : length(isInit0Lst)
		# cRatio = cRatioLst[iRatio];
		# beta = betaLst[iBeta];
		# cFerroRatio = cFerroRatioLst[iFerroRatio];
		# sgnArea = sgnAreaLst[iSgnArea];
		# sgnPerim = sgnPerimLst[iSgnPerim];
		# cArea, cPerim, cFerro = paramsFromBetaRatio( beta, cRatio, cFerroRatio, sgnArea, sgnPerim );
		
		# cAreaLst[iBeta, iRatio, iFerroRatio, iSgnArea, iSgnPerim] = cArea;
		# cPerimLst[iBeta, iRatio, iFerroRatio, iSgnArea, iSgnPerim] = cPerim;
		# cFerroLst[iBeta, iRatio, iFerroRatio, iSgnArea, iSgnPerim] = cFerro;
	# end
	
	# return cAreaLst, cPerimLst, cFerroLst;
# end

function paramsFromBetaRatio( beta, cRatio, cFerroRatio, sgnArea, sgnPerim )
	cRatioSq = sqrt(cRatio);
	cArea = sgnArea * beta * cRatioSq;
	cPerimAbs = beta / cRatioSq;
	cPerim = sgnPerim * cPerimAbs;
	cFerro = cPerimAbs * cFerroRatio;
	
	return cArea, cPerim, cFerro; 
end

function paramsLstFromBetaRatio( betaLst, cRatioLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst )
	cAreaLst = zeros( length(betaLst), length(cRatioLst), length(cFerroRatioLst), length(sgnAreaLst), length(sgnPerimLst) );
	cPerimLst = similar(cAreaLst);
	cFerroLst = similar(cAreaLst);
	
	for iBeta = 1 : length(betaLst), iRatio = 1 : length(cRatioLst), iFerroRatio = 1 : length(cFerroRatioLst), iSgnArea = 1 : length(sgnAreaLst), iSgnPerim = 1 : length(sgnPerimLst)
		cRatio = cRatioLst[iRatio];
		beta = betaLst[iBeta];
		cFerroRatio = cFerroRatioLst[iFerroRatio];
		sgnArea = sgnAreaLst[iSgnArea];
		sgnPerim = sgnPerimLst[iSgnPerim];
		cArea, cPerim, cFerro = paramsFromBetaRatio( beta, cRatio, cFerroRatio, sgnArea, sgnPerim );
		
		cAreaLst[iBeta, iRatio, iFerroRatio, iSgnArea, iSgnPerim] = cArea;
		cPerimLst[iBeta, iRatio, iFerroRatio, iSgnArea, iSgnPerim] = cPerim;
		cFerroLst[iBeta, iRatio, iFerroRatio, iSgnArea, iSgnPerim] = cFerro;
	end
	
	return cAreaLst, cPerimLst, cFerroLst;
end

struct RadiusParamsGroup <: ParamsGroup
	cRadiusLst::Vector{Float64};
	cAngleLst::Vector{Float64};
	cFerroRatioLst::Vector{Float64};
	sgnAreaLst::Vector{Int64};
	sgnPerimLst::Vector{Int64};
	
	groupParamsLst::Vector{<:Vector};
	paramsLst::Vector{<:Array{Float64}};
	
	function RadiusParamsGroup( cRadiusLst, cAngleLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst )
		cAreaLst = zeros( length(cRadiusLst), length(cAngleLst), length(cFerroRatioLst), length(sgnAreaLst), length(sgnPerimLst) );
		cPerimLst = similar(cAreaLst);
		cFerroLst = similar(cAreaLst);
		
		for iRadius = 1 : length(cRadiusLst), iAngle = 1 : length(cAngleLst), iCFerroRatio = 1 : length(cFerroRatioLst), iSgnArea = 1 : length(sgnAreaLst), iSgnPerim = 1 : length(sgnPerimLst)
			cAreaLst[iRadius,iAngle,iCFerroRatio,iSgnArea,iSgnPerim] = sgnAreaLst[iSgnArea] * cRadiusLst[iRadius] * cos(cAngleLst[iAngle]);
			cPerimLst[iRadius,iAngle,iCFerroRatio,iSgnArea,iSgnPerim] = sgnPerimLst[iSgnPerim] * cRadiusLst[iRadius] * sin(cAngleLst[iAngle]);
			cFerroLst[iRadius,iAngle,iCFerroRatio,iSgnArea,iSgnPerim] = cFerroRatioLst[iCFerroRatio] * cRadiusLst[iRadius];
		end
		
		paramsLst = [cAreaLst, cPerimLst, cFerroLst];
		groupParamsLst = [cRadiusLst, cAngleLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst];
		
		new( cRadiusLst, cAngleLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst, groupParamsLst, paramsLst );
	end
end

function getGroupName( grpType::Type{RadiusParamsGroup} )
	return "RadiusParamsGroup";
end

function getGroupAttr( groupType::Type{RadiusParamsGroup} )
	attrLst = ["cRadius", "cAngle", "cFerroRatio", "sgnAreaLst", "sgnPerimLst"];
	return attrLst;
end

function getGroupValLst( grp::RadiusParamsGroup; digs = 3 )
	valLst = Any[];
	append!( valLst, summarizeArrAttrRound.( [grp.cRadiusLst, grp.cAngleLst, grp.cFerroRatioLst]; digs = digs ), summarizeArrBoolAttr.([grp.sgnAreaLst, grp.sgnPerimLst]) );
	return valLst;
end

struct CartesianParamsGroup <: ParamsGroup
	cAreas::Vector{Float64};
	cPerims::Vector{Float64};
	cFerros::Vector{Float64};
	
	groupParamsLst::Vector{<:Vector};
	paramsLst::Vector{<:Array{Float64}};
	
	function CartesianParamsGroup( cAreas, cPerims, cFerros )
		groupParamsLst = [cAreas,cPerims,cFerros];
		cAreaLst = zeros( length(cAreas), length(cPerims), length(cFerros) );
		cPerimLst = similar(cAreaLst);
		cFerroLst = similar(cAreaLst);
		for iA = 1 : length(cAreas), iP = 1 : length(cPerims), iF = 1 : length(cFerros)
			cAreaLst[iA,iP,iF] = cAreas[iA];
			cPerimLst[iA,iP,iF] = cPerims[iP];
			cFerroLst[iA,iP,iF] = cFerros[iF];
		end
		paramsLst = [cAreaLst, cPerimLst,cFerroLst];
		
		new( cAreas, cPerims, cFerros, groupParamsLst, paramsLst );
	end
end

function getGroupName( grpType::Type{CartesianParamsGroup} )
	return "CartesianParamsGroup";
end

function getGroupAttr( groupType::Type{CartesianParamsGroup} )
	attrLst = ["cAreas", "cPerims", "cFerros"];
	return attrLst;
end

function getGroupValLst( grp::CartesianParamsGroup; digs = 3 )
	valLst = Any[];
	append!( valLst, summarizeArrAttrRound.( [grp.cAreas, grp.cPerims, grp.cFerros]; digs = digs ) );
	return valLst;
end

function getAttrValGroupsSummarized( divNumLst::Vector{Int64}, itNumLst::Vector{Int64}, isInit0Lst::Vector{Bool}, paramsGrpLst::Vector{ParamsGroup} )
	attrLstFLst = ["div", "it", "isInit0"];
	valLstFLst = Any[];
	append!( valLstFLst, summarizeArrAttr.( [divNumLst, itNumLst] ), summarizeArrBoolAttr.([isInit0Lst]) );
	for iGrp = 1 : length(paramsGrpLst)
		append!( attrLstFLst, getGroupAttr( paramsGrpLst[iGrp] ) );
		append!( valLstFLst, getGroupValLst( paramsGrpLst[iGrp] ) );
	end
	
	return attrLstFLst, valLstFLst;
end

function getParamLstFromGroupLst( grpLst::Vector{<:ParamsGroup} )
	lnGrpLst = length(grpLst);
	cAreaLst = Vector{Array{Float64}}(undef, lnGrpLst);
	cPerimLst = similar(cAreaLst);
	cFerroLst = similar(cAreaLst);
	
	for iGrp = 1 : lnGrpLst
		cAreaLst[iGrp], cPerimLst[iGrp], cFerroLst[iGrp] = paramsLstFromGroup( grpLst[iGrp] );
	end
	
	return cAreaLst, cPerimLst, cFerroLst;
end

function runLoopMC_withParamsGroup( updaterType::(Type{T} where T <: LoopsUpdater), itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, isInit0Lst::Vector{Bool}, paramsGrpLst::Vector{ParamsGroup}; fMod = "", betaBase = 1 )
	# lnGrpLst = length(ParamsGroup);
	# cAreaLst = Vector{Array{Float64}}(undef, lnGrpLst);
	# cPerimLst = similar(cAreaLst);
	# cFerroLst = similar(cAreaLst);
	
	# for iGrp = 1 : lnGrpLst
		# cAreaLst[iParamLst], cPerimLst[iParamLst], cFerroLst[iParamLst] = paramsLstFromGroup( paramsGrpLst[iGrp] );
	# end
	cAreaLst, cPerimLst, cFerroLst = getParamLstFromGroupLst( paramsGrpLst );
	
	runLoopMC_withParamsBase( updaterType, itNumLst, divNumLst, isInit0Lst, cAreaLst, cPerimLst, cFerroLst; fMod = fMod );
end

function runLoopMC_withParamsBase( updaterType::Type{<:LoopsUpdater}, itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, isInit0Lst::Vector{Bool}, cAreaLst::Vector{Array{Float64}}, cPerimLst::Vector{Array{Float64}}, cFerroLst::Vector{Array{Float64}}; fMod = "" )
	lnParamLst = length(cAreaLst);
	@time for itNum in itNumLst, divNum in divNumLst, isInit0 in isInit0Lst, iParamLst = 1:lnParamLst
		for idParam in CartesianIndices(cAreaLst[iParamLst])
			cArea = cAreaLst[iParamLst][idParam];
			cPerim = cPerimLst[iParamLst][idParam];
			cFerro = cFerroLst[iParamLst][idParam];
			
			@time fNameSmart = loops_MC_methods( divNum, itNum; cPerim = cPerim, cArea = cArea, cFerro = cFerro, fMod = fMod, isInit0 = isInit0, updaterType = updaterType );
			GC.gc()
		end
	end
end

function genFNameLstLoopMC( updaterType::Type{<:LoopsUpdater}, itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, isInit0Lst::Vector{Bool}, paramsGrpLst::Vector{ParamsGroup}; fMain::String, fMod = "", rndDigits::Int64 = 3 )
	lnGrpLst = length( paramsGrpLst );
	cAreaLsts, cPerimLsts, cFerroLsts = getParamLstFromGroupLst( paramsGrpLst );
	fNameLst = Vector{String}(undef,0);
	# fNameArr = [ Array{String}(undef,size(cAreaLsts[iGrp])) for iIt = 1:length(itNumLst), iDiv = 1:length(divNumLst), iInit0 = 1 : length(isInit0Lst), iGrp = 1 : lnGrpLst ];
	fNameArr = [ Array{String}( undef, length(itNumLst), length(divNumLst), length(isInit0Lst), size(cAreaLsts[iGrp])... ) for iGrp = 1 : lnGrpLst ];
	for iIt = 1 : length(itNumLst), iDiv = 1 : length(divNumLst), iInit0 = 1 : length(isInit0Lst), iParamLst = 1 : length(paramsGrpLst)
		cAreaLst = cAreaLsts[iParamLst];
		cPerimLst = cPerimLsts[iParamLst];
		cFerroLst = cFerroLsts[iParamLst];
		for idParam in CartesianIndices(cAreaLst)
			itNum = itNumLst[iIt];
			divNum = divNumLst[iDiv];
			isInit0 = isInit0Lst[iInit0];
			cArea = cAreaLst[idParam];
			cPerim = cPerimLst[idParam];
			cFerro = cFerroLst[idParam];
			
			attrLst, valLst = getAttrValLstLoopsMC( divNum, itNum, cArea, cPerim; cFerro = cFerro );
			fModFull = getFModLoopsMC( fMod, updaterType, isInit0 );
			fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModFull );
			
			push!( fNameLst, fName );
			# fNameArr[iIt,iDiv,iInit0,iParamLst][idParam] = fName;
			fNameArr[iParamLst][iIt,iDiv,iInit0,idParam] = fName;
			GC.gc()
		end
	end
	
	fLstMain = fMain * "_fLst";
	fLstJld2Main = fLstMain * "Jld2";
	fLstJld2MasterMain = fLstMain * "Jld2" * "Master";
	attrLstFLst, valLstFLst = getAttrValGroupsSummarized( divNumLst, itNumLst, isInit0Lst, paramsGrpLst );
	fModFLst = getFModLoopsMC( fMod, updaterType );
	fLstName = fNameFunc( fLstMain, attrLstFLst, valLstFLst, ".txt"; fMod = fModFLst );
	fLstJld2MasterName = fNameFunc( fLstJld2MasterMain, attrLstFLst, valLstFLst, jld2Type; fMod = fModFLst );
	fLstJld2NameLst = Vector{String}(undef, length(paramsGrpLst));
	for iGrp = 1 : length(paramsGrpLst)
		attrLstFLst, valLstFLst = getAttrValGroupsSummarized( divNumLst, itNumLst, isInit0Lst, paramsGrpLst[[iGrp]] );
		fLstJld2NameLst[iGrp] = fNameFunc( fLstJld2Main, attrLstFLst, valLstFLst, jld2Type; fMod = fModFLst );
		save( fLstJld2NameLst[iGrp], "fNameArr", fNameArr[iGrp] );
	end
	writedlm( fLstName, fNameLst );
	save( fLstJld2MasterName, "fLstNameLst", fLstJld2NameLst, "fNameArrGroups", fNameArr );
	
	open(dirLog * fNameFileLstLst, "w") do io
		println(io, fLstName);
	end
	open( dirLog * fNameFileLstJld2Lst, "w" ) do io
		println(io, fLstJld2MasterName);
	end
	
	return fLstName;
end

function saveParamsLoopMC( updaterType::(Type{T} where T <: LoopsUpdater), itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, isInit0Lst::Vector{Bool}, paramsGrpLst::Vector{ParamsGroup}; fMod = "", rndDigits = 3 )
	cAreaLst, cPerimLst, cFerroLst = getParamLstFromGroupLst( paramsGrpLst );
	cAreaLstStr = (a->string.(a)).(cAreaLst);
	cPerimLstStr = (a->string.(a)).(cPerimLst);
	cFerroLstStr = (a->string.(a)).(cFerroLst);
	
	fParamMain = "loopsMC_params";
	attrLstParamLst, valLstParamLst = getAttrValGroupsSummarized( divNumLst, itNumLst, isInit0Lst, paramsGrpLst );
	fModFull = getFModLoopsMC( fMod, updaterType );
	
	grpNameLst = getGroupName.( paramsGrpLst );
	grpParamsNameLst = getGroupAttr.( paramsGrpLst );
	grpParamsLst = getGroupParamsLst.( paramsGrpLst );
	
	fParamName = fNameFunc( fParamMain, attrLstParamLst, valLstParamLst, jld2Type; fMod = fModFull );
	save( fParamName, "itNumLst", itNumLst, "divNumLst", divNumLst, "isInit0Lst", isInit0Lst, "grpNameLst", grpNameLst, "grpParamsNameLst", grpParamsNameLst, "grpParamsLst", grpParamsLst, "cAreaLst", cAreaLst, "cAreaLstStr", cAreaLstStr, "cPerimLst", cPerimLst, "cPerimLstStr", cPerimLstStr, "cFerroLst", cFerroLst, "cFerroLstStr", cFerroLstStr );
	
	open( dirLog * fNameSaveParamsLst, "w" ) do io
		println( io, fParamName );
	end
	
	return fParamName;
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
		# cRatioSq = sqrt(cRatio);
		# cArea = sgnArea * beta * cRatioSq;
		# cPerimAbs = beta / cRatioSq;
		# cPerim = sgnPerim * cPerimAbs;
		# cFerro = cPerimAbs * cFerroRatio;
		cArea, cPerim, cFerro = paramsFromBetaRatio( beta, cRatio, cFerroRatio, sgnArea, sgnPerim );
		
		# attrLst = cFerro == 0 ? attrLstLoops : attrLstLoopsFerro;
		# valLst = Any[divNum, itNum, cArea, cPerim, betaBase];
		# if cFerro != 0
			# push!( valLst, cFerro );
		# end
		attrLst, valLst = getAttrValLstLoopsMC( divNum, itNum, cArea, cPerim; beta = betaBase, cFerro = cFerro );
		# fModFull = fModWithMethod;
		# if isInit0
			# fModFull = Utils.strAppendWith_( fModWithMethod, "isInit0" );
		# end
		fModFull = getFModLoopsMC( fMod, updaterType, isInit0 );
		fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModFull );
		
		push!( fNameLst, fName );
		fNameArr[iIt,iDiv,iBeta,iRatio,iFerroRatio,iSgnArea,iSgnPerim,iInit0] = fName;
		GC.gc()
	end
	
	fLstMain = fMain * "_fLst";
	fLstJld2Main = fLstMain * "Jld2";
	attrLstFLst = ["div", "it", "beta", "cRatio", "cFerroRatio", "isInit0", "sgnArea", "sgnPerim"];
	valIsInit0 = length(isInit0Lst) > 1 ? 2 : isInit0Lst[1];
	valSgnArea = length(sgnAreaLst) > 1 ? 2 : sgnAreaLst[1];
	valSgnPerim = length(sgnPerimLst) > 1 ? 2 : sgnPerimLst[1];
	# valLstFLst = Any[divNumLst[[1,end]],itNumLst[[1,end]],betaLst[[1,end]],round.(cRatioLst[[1,end]]; digits=rndDigits),cFerroRatioLst[[1,end]], isInit0Lst[[1,end]], sgnAreaLst[[1,end]], sgnPerimLst[[1,end]]];
	valLstFLst = Any[divNumLst[[1,end]],itNumLst[[1,end]],betaLst[[1,end]],round.(cRatioLst[[1,end]]; digits=rndDigits),cFerroRatioLst[[1,end]], valIsInit0, valSgnArea, valSgnPerim];
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

function genFNameLstLoopMC_old( updaterType::(Type{T} where T <: LoopsUpdater), itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, betaLst::Vector{Float64}, cRatioLst::Vector{Float64}, cFerroRatioLst::Vector, sgnAreaLst::Vector{Int64}, sgnPerimLst::Vector{Int64}, isInit0Lst::Vector{Bool}; fMain::String, fMod = "", betaBase = 1, isFModMethod = true, rndDigits = 3 )
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
		# cRatioSq = sqrt(cRatio);
		# cArea = sgnArea * beta * cRatioSq;
		# cPerimAbs = beta / cRatioSq;
		# cPerim = sgnPerim * cPerimAbs;
		# cFerro = cPerimAbs * cFerroRatio;
		cArea, cPerim, cFerro = paramsFromBetaRatio( beta, cRatio, cFerroRatio, sgnArea, sgnPerim );
		
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
	
	fLstMain = fMain * "_fLst";
	fLstJld2Main = fLstMain * "Jld2";
	attrLstFLst = ["div", "it", "beta", "cRatio", "cFerroRatio", "isInit0", "sgnArea", "sgnPerim"];
	valIsInit0 = length(isInit0Lst) > 1 ? 2 : isInit0Lst[1];
	valSgnArea = length(sgnAreaLst) > 1 ? 2 : sgnAreaLst[1];
	valSgnPerim = length(sgnPerimLst) > 1 ? 2 : sgnPerimLst[1];
	# valLstFLst = Any[divNumLst[[1,end]],itNumLst[[1,end]],betaLst[[1,end]],round.(cRatioLst[[1,end]]; digits=rndDigits),cFerroRatioLst[[1,end]], isInit0Lst[[1,end]], sgnAreaLst[[1,end]], sgnPerimLst[[1,end]]];
	valLstFLst = Any[divNumLst[[1,end]],itNumLst[[1,end]],betaLst[[1,end]],round.(cRatioLst[[1,end]]; digits=rndDigits),cFerroRatioLst[[1,end]], valIsInit0, valSgnArea, valSgnPerim];
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

function runLoopMC_withParams( updaterType::(Type{T} where T <: LoopsUpdater), itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, betaLst::Vector{Float64}, cRatioLst::Vector{Float64}, cFerroRatioLst::Vector, sgnAreaLst::Vector{Int64}, sgnPerimLst::Vector{Int64}, isInit0Lst::Vector{Bool}; fMod = "", betaBase = 1 )
	@time for itNum in itNumLst, divNum in divNumLst, beta in betaLst, cRatio in cRatioLst, cFerroRatio in cFerroRatioLst, sgnArea in sgnAreaLst, sgnPerim in sgnPerimLst, isInit0 in isInit0Lst
		# cRatioSq = sqrt(cRatio);
		# cArea = sgnArea * beta * cRatioSq;
		# cPerimAbs = beta / cRatioSq;
		# cPerim = sgnPerim * cPerimAbs;
		# cFerro = cPerimAbs * cFerroRatio;
		cArea, cPerim, cFerro = paramsFromBetaRatio( beta, cRatio, cFerroRatio, sgnArea, sgnPerim );
		# @time fNameSmart = methodLoops( divNum, itNum; cPerim = cPerim, cArea = cArea, cFerro = cFerro, beta = betaBase, fMod = fModSmart, isInit0 = isInit0 );
		@time fNameSmart = loops_MC_methods( divNum, itNum; cPerim = cPerim, cArea = cArea, cFerro = cFerro, beta = betaBase, fMod = fMod, isInit0 = isInit0, updaterType = updaterType );
		GC.gc()
	end
end

function genFNameLstLoopMC_old( updaterType::(Type{T} where T <: LoopsUpdater), itNumLst::Vector{Int64}, divNumLst::Vector{Int64}, betaLst::Vector{Float64}, cRatioLst::Vector{Float64}, cFerroRatioLst::Vector, sgnAreaLst::Vector{Int64}, sgnPerimLst::Vector{Int64}, isInit0Lst::Vector{Bool}; fMain::String, fMod = "", betaBase = 1, isFModMethod = true, rndDigits = 3 )
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
		# cRatioSq = sqrt(cRatio);
		# cArea = sgnArea * beta * cRatioSq;
		# cPerimAbs = beta / cRatioSq;
		# cPerim = sgnPerim * cPerimAbs;
		# cFerro = cPerimAbs * cFerroRatio;
		cArea, cPerim, cFerro = paramsFromBetaRatio( beta, cRatio, cFerroRatio, sgnArea, sgnPerim );
		
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
	
	fLstMain = fMain * "_fLst";
	fLstJld2Main = fLstMain * "Jld2";
	attrLstFLst = ["div", "it", "beta", "cRatio", "cFerroRatio", "isInit0", "sgnArea", "sgnPerim"];
	valIsInit0 = length(isInit0Lst) > 1 ? 2 : isInit0Lst[1];
	valSgnArea = length(sgnAreaLst) > 1 ? 2 : sgnAreaLst[1];
	valSgnPerim = length(sgnPerimLst) > 1 ? 2 : sgnPerimLst[1];
	# valLstFLst = Any[divNumLst[[1,end]],itNumLst[[1,end]],betaLst[[1,end]],round.(cRatioLst[[1,end]]; digits=rndDigits),cFerroRatioLst[[1,end]], isInit0Lst[[1,end]], sgnAreaLst[[1,end]], sgnPerimLst[[1,end]]];
	valLstFLst = Any[divNumLst[[1,end]],itNumLst[[1,end]],betaLst[[1,end]],round.(cRatioLst[[1,end]]; digits=rndDigits),cFerroRatioLst[[1,end]], valIsInit0, valSgnArea, valSgnPerim];
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
	cAreaLst = zeros( length(betaLst), length(cRatioLst), length(cFerroRatioLst), length(sgnAreaLst), length(sgnPerimLst) );
	cPerimLst = similar(cAreaLst);
	cFerroLst = similar(cAreaLst);
	for iIt = 1 : length(itNumLst), iDiv = 1 : length(divNumLst), iBeta = 1 : length(betaLst), iRatio = 1 : length(cRatioLst), iFerroRatio = 1 : length(cFerroRatioLst), iSgnArea = 1 : length(sgnAreaLst), iSgnPerim = 1 : length(sgnPerimLst), iInit0 = 1 : length(isInit0Lst)
		itNum = itNumLst[iIt];
		cRatio = cRatioLst[iRatio];
		beta = betaLst[iBeta];
		cFerroRatio = cFerroRatioLst[iFerroRatio];
		sgnArea = sgnAreaLst[iSgnArea];
		sgnPerim = sgnPerimLst[iSgnPerim];
		# cRatioSq = sqrt(cRatio);
		# cArea = sgnArea * beta * cRatioSq;
		# cPerimAbs = beta / cRatioSq;
		# cPerim = sgnPerim * cPerimAbs;
		# cFerro = cPerimAbs * cFerroRatio;
		cArea, cPerim, cFerro = paramsFromBetaRatio( beta, cRatio, cFerroRatio, sgnArea, sgnPerim );
		
		cAreaLst[iBeta, iRatio, iFerroRatio, iSgnArea, iSgnPerim] = cArea;
		cPerimLst[iBeta, iRatio, iFerroRatio, iSgnArea, iSgnPerim] = cPerim;
		cFerroLst[iBeta, iRatio, iFerroRatio, iSgnArea, iSgnPerim] = cFerro;
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
