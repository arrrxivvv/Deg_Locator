const attrLstBase = ["N", "param_divide", "instanceNum", "seed"];
const attrLstDim = vcat(["dim"],attrLstBase);

const varNamePosLoc = "posLocLst";
const varNameNegLoc = "negLocLst";
const varNamePosN = "posNlst";
const varNameNegN = "negNlst";
const varNameLocDensity = "locDensity";
const varNameLocCorrAvg = "locCorrLstAvg";
const varNameLocCorrStd = "locCorrLstStd";
const varNameParityArr = "parity2dArr";
const varNameZacAvg = "zacAvgArr";
const varNameZacCorr = "parityCorr";
const varNameZacCorrAvg = "zacCorrAvg";
const varNameZacCorrStd = "zacCorrStd";

const fNameDegGOE = "deg_GOE_3d_full";
const fNameLocDensityGOE = "locDenstiy_GOE";
const fNameLocCorr = "locCorr";
const fNameZacArr = "parityArr";
const fNameZacAvg = "zacAvgArr";
const fNameZacCorr = "parityCorr_GOE";
const fNameZacCorrAvg = "parityCorrAvg";
const fNameZacCorrStd = "parityCorrStd";

function fAttrValFunc( N, param_divide, itNum, seed, dim; alpha = 0, scale = 1 )
	attrLst = deepcopy(attrLstBase);
	valLst = [N, param_divide, itNum, seed];
	optAttrLst = ["scale","alpha"];
	optValLst = [scale,alpha];
	optDefLst = [1,0];
	if !isnothing(dim)
		attrLst = insert!( attrLst, 1, "dim" );
		valLst = insert!( valLst, 1, dim );
	end
	for io = 1 : length(optValLst)
		if optValLst[io] != optDefLst[io]
			push!( attrLst, optAttrLst[io] );
			push!( valLst, optValLst[io] );
		end
	end
	return attrLst, valLst;
end

function fileNameAttrFunc( N, param_divide, num_it, seedFedStart; dim = nothing )
	fileNameAttr = string( "N_", N, "_param_divide_", param_divide, "_instanceNum_", num_it, "_seed_", seedFedStart );
	if !isnothing(dim)
		fileNameAttr = string( "dim_", dim, "_", fileNameAttr );
	end
	return fileNameAttr;
end

function fNameAttrLstFunc( N, param_divide, itNum, seed; dim = nothing, alpha = nothing )
	attrLst = deepcopy( attrLstBase );
	valLst = [N, param_divide, itNum, seed];
	if !isnothing(dim)
		insert!( attrLst, 1, "dim" );
		insert!( valLst, 1, dim );
	end
	if !isnothing(alpha)
		append!( attrLst, ["alpha"] );
		append!( valLst, alpha );
	end
	return attrLst, valLst;
end

function fAttrOptLstFunc( N, param_divide, itNum, seed; dim = nothing, scale = degOptDefaultLst[1], ratio = degOptDefaultLst[2], alpha = degOptDefaultLst[3], enumSaveMem = memNone, thresNM = rtFndDefaultLst[1], thresEDeg = rtFndDefaultLst[2] )
	attrLst, valLst = fNameAttrLstFunc( N, param_divide, itNum, seed; dim = dim );
	optValLst = [scale, ratio, alpha];
	for ii = eachindex(degOptAttrLst)
		if optValLst[ii] != degOptDefaultLst[ii]
			push!( attrLst, deOptAttrLst[ii] );
			push!( valLst, optValLst[ii] );
		end
	end
	rtFndValLst = [thresNM, thresEDeg];
	if enumSaveMem == rootFind
		for ii = eachindex(rtFndAttrLst)
			push!( attrLst, rtFndAttrLst[ii] );
			push!( valLst, rtFndValLst[ii] );
		end
	end
	return attrLst, valLst;
end

function fNameFunc( fNameMain, attrLst, valLst, fExt; fMod = "" )
	fName = fNameMain;
	for ii = 1 : length(attrLst)
		fName = fName * "_" * attrLst[ii] * "_" * string( valLst[ii] );
	end
	if fMod != ""
		fName = fName * "_" * fMod;
	end
	fName = fName * fExt;
	return fName;
end

# function fileNameAttrFunc( N, param_divide, num_it, seedFedStart )
	# fileNameAttr = string( "N_", N, "_param_divide_", param_divide, "_instanceNum_", num_it, "_seed_", seedFedStart );
	# return fileNameAttr;
# end
