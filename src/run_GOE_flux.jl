using DegLocatorDiv
using Utils
using FilenameManip
using JLD2

dirLog = "./log/";
fNameFileLstLst = "fNameFileLstLst.txt";
fNameFileLstJld2Lst = "fNameFileLstJld2Lst.txt";
fNameSaveParamsLst = "fNameSaveParamsLst.txt";

nlst = [10:5:60;];

nDim = 3;

divLstBase = [64,64];
divLstBase = [128,128,16];
divLst = copy( divLstBase );

divLstLst = zeros(Int64, nDim, length(nlst));

resBase = 16;
nBase = 10;

itNum = 10;

seed = 1000;

fMod = "";

itNumStop = 1;

# isFileNameOnly = false;
isFileNameOnly = true;

fNameArr = Vector{String}(undef, length(nlst));
fMainBase = "deg_GOE";
fNameJldArrMain = fMainBase * "_" * "fNameArr";
fMainParam = fMainBase * "_" * "params";

attrLstJldArr = ["nLst","divBase"];
valLstJldArr = [nlst[1,end],divLstBase[1]];

fNameJldArrName = fNameFunc( fNameJldArrMain, attrLstJldArr, valLstJldArr, jld2Type );
fNameParam = fNameFunc( fMainParam, attrLstJldArr, valLstJldArr, jld2Type );

# for n in nlst
for ii = 1 : length(nlst)
	n = nlst[ii];
	divLst .= Int64.( floor.( sqrt( n / nBase ) .* divLstBase ) );
	divLstLst[:,ii] .= divLst;
	
	println( "n = $n, res = $(divLst)" );
	with_logger( errLogger )do 
		# @time divB_profile_flux( n, divLst, itNum, seed; enumSaveMem = memEig, nDim = nDim );
		# @time divB_profile_GOE_layered( n, divLst, itNum, 1000; fMod = "" )
		
		# @time zakArr_corr_GOE_from_file( n, divLst, itNum, seed; fMod = fMod, dim = nDim, itNumStop = itNumStop )
		# @time DegLocatorDiv.zakArr_corr_FFT_GOE_from_file( n, divLst, itNum, seed; fMod = fMod, dim = nDim, itNumStop = itNumStop )
		@time fNameArr[ii] = DegLocatorDiv.zakArr_corr_FFT_GOE_from_file( n, divLst, itNum, seed; fMod = fMod, dim = nDim, isFileNameOnly = isFileNameOnly );
		# @time deg_GOE3_zak_resave( n, divLst, itNum, seed; fMod = fMod );
	end
end

save( fNameJldArrName, "fNameArr", fNameArr );
save( fNameParam, "nLst", nlst, "divLstLst", divLstLst, "itNum", itNum );

open( dirLog * fNameFileLstJld2Lst, "w" ) do io
	println(io, fNameJldArrName);
end

open( dirLog * fNameSaveParamsLst, "w" ) do io
	println(io, fNameParam);
end

