using DegLocatorDiv
using Utils

nlst = [10:5:60;];

divLstBase = [80, 64, 64];
divLst = copy( divLstBase );

resBase = 16;
nBase = 10;

itNum = 100;

seed = 1000;
nDim = 3;

fMod = "flux";

for n in nlst
	divLst .= Int64.( floor.( sqrt( n / nBase ) .* divLstBase ) );
	
	println( "n = $n, res = $(divLst[end])" );
	with_logger( errLogger )do 
		# @time divB_profile_flux( n, divLst, itNum, seed; enumSaveMem = memEig );
		# @time distillLocsFromFile( itNum, n, divLst, seed; dim = nDim, isNewFlux = true, fMod = fMod );
		@time deltaN_avg_lst_fromFile( [n], itNum, divLst, seed; dim = nDim, fileNameMod = "_"*fMod, opt = "distill2" );
		@time deltaN_var_lst_fromFile( [n], itNum, divLst, seed; dim = nDim, fileNameMod = "_"*fMod, opt = "distill2" );
	end
end
