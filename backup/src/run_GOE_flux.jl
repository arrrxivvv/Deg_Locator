using DegLocatorDiv
using Utils

nlst = [10:5:60;];

nDim = 3;

divLstBase = [64,64];
divLstBase = [128,128,16];
divLst = copy( divLstBase );

resBase = 16;
nBase = 10;

itNum = 10;

seed = 1000;

fMod = "";

itNumStop = 1;

for n in nlst
	divLst .= Int64.( floor.( sqrt( n / nBase ) .* divLstBase ) );
	
	println( "n = $n, res = $(divLst)" );
	with_logger( errLogger )do 
		# @time divB_profile_flux( n, divLst, itNum, seed; enumSaveMem = memEig, nDim = nDim );
		# @time divB_profile_GOE_layered( n, divLst, itNum, 1000; fMod = "" )
		
		# @time zakArr_corr_GOE_from_file( n, divLst, itNum, seed; fMod = fMod, dim = nDim, itNumStop = itNumStop )
		@time deg_GOE3_zak_resave( n, divLst, itNum, seed; fMod = fMod );
	end
end
