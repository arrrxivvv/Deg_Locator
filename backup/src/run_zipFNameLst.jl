using Utils
using DegLocatorDiv

function printFLst( mSzLst, itNumLst, seedLst, fMainLst, fModLst, fMethodModLst; nDim = 3, divLstBase = [80, 16, 16] )
	zipFLstName = "zipFNameLst.txt";
	zipNameIO = open( zipFLstName, "w" );
	resBase = 16;
	resTop = 2*resBase;
	nBase = 10; 
	nTop = 50;

	# divLstBase = [80, 32, 32];
	divLst = copy(divLstBase);
	
	isSqrt = true;
	isFirst = true;

	for mSz in mSzLst, itNum in itNumLst, seed in seedLst, fMain in fMainLst, fMod in fModLst, fMethodMod in fMethodModLst
		if isSqrt
			divLst .= Int64.( floor.( sqrt(mSz / nBase) .* divLstBase ) );
		else
			res = Int64( floor( (resTop - resBase) / (nTop - nBase) * (mSz - nBase) ) ) + resBase;
			divLst[2] = res;
			divLst[3] = res;
		end
		attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim = nDim );
		fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = [fMethodMod, fMod] );
		
		if isFirst
			isFirst = false;
		else
			write( zipNameIO, "\n" );
		end
		write( zipNameIO, fName );
	end

	close( zipNameIO );
end

fMainLst = ["deg"];
itNumLst = [100];
seedLst = [1000];
nDim = 2;
fExt = jld2Type;
divLstBase = [64,64];

fModLst = [""];
fMethodModLst = ["flux"];

mSzLst = [10:5:60;];

printFLst( mSzLst, itNumLst, seedLst, fMainLst, fModLst, fMethodModLst; nDim = nDim, divLstBase = divLstBase );