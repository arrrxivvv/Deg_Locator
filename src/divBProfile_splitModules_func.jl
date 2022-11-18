using Utils
using Statistics
using JLD2
using Logging; using LoggingExtras;

function divB_profile_new( mSz, divLst, itNum, seedFed; nDim = 3 )
	if seedFed > 0
		Random.seed!(seedFed);
	end
	minNum = 0;
	maxNum = 2*pi;
	paramsFull = degParamsInit( mSz, divLst, minNum, maxNum, nDim );
	matsFull = matsGridHThreaded( paramsFull, threaded_zeros(ComplexF64,mSz,mSz) );
	degBerrysFull = degBerrysInit( paramsFull, matsFull; isFullInit = true );
	
	totalNumLst = zeros(Int64,itNum);
	
	for it = 1 : itNum
		print( "\r Iteration: " * string(it) * " / " * string(itNum) );
		Hlst = DegLocatorDiv.HlstFunc(H_GUE,nDim,mSz);
		HmatFun = (H,xLst) -> Hmat_3comb_ratio!( H, xLst, Hlst );
		startNextEigen( matsFull );
		
		@info("Eigen:")
		Utils.@timeInfo eigenAll( matsFull; HmatFun = HmatFun );

		@info("Link:")
		Utils.@timeInfo DegLocatorDiv.linksCalcAll( degBerrysFull );
		
		@info("Bfield:")
		Utils.@timeInfo DegLocatorDiv.BfieldCalcAll( degBerrysFull );
		
		@info("DivB:")
		Utils.@timeInfo DegLocatorDiv.divBCalcAll( degBerrysFull );
		
		totalNumLst[it] = sum(sum( (x->abs.(x).>1e-9).(degBerrysFull.divBLst) ));
	end
	return totalNumLst;
end
