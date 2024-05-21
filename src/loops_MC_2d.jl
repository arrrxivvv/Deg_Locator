# Loops_MC module

function loops_MC_2d( divNum = 64, itNum = 10000; cArea = 0, cPerim = 1 )
	initializer = ConstantInitializer();
	updaterType = SingleUpdater;
	flipChecker = NeighborFlipChecker( cArea, cPerim );
	
	fName = loops_MC_method_Base( divNum, itNum; updaterType, flipChecker, initializer, nDim = 2 );
	
	
end
