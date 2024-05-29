#Loops_MC module

function SwitchingUpdater{T_tuple}( params ) where {T_tuple}
	updaterLst = ntuple( ii -> T_tuple.parameters[ii](params), length(T_tuple.parameters) );
	
	SwitchingUpdater( updaterLst );
end

function getUpdaterFMod( updaterType::Type{<:SwitchingUpdater} )
	return join( getUpdaterFMod.( (updaterType.parameters[1]).parameters ), "_" );
end

function updateLoops( updater::SwitchingUpdater, flipChecker::SwitchingFlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	if getSwitchLstLen( updater ) != getSwitchLstLen( flipChecker )
		error( "Switching updater and flipchecker list length not matched" );
	end
	updaterNow = updater.updaterTup[updater.counterRef[]];
	flipCheckerNow = flipChecker.flipCheckerTup[updater.counterRef[]];
	
	updateLoops( updaterNow, flipCheckerNow, BfieldLst, linkLst, linkFerroLst, params );
	
	updater.counterRef[] = mod( updater.counterRef[], getSwitchLstLen(updater) );
	updater.counterRef[] += 1;
end

function getSwitchLstLen( updater::SwitchingUpdater )
	return length(updater.updaterTup);
end

function updateLoops( updater::LoopsUpdater, flipChecker::FlipChecker, flipProposer::SwitchingFlipProposer, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	flipProposerNow = flipProposer.flipProposerTup[flipProposer.counterRef[]];
	
	updateLoops( updater, flipChecker, flipProposerNow, BfieldLst, linkLst, linkFerroLst, params );
	
	flipProposer.counterRef[] = mod( flipProposer.counterRef[], getSwitchLstLen(flipProposer) );
	flipProposer.counterRef[] += 1;
end





function getUpdaterFMod( updaterType::Type{SingleUpdater} )
	return "upSingle";
end

function updateLoops( updater::SingleUpdater, flipChecker::FlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	pos = rand(params.posLst);
	dim = rand(1:params.nDimB);
	
	# if flipCheck( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst )
		# flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = dim );
	# end
	flipCheckDoIt( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
end

function updateLoops( updater::SingleUpdater, flipChecker::FlipChecker, flipProposer::FlipProposer, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	pos = rand(params.posLst);
	dim = rand(1:params.nDimB);
	
	flipCheckDoIt( flipChecker, flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
end


function getUpdaterFMod( updaterType::Type{ABUpdater} )
	return "upAB";
end

function updateLoops( updater::ABUpdater, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	for dim = 1 : params.nDim
		for iAB = 1 : 2
			Threads.@threads for pos in updater.posABLst[iAB,dim]
				# @time begin
				iArea = BfieldLst[dim][pos] + 1;
				iL = 1;
				for iLnkDim in 1 : params.nDimLayer
					dimLink = params.linkDimLst[dim][iLnkDim];
					dimLinkSh = params.linkDimShLst[dim][iLnkDim];
					iL += linkLst[dimLink][pos];
					iL += linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
				end
				iLFerro = 1;
				for iLnkDim = 1 : params.nDimLayer
					dimLink = params.linkDimLst[dim][iLnkDim];
					dimLinkSh = params.linkDimShLst[dim][iLnkDim];
					iLFerro += linkFerroLst[iLnkDim,dim][pos];
					iLFerro += linkFerroLst[iLnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]];
				end
				# end
				# @time begin
				pSwitchRand = rand();
				if pSwitchRand < updater.pFlipLst[iLFerro, iL, iArea];
					flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = dim );
				end
				# end
				# @infiltrate
			end
		end
	end
end



function getUpdaterFMod( updaterType::Type{AB2dUpdater} )
	return "upAB2d";
end

function updateLoops( updater::AB2dUpdater, flipChecker::FlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	if params.nDim != 2
		error( "nDim != 2" );
	end
	
	dim = 1;
	for iAB = 1 : 2
		for pos in updater.posABLst[iAB]
			# if flipCheck( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst )
				# flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = dim );
			# end
			flipCheckDoIt( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
		end
	end
end



function updateLoops( updater::StaggeredCubeUpdater, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	for iAdv = 1 : params.nDim+1
		rand!( updater.randIDimLst, updater.iDimLst );
		rand!( updater.randIShLst, updater.iIsShLst );
		Threads.@threads for idStag in updater.idStagLst
			posCube = updater.posStagCubeLst[iAdv][idStag];
			pos = updater.posShOrNotLst[updater.randIDimLst[idStag]][updater.randIShLst[idStag]][posCube];
			iArea = BfieldLst[updater.randIDimLst[idStag]][pos] + 1;
			iL = 1;
			for iLnkDim in 1 : params.nDimLayer
				dimLink = params.linkDimLst[updater.randIDimLst[idStag]][iLnkDim];
				dimLinkSh = params.linkDimShLst[updater.randIDimLst[idStag]][iLnkDim];
				iL += linkLst[dimLink][pos];
				iL += linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
			end
			iLFerro = 1;
			for iLnkDim = 1 : params.nDimLayer
				dimLink = params.linkDimLst[updater.randIDimLst[idStag]][iLnkDim];
				dimLinkSh = params.linkDimShLst[updater.randIDimLst[idStag]][iLnkDim];
				iLFerro += linkFerroLst[iLnkDim,updater.randIDimLst[idStag]][pos];
				iLFerro += linkFerroLst[iLnkDim,updater.randIDimLst[idStag]][params.posLstShLst[dimLinkSh,1][pos]];
			end
			pSwitch = rand();
			pFlip = updater.pFlipLst[iLFerro, iL, iArea];
			if pSwitch < pFlip
				flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = updater.randIDimLst[idStag] );
			end
			# @infiltrate
		end
		# @infiltrate
	end
end

function getUpdaterFMod( updaterType::Type{StaggeredCubeUpdater} )
	return "upStagCube";
end




function getUpdaterFMod( updaterType::Type{StaggeredCubeUpdaterBase} )
	return "upStagCube";
end

function updateLoops( updater::StaggeredCubeUpdaterBase, flipChecker::FlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	for iAdv = 1 : params.nDim+1
		rand!( updater.randIDimLst, updater.iDimLst );
		rand!( updater.randIShLst, updater.iIsShLst );
		Threads.@threads for idStag in updater.idStagLst
			posCube = updater.posStagCubeLst[iAdv][idStag];
			pos = updater.posShOrNotLst[updater.randIDimLst[idStag]][updater.randIShLst[idStag]][posCube];
			
			# if flipCheck( flipChecker, params, updater.randIDimLst[idStag], pos, BfieldLst, linkLst, linkFerroLst )
				# flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = updater.randIDimLst[idStag] );
			# end
			flipCheckDoIt( flipChecker, params, updater.randIDimLst[idStag], pos, BfieldLst, linkLst, linkFerroLst );
		end
	end
end

function updateLoops( updater::StaggeredCubeUpdaterBase, flipChecker::FlipChecker, flipProposer::FlipProposer, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	for iAdv = 1 : params.nDim+1
		rand!( updater.randIDimLst, updater.iDimLst );
		rand!( updater.randIShLst, updater.iIsShLst );
		Threads.@threads for idStag in updater.idStagLst
			posCube = updater.posStagCubeLst[iAdv][idStag];
			pos = updater.posShOrNotLst[updater.randIDimLst[idStag]][updater.randIShLst[idStag]][posCube];
			
			flipCheckDoIt( flipChecker, flipProposer, params, updater.randIDimLst[idStag], pos, BfieldLst, linkLst, linkFerroLst );
		end
	end
end



function getUpdaterFMod( updaterType::Type{CubeUpdater} )
	return "cubeOnly";
end

function updateLoops( updater::CubeUpdater, flipChecker::CubeFlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	posCube = rand( params.posLst );
	dimDummy = 1;
	
	flipCheckDoIt( flipChecker, params, dimDummy, posCube, BfieldLst, linkLst, linkFerroLst );
	# if flipCheck( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst )
		# flipBfieldCubeAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos );
	# end
end




function getUpdaterFMod( updaterType::Type{CubeStaggeredCubeUpdater} )
	return "cubeOnlyStaggeredCube";
end

function updateLoops( updater::CubeStaggeredCubeUpdater, flipChecker::CubeFlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	dimDummy = 1;
	for iAdv = 1 : params.nDim+1
		Threads.@threads for idStag in updater.idStagLst
			posCube = updater.posStagCubeLst[iAdv][idStag];
			
			flipCheckDoIt( flipChecker, params, dimDummy, posCube, BfieldLst, linkLst, linkFerroLst );
			# if flipCheck( flipChecker, params, dimDummy, posCube, BfieldLst, linkLst, linkFerroLst )
				# flipBfieldCubeAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = posCube );
			# end
		end
	end
end


