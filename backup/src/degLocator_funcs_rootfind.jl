using LinearAlgebra
using Arpack
using NelderMeadTest

function locate_rootFind( degObj, HmatFun )
	locCount = zeros(Int64, degObj.N-1);
	locLst = Vector{Matrix{Float64}}(undef,degObj.N-1);
	# locLst = [ [] for iN = 1 : degObj.N-1 ];
	cntNM = 0;
	for iN = 1 : degObj.N-1
		tFull = @timed begin
		# print("\e[s");
		Threads.@threads for idCart in degObj.posLst
		# for idCart in degObj.posLst
			# Threads.@threads for iSimp = 1 : degObj.lnSimpAll
			for iSimp = 1 : degObj.lnSimpAll
				# print("\r", idCart, ", iSimp = ", iSimp, ", ")
				for iVal = 1 : degObj.param_dim+1
					# @view( degObj.ptLst[:,iVal] ) .= degObj.param_mesh[degObj.shSimpCartLst[iVal,iSimp][idCart]];
					@view( degObj.ptLst[:,iVal] ) .= degObj.shMeshLst[iVal,iSimp][idCart];
					broadcastAssign!( degObj.eLstTmp, degObj.Elst[degObj.shSimpCartLst[iVal,iSimp][idCart]] );
					degObj.valLst[iVal] = degObj.eLstTmp[iN+1] - degObj.eLstTmp[iN];
				end
				# tFull = @timed 
				cntNM = nmOpt!( (xLst -> eValDiff!( getThrInst( degObj.Htmp ), getThrInst( degObj.eLstTmp ), HmatFun, xLst, iN, degObj )), degObj.param_dim, getThrInst(degObj.ptLst), getThrInst(degObj.valLst), getThrInst(degObj.ixLst), getThrInst(degObj.ptLstEnd), getThrInst(degObj.ptLstTmp), getThrInst(degObj.valLstTmp), degObj.thresNM );
				# @info( timeMemStr( tFull.time, tFull.bytes ) )
				@view( degObj.locLstFromSimp[idCart][:,iSimp] ) .= wrapCoorArr!( @view( degObj.ptLstEnd[:,1] ), degObj.param_max );
				degObj.eDiffLst[idCart,iSimp] = degObj.valLst[1];
				# @infiltrate cond = idCart == CartesianIndex(22,1,1)
			end 
		end
		end
		@info( timeMemStr( tFull.time, tFull.bytes ) )
		# @infiltrate
		locCount[iN], locLst[iN] = distill_rootFind_n( degObj );
	end
	posCount = locCount;
	negCount = copy(locCount);
	posLocLst = locLst;
	negLocLst = deepcopy(locLst);
	return posCount, negCount, posLocLst, negLocLst;
end

function locate_rootFind_full( degObj, HmatFun )
	return locate_rootFind_base( degObj, HmatFun, nmOptDegFull! );
end

function locate_rootFind_Lanczos( degObj, HmatFun )
	return locate_rootFind_base( degObj, HmatFun, nmOptDegLanczos! );
end

function locate_rootFind_base( degObj, HmatFun, nmOptFun! )
	locCount = zeros(Int64, degObj.N-1);
	locLst = Vector{Matrix{Float64}}(undef,degObj.N-1);
	# locLst = [ [] for iN = 1 : degObj.N-1 ];
	cntNM = 0;
	for iN = 1 : degObj.N-1
		tFull = @timed begin
		# print("\e[s");
		# Threads.@threads for idCart in degObj.posLst
		for idCart in degObj.posLst
			# Threads.@threads for iSimp = 1 : degObj.lnSimpAll
			for iSimp = 1 : degObj.lnSimpAll
				# print("\r", idCart, ", iSimp = ", iSimp, ", ")
				for iVal = 1 : degObj.param_dim+1
					# @view( degObj.ptLst[:,iVal] ) .= degObj.param_mesh[degObj.shSimpCartLst[iVal,iSimp][idCart]];
					@view( degObj.ptLst[:,iVal] ) .= degObj.shMeshLst[iVal,iSimp][idCart];
					broadcastAssign!( degObj.eLstTmp, degObj.Elst[degObj.shSimpCartLst[iVal,iSimp][idCart]] );
					degObj.valLst[iVal] = degObj.eLstTmp[iN+1] - degObj.eLstTmp[iN];
				end
				tFull = @timed 				cntNM = nmOptFun!( HmatFun, iN, degObj );
				@info( timeMemStr( tFull.time, tFull.bytes ) )
				# @infiltrate
				@view( degObj.locLstFromSimp[idCart][:,iSimp] ) .= wrapCoorArr!( @view( degObj.ptLstEnd[:,1] ), degObj.param_max );
				degObj.eDiffLst[idCart,iSimp] = degObj.valLst[1];
				# @infiltrate cond = idCart == CartesianIndex(22,1,1)
			end
		end
		end
		@info( timeMemStr( tFull.time, tFull.bytes ) )
		# @infiltrate
		locCount[iN], locLst[iN] = distill_rootFind_n( degObj );
	end
	posCount = locCount;
	negCount = copy(locCount);
	posLocLst = locLst;
	negLocLst = deepcopy(locLst);
	return posCount, negCount, posLocLst, negLocLst;
end

function nmOptDegFull!( HmatFun, iN, degObj )
	fun = (xLst -> eValDiff!( getThrInst( degObj.Htmp ), getThrInst( degObj.eLstTmp ), HmatFun, xLst, iN, degObj ));
	
	return nmOptDegBase!( fun, degObj );
end

function nmOptDegLanczos!( HmatFun, iN, degObj )
	fun = ((xLst, aux, auxOut) -> degLanczos!( getThrInst(degObj.Htmp), HmatFun, xLst, aux,iN, degObj, auxOut ));
	funInit = ((xLst, auxOut) -> eValDiffAux!( getThrInst(degObj.Htmp), getThrInst( degObj.eLstTmp ), getThrInst( degObj.vLstTmp ), HmatFun, xLst, iN, degObj, auxOut ));
	
	return nmOptDegBase!( fun, degObj; yesAux = true, funInit = funInit );
end

function nmOptDegBase!( fun, degObj; funInit = nothing, yesAux = false )
	return nmOpt!( fun, degObj.param_dim, getThrInst(degObj.ptLst), getThrInst(degObj.valLst), getThrInst(degObj.ixLst), getThrInst(degObj.ptLstEnd), getThrInst(degObj.ptLstTmp), getThrInst(degObj.valLstTmp), degObj.thresNM; funInit! = funInit, yesAux = yesAux, auxLst = getThrInst( degObj.auxLst ), auxLstTmp = getThrInst( degObj.auxLstTmp ) );
end

function locate_rootFind_arrAsInd( degObj )
	for iN = 1 : degObj.N-1
		for idCart = 1 : degObj.posLst
			for iSimp = 1 : degObj.lnSimp+degObj.lnExtra
				for id = 1 : degObj.param_dim
					@view(degObj.ptIndLst[id,:]) .= idCart[id];
					degObj.ptIntLst[id,1] = wrapIntInd( degObj.ptIndLst[id,1] + degObj.shPt1[iSimp] );
				end
				for id = 1 : degObj.param_dim
					degObj.ptIntLst[id,id+1] = wrapIntInd( degObj.ptIndLst[id,id+1] + degObj.shPtOther[iSimp,id], degObj.param_divide[id] );
				end
				for iVal = 1 : degObj.param_Dim+1
					@view( degObj.ptLst[Threads.threadid()][:,iVal] ) .= degObj.param_mesh[@view(degObj.ptIntLst[:,iVal])...];
				end
			
				for iVal = 1 : degObj.param_Dim+1
					degObj.eLstTmp = degObj.Elst[@view(degObj.ptIntLst[:,iVal])...];
					degObj.valLst[iVal] = abs( degObj.eLstTmp[iN+1] - degObj.eLstTmp[iN] );
				end
				fun = (xLst -> eValDiff!( getThrInst( degObj.Htmp ), getThrInst( degObj.eLstTmp ), HmatFun, xLst, iN, degObj ));
				nmOpt!( fun, degObj.param_dim, getThrInst(degObj.ptLst), getThrInst(degObj.valLst), getThrInst(degObj.ixLst), getThrInst(degObj.ptLstEnd), getThrInst(degObj.ptLstTmp), getThrInst(degObj.valLstTmp), degObj.thresNM );
				@view( degObj.locLstFromSimp[posLocLst][:,iSimp] )  .= @view( degObj.ptLst[:,1] );
				degObj.eDiffLst[posLst,iSimp] = valLst[1];
				# @infiltrate
			end 
		end
		distill_rootFind_n( degObj, iN );
	end
end

function distill_rootFind_n( degObj )
	bordArr = zeros( Int64, degObj.param_dim, 2^degObj.param_dim );
	bordId = 1;
	shNum = 1;
	shId = 1;
	isDup = false;
	diffArr = zeros(degObj.param_dim);
	absDiff = 0;
	locLstN = [];
	locInt = zeros( Int64, degObj.param_dim );
	for ii in degObj.posLst
		for iSimp = 1 : degObj.lnSimpAll
			if degObj.eDiffLst[ii, iSimp] <= degObj.thresEDeg
				bordArr .= 0;
				loc = @view( 				degObj.locLstFromSimp[ii][:,iSimp] );
				locInt .= floor.( loc ./ degObj.param_step ) .+ 1;
				locRemain = loc .- (locInt .- 1) .* degObj.param_step;
				shNum = 1;
				for id = 1 : degObj.param_dim
					if locRemain[id] > degObj.thresDegCollision && locRemain[id] < degObj.param_step[id] - degObj.thresDegCollision
						continue
					else
						if locRemain[id] < degObj.thresDegCollision
							bndId = -1;
						elseif locRemain[id] > degObj.thresDegCollision
							bndId = 1;
						end
						for iSh = 1 : shNum
							@view( bordArr[:,shNum+iSh] ) .= @view( bordArr[:,iSh] );
							bordArr[id,shNum+iSh] = bndId;
						end
						shNum = shNum * 2;
					end
				end
				for id = 1 : degObj.param_dim
					@view( bordArr[id,:] ) .+= locInt[id];
				end
				wrapIndArr!( bordArr, degObj.param_divide );
				isDup = false;
				for iSh = 1 : shNum
					shId = getLinInd( @view(bordArr[:,iSh]), degObj.param_divide, degObj.param_dim );
					for iHash = 1 : length(degObj.locLstHash[shId])
						diffArr .= loc .- degObj.locLstHash[shId][iHash];
						absDiff = norm(diffArr);
						if absDiff <= degObj.thresDegCollision
							isDup = true;
							break;
						end
					end
				end
				# @infiltrate
				if !isDup
					locIntLin = getLinInd( locInt, degObj.param_divide, degObj.param_dim );
					push!( degObj.locLstHash[locIntLin], copy(loc) );
					push!( locLstN, copy(loc) );
				end
			end
		end
	end
	# @infiltrate
	return length(locLstN), arrLstToArr( locLstN );
end

function hash_locs( degObj, loc )
	
end

function eValDiff!( H, eLst, HmatFun, xLst, n, degObj )
	HmatFun(H, xLst);
	vLst0 = similar( H, 1, 0 );
	eigenZheevrInThr!( H, eLst, vLst0, degObj; jobz = 'N' );
	return abs( eLst[n+1] - eLst[n] );
end

function eValDiffAux!( H, eLst, vLst, HmatFun, xLst, n, degObj, auxOut )
	HmatFun(H, xLst);
	eigenZheevrInThr!( H, eLst, vLst, degObj );
	Ediff = abs( eLst[n+1] - eLst[n] );
	Eavg = ( eLst[n+1] + eLst[n] ) / 2;
	v0 = vLst[:,n];
	auxOut.Eavg = Eavg;
	auxOut.v0 .= v0;
	return Ediff; #, auxLanczos( Eavg, v0 );
end

function degLanczos!( Htmp, HmatFun, xLst, auxE_v_lst, n, degObj, auxOut )
	HmatFun(Htmp, xLst);
	
	nev = 2;
	
	eLstReal, idPosLst, idNegLst, vLst = eigsWithId!( Htmp, auxE_v_lst.v0, auxE_v_lst.Eavg );
	
	if ( isempty(idPosLst) || isempty(idNegLst) ) && nev <= size(Htmp,1)
		nev += 2;
		# with_logger(error_only_logger) do
		eLstReal, idPosLst, idNegLst, vLst = eigsWithId!( Htmp, auxE_v_lst.v0, auxE_v_lst.Eavg; nev = nev );
		# end
	end
	# @info("I got out")
	# what if this still does not work? 
	# @infiltrate
	
	# @infiltrate cond = isempty(idPosLst) || isempty(idNegLst);
	
	if isempty(idPosLst) || isempty(idNegLst)
		dE = eValDiffAux!( Htmp, getThrInst(degObj.eLstTmp), getThrInst(degObj.vLstTmp), HmatFun, xLst, n, degObj, auxOut );
	else
		ePos, idPos = findmax( id->eLstReal[id], idPosLst );
		eNeg, idNeg = findmax( id->eLstReal[id], idNegLst );
		
		dE = ePos - eNeg;
		Eavg = (ePos + eNeg) / 2;
		v0 = vLst[:,idNeg];
		auxOut.Eavg = Eavg;
		auxOut.v0 .= v0;
	end
	
	return dE; #, ePos, eNeg;
end

function eigsWithId!( Hmat, v0, Eavg; nev=2 )
	eLst, vLst = eigs( Hmat; nev=nev, v0=v0, sigma=Eavg );
	eLstReal = real.(eLst);
	
	idPosLst = findall(x->x>Eavg,eLstReal);
	idNegLst = findall(x->x<=Eavg,eLstReal);
	
	return eLstReal, idPosLst, idNegLst, vLst;
end

mutable struct auxLanczos
	Eavg::Float64;
	v0::Vector{ComplexF64};
end
