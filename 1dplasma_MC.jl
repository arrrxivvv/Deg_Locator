using Utils
using LinearAlgebra
using JLD
using Logging; using LoggingExtras;
using Random
using NPZ

function updateNarr( nArrPol, Uthis, T )
	param_divide_num = length(nArrPol[1]);
	nArrPolNext = deepcopy(nArrPol);
	for idPol = 1:2
		nArrPolNext[idPol][rand(1:param_divide_num)] += 1;
		idRm = findall(!iszero,nArrPolNext[idPol])[rand(1:end)];
		nArrPolNext[idPol][idRm] -= 1;
	end
	
	Unext = UofNarr( nArrPolNext );
	
	if Unext <= Uthis
		nArrPol = nArrPolNext;
	elseif ( ( rnd = rand() ) <= exp( - (Unext - Uthis) / T ) )
		# @info( rnd, ", ", exp( - (Unext - Uthis) / T ) );
		nArrPol = nArrPolNext;
	end
	
	return nArrPol, Unext;
end

function EofNarr( nArrPol )
	nArr = nArrPol[1] - nArrPol[2];
	param_divide_num = length(nArr);
	Earr = zeros( param_divide_num );
	idHalf = Int64( param_divide_num / 2 );
	EthisBase = zeros( param_divide_num );
	EthisBase[1:idHalf] .= 1;
	EthisBase[idHalf+1:end] .= -1;
	for loc = 1 : param_divide_num
		Earr += nArr[loc] * circshift( EthisBase, loc-1 );
	end
	
	return Earr;
end

function UofNarr( nArrPol )
	return sum( EofNarr(nArrPol).^2 );
end

function recordArr!( nArrAllPol, Uall, nArrPol, Uthis, idTr )
	for idPol in 1:2
		nArrAllPol[idPol][idTr,:] = nArrPol[idPol];
	end
	Uall[idTr] = Uthis;
end

function onedPlasmaMC( param_divide_num, trialNum, Np, T )
	nPArr = zeros( Int64, param_divide_num );
	nnArr = zeros( Int64, param_divide_num );
	nArrPol = [nPArr, nnArr];
	
	nPArrAll = zeros( Int64, (trialNum, param_divide_num) );
	nnArrAll = zeros( Int64, (trialNum, param_divide_num) );
	nArrAllPol = [nPArrAll, nnArrAll];
	Uall = zeros( trialNum );
	
	for n in 1:Np, idPol = 1:2
		nArrPol[idPol][rand(1:param_divide_num)] += 1;
	end
	Uthis = UofNarr( nArrPol );
	
	idTr = 1;
	recordArr!(nArrAllPol, Uall, nArrPol, Uthis, idTr);
	
	for idTr = 2 : trialNum
		nArrPol, Uthis = updateNarr( nArrPol, Uthis, T );
		recordArr!(nArrAllPol, Uall, nArrPol, Uthis, idTr);
	end
	
	fileNameMain = "MC";
	fileNameAttr = string( "param_divide_num", "_", param_divide_num, "trialNum", "_", trialNum, "Np", "_", Np, "T", "_", T );
	fileNameNarr = string( fileNameMain, "_Narr", "_", fileNameAttr, npyType );
	fileNameU = string( fileNameMain, "_U", "_", fileNameAttr, npyType );
	
	fileNameJLD = string( fileNameMain, "_", fileNameAttr, jldType );
	
	save( fileNameJLD, "nArrAllPol", nArrAllPol, "Uall", Uall );
	
	nArrAllArr = arrLstToArr( nArrAllPol );
	
	npzwrite( fileNameNarr, nArrAllArr );
	npzwrite( fileNameU, Uall );
	
	return nArrAllPol, Uall;
end

function onedPlasmaClstFromFile( param_divide_num, trialNum, Np, T )
	fileNameMain = "MC";
	fileNameAttr = string( "param_divide_num", "_", param_divide_num, "trialNum", "_", trialNum, "Np", "_", Np, "T", "_", T );	
	fileNameJLD = string( fileNameMain, "_", fileNameAttr, jldType );
	
	data = load( fileNameJLD );
	Uall = data["Uall"];
	nArrAllPol = data["nArrAllPol"];
	
	nArrAll = nArrAllPol[1] - nArrAllPol[2];
	
	cLst = zeros( Int64, param_divide_num );
	cLstSqAv = zeros( param_divide_num );
	for idTr = 1 : trialNum
		nArr = nArrAll[idTr, :];
		for loc = 1 : param_divide_num-1
			cLst[loc+1] = cLst[loc] + nArr[loc];
		end
		cLstSqAv += cLst.^2;
	end
	cLstSqAv = cLstSqAv / trialNum;
	
	fileNameClstNpy = string( fileNameMain, "_Clst_", fileNameAttr, npyType );
	npzwrite( fileNameClstNpy, cLstSqAv );
end

function onedPlasmaFfunc( param_divide_num, trialNum, Np, T )
	fileNameMain = "MC";
	fileNameAttr = string( "param_divide_num", "_", param_divide_num, "trialNum", "_", trialNum, "Np", "_", Np, "T", "_", T );	
	fileNameJLD = string( fileNameMain, "_", fileNameAttr, jldType );
	
	data = load( fileNameJLD );
	Uall = data["Uall"];
	nArrAllPol = data["nArrAllPol"];
	
	fFuncPol = [ zeros( param_divide_num, param_divide_num ) for idPol = 1 : 4 ];
	
	idPol = 1;
	for pol1 = 1:2, pol2 = 1:2
		for idTr = 1 : trialNum
			nArrP = nArrAllPol[pol1][idTr,:];
			nArrN = nArrAllPol[pol2][idTr,:];
			
			fFuncThis = nArrP .* nArrN';
			if pol1 == pol2
				fFuncThis -= sqrt.(Diagonal(fFuncThis));
			end
			
			fFuncPol[idPol] += fFuncThis
			
		end
		idPol += 1;
	end
	fFuncPol /= ( trialNum * Np^2 );
	
	fileNameFfuncNpy = string( fileNameMain, "_Ffunc_",  fileNameAttr, npyType );
	npzwrite( fileNameFfuncNpy, arrLstToArr( fFuncPol ) );
end

function onedPlasmaProfile( param_divide_num, trialNum, Np, Tlst )
	for T in Tlst
		onedPlasmaMC( param_divide_num, trialNum, Np, T );
		onedPlasmaClstFromFile( param_divide_num, trialNum, Np, T );
		onedPlasmaFfunc( param_divide_num, trialNum, Np, T );
	end
end

