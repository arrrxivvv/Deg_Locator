using Utils
using DegLocatorDiv
using JLD
using NPZ

function getShLst( fFunc, dimSh, sh )
	dim = length( size( fFunc ) );
	shLst = zeros( Int64, dim );
	shLst[dimSh] = sh;
	return shLst;
end

function shDim( fFunc, dimSh, sh )
	circshift( fFunc, getShLst( fFunc, dimSh, sh ) );
end

function reverseFfunc( fFunc, dimSh )
	fFuncFlip = reverse( fFunc; dims=dimSh );
	shDim( fFunc, dimSh, 1 );
	return fFuncFlip;
end

function piShFfunc( fFunc, dimSh )
	ln = size( fFunc, dimSh );
	halfDivide = Int64( ln / 2 );
	return shDim( fFunc, dimSh, halfDivide );
end

function piShReverseFfunc( fFunc, dimSh )
	fFuncReverse = reverseFfunc(fFunc, dimSh);
	return piShFunc( fFuncreverse, dimsh );
end

function flipFfuncSimul( fFunc )
	fFuncReverse = reverseFfunc( fFunc, 1 );
	fFuncReverse = reverseFfunc( fFuncReverse, 2 );
	
	fFuncPiSh = piShFfunc( fFunc, 1 );
	fFuncPiSh = piShFfunc( fFuncPiSh, 2 );
	
	fFuncPiShFlip = piShFfunc( fFuncReverse, 1 );
	fFuncPiShFlip = piShFfunc( fFuncPiShFlip, 2 );
	
	flipLst = [fFunc, fFuncReverse, fFuncPiSh, fFuncPiShFlip];
	
	return flipLst, length(flipLst);
end

function flipFfuncIndv( fFunc )
	flipLst = [fFunc];
	numDim = length( size(fFunc) );
	for dim = 1 : numDim
		flipLstTmp = [];
		for fFuncTmp in flipLst
			fFuncReverse = reverseFfunc( fFuncTmp, dim );
			fFuncPiSh = piShFfunc( fFuncTmp, dim );
			fFuncPiShFlip = piShFfunc( fFuncReverse, dim );
			append!( flipLstTmp, [fFuncReverse, fFuncPiSh, fFuncPiShFlip] );
		end
		append!( flipLst, flipLstTmp );
	end
	
	return flipLst, length(flipLst);
end

function flipFfuncNone( fFunc )
	return [fFunc], 1;
end

# function flipFfunc( fFunc, param_divide_num )
	# fFuncFlip = reverse( fFunc; dims=1 );
	# fFuncFlip = reverse( fFuncFlip; dims=2 );
	# fFuncFlip = circshift( fFuncFlip, (1,1) );
	
	# halfDivide = Int64( param_divide_num );
	# fFuncFlip2 = circshift( fFunc, (halfDivide, halfDivide) );
	# fFuncFlip3 = circshift( fFuncFlip, (halfDivide, halfDivide) );
	
	# return [fFuncFlip, fFuncFlip2, fFuncFlip3];
# end

function locs12ToFFunc!( locs1, locs2, fFuncTemp )
	fFuncTemp .= 0;
	locs1z = locs1[:,1];
	locs2z = locs2[:,1];
	for ix in locs1z, iy in locs2z
		fFuncTemp[ix,iy] += 1;
	end
end

function locs12ToFDiffAbs!( locs1, locs2, fFuncTemp )
	fFuncTemp .= 0;
	ln = size( fFuncTemp, 1 );
	for i1 in 1:size(locs1,1), i2 in 1:size(locs2,1)
		ix = locs1[i1,:];
		iy = locs2[i2,:];
		iDiff = Int64.( abs.( ix .- iy ) ) .+ 1;
		# iDiff = mod.( Int64.( ( ix .- iy ) ) , ln ) .+ 1;
		fFuncTemp[ Tuple(iDiff)... ] += 1;
	end
end

function locs12ToFDiffMod!( locs1, locs2, fFuncTemp )
	fFuncTemp .= 0;
	ln = size( fFuncTemp, 1 );
	for i1 in 1:size(locs1,1), i2 in 1:size(locs2,1)
		ix = locs1[i1,:];
		iy = locs2[i2,:];
		# iDiff = Int64.( abs.( ix .- iy ) ) .+ 1;
		iDiff = mod.( Int64.( ( ix .- iy ) ) , ln ) .+ 1;
		fFuncTemp[ Tuple(iDiff)... ] += 1;
	end
end

function locs12ToFDiffFull!( locs1, locs2, fFuncTemp )
	fFuncTemp .= 0;
	ln = size( fFuncTemp, 1 );
	param_divide_num = Int64( (ln+1)/2 );
	for i1 in 1:size(locs1,1), i2 in 1:size(locs2,1)
		ix = locs1[i1,:];
		iy = locs2[i2,:];
		iDiff = Int64.( ( ix .- iy ) ) .+ param_divide_num;
		fFuncTemp[ Tuple(iDiff)... ] += 1;
	end
end

function locs12ToFDiff!( locs1, locs2, fFuncTemp )
	fFuncTemp .= 0;
	ln = size( fFuncTemp, 1 );
	for i1 in 1:size(locs1,1), i2 in 1:size(locs2,1)
		ix = locs1[i1,:];
		iy = locs2[i2,:];
		# iDiff = Int64.( abs.( ix .- iy ) ) .+ 1;
		iDiff = mod.( Int64.( ( ix .- iy ) ) , ln ) .+ 1;
		# println( typeof(iDiff), ' ', size(iDiff), ' ', iDiff );
		fFuncTemp[ Tuple(iDiff)... ] += 1;
	end
end

function jointStatFromLocs( param_divide_num, locs1, locs2 )
	fFunc1 = zeros( param_divide_num );
	fFunc2 = zeros( param_divide_num );
	for ix in locs1
		fFunc1[ix] += 1;
	end
	for iy in locs2
		fFunc2[iy] += 1;
	end
	fFunc2d = fFunc1 .* fFunc2';
	
	return fFunc2d;
end

function fFunc_stat_fine( posLocLst, negLocLst, param_divide_num, N, num_it )
	locLst = [ posLocLst, negLocLst ];
	
	idPosNegNum = 4;
	fFuncLst = zeros( N, idPosNegNum, param_divide_num, param_divide_num );	
	
	for it = 1 : num_it
		for n = 1:N
			idPosNeg = 1;
			for idPol1 = 1:2
				locs1 = locLst[idPol1][it][n][:,1];
				for idPol2 = 1:2
					locs2 = locLst[idPol2][it][n][:,1];
					fFuncTmp = joinStatFromLocs( param_divide_num, locs1, locs2 );
					if n>1
						locs1Prev = locLst[idPol1][it][n-1][:,1];
						locs2Prev = locLst[idPol2][it][n-1][:,1];
						fFunc12 = jointStatFromLocs( param_divide_num, locs1Prev, locs2Prev );
						fFunc21 = jointStatFromLocs( param_divide_num, locs2Prev, locs2Prev );
					end
					fFuncLst[ n, idPosNeg, :, : ] += fFuncTemp;
					idPosNeg += 1;
				end
			end
		end
	end
	
	return fFuncLst;
end

function fFunc_stat_base( fFuncFunc!, flipFfunc, posLocLst, negLocLst, fFuncLn, N, num_it; idPolEnd = 2 )
	locLst = [ posLocLst, negLocLst ];
	
	idPosNegNum = 4;
	fFuncLst = zeros( N, idPosNegNum, fFuncLn, fFuncLn );
	fFuncTemp = zeros( fFuncLn, fFuncLn );
	
	for n = 1:N
		for it = 1 : num_it
			idPosNeg = 1;	
			for idPol1 = 1:idPolEnd
				locs1 = locLst[idPol1][it][n];
				for idPol2 = 1:idPolEnd
					locs2 = locLst[idPol2][it][n];
					fFuncTemp .= 0;
					fFuncFunc!( locs1, locs2, fFuncTemp );
					# fFuncLst[ n, idPosNeg, :, : ] += fFuncTemp;
					fFuncFlipLst, symFact = flipFfunc( fFuncTemp );
					# for fFuncFlip in fFuncFlipLst
						# fFuncLst[ n, idPosNeg, :, : ] += fFuncFlip;
					# end
					fFuncLst[ n, idPosNeg, :, : ] += sum(fFuncFlipLst) ./ symFact;
					idPosNeg += 1;
				end
			end
		end
		# symFact = 4;
		fFuncLst[ n, :, :, : ] = fFuncLst[ n, :, :, : ] ./ num_it; #  ./ symFact;
	end
	
	return fFuncLst;
end

function fFunc_stat_base_nd( fFuncFunc!, flipFfunc, posLocLst, negLocLst, fFuncLnLst, N, num_it; idPolEnd = 2 )
	locLst = [ posLocLst, negLocLst ];
	
	idPosNegNum = 4;
	fFuncLst = zeros( N, idPosNegNum, fFuncLnLst... );
	fFuncTemp = zeros( fFuncLnLst... );
	
	for n = 1:N
		for it = 1 : num_it
			idPosNeg = 1;	
			for idPol1 = 1:idPolEnd
				locs1 = locLst[idPol1][it][n];
				for idPol2 = 1:idPolEnd
					locs2 = locLst[idPol2][it][n];
					fFuncTemp .= 0;
					fFuncFunc!( locs1, locs2, fFuncTemp );
					# fFuncLst[ n, idPosNeg, :, : ] += fFuncTemp;
					fFuncFlipLst, symFact = flipFfunc( fFuncTemp );
					# for fFuncFlip in fFuncFlipLst
						# fFuncLst[ n, idPosNeg, :, : ] += fFuncFlip;
					# end
					fFuncLst[ n, idPosNeg, .. ] += sum(fFuncFlipLst) ./ symFact;
					idPosNeg += 1;
				end
			end
		end
		# symFact = 4;
		fFuncLst[ n, .. ] = fFuncLst[ n, .. ] ./ num_it; #  ./ symFact;
	end
	
	return fFuncLst;
end

function fFunc_stat( posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = 2 )
	fFuncLst = fFunc_stat_base( locs12ToFFunc!, flipFfuncSimul, posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = idPolEnd );
	return fFuncLst;
end

function fFunc_diff_stat( posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = 2 )
	fFuncLst = fFunc_stat_base( locs12ToFDiff!, flipFfuncNone, posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = idPolEnd );
	return fFuncLst;
end

# function fFunc_diff_stat( posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = 2 )
	# fFuncLst = fFunc_stat_base( locs12ToFDiffFull!, flipFfuncNone, posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = idPolEnd );
	# return fFuncLst;
# end


function fFunc_diff_full_stat( posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = 2 )
	fFuncLn = 2*param_divide_num - 1;
	fFuncLst = fFunc_stat_base( locs12ToFDiffFull!, flipFfuncNone, posLocLst, negLocLst, fFuncLn, N, num_it; idPolEnd = idPolEnd );
	return fFuncLst;
end

function fFunc_diff_full_stat_3d( posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = 2 )
	fFuncLn = 2*param_divide_num - 1;
	param_dim = 3;
	fFuncLnLst = fFuncLn .* ones( Int64, param_dim );
	fFuncLst = fFunc_stat_base_nd( locs12ToFDiffFull!, flipFfuncNone, posLocLst, negLocLst, fFuncLnLst, N, num_it; idPolEnd = idPolEnd );
	return fFuncLst;
end

function fFunc_stat_from_file_base( fFunc_statFunc, oFileNameMain, N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = "" )
	fileNameAttr = fileNameAttrFunc( N, param_divide, num_it, seedFedStart );
	if fileNameMod != ""
		fileNameAttr = string( fileNameAttr, "_", fileNameMod );
	end
	fileNameMain = "deg";
	if param_dim == 2 || param_dim == 3
		fileNameAttr = string( "dim_", param_dim, "_", fileNameAttr );
	end
	oFileName = string( fileNameMain, "_", fileNameAttr );
	varFileName = string( oFileName, jldType );
	@time data = load(varFileName);
	
	param_divide = Int64.(param_divide);
	if param_divide isa Array{Int64}
		param_divide_num = param_divide[1];
	else
		param_divide_num = param_divide;
	end
	
	idPolEnd = 2;
	if param_dim == 2
		idPolEnd = 1;
	end
	
	@time fFuncLst = fFunc_statFunc( data["posLocLst"], data["negLocLst"], param_divide_num, N, num_it; idPolEnd = idPolEnd );
	
	# oFileNameMain = "ffunc";
	if fileOutNameMod != ""
		fileNameAttr = string( fileNameAttr, "_", fileOutNameMod );
	end
	fileNameNpy = string( oFileNameMain, "_",  fileNameAttr, npyType );
	npzwrite( fileNameNpy, fFuncLst );
	
	return fFuncLst;
end

function fFunc_stat_from_file( N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = "" )
	oFileName = "ffunc";
	fFuncLst = fFunc_stat_from_file_base( fFunc_stat, oFileName, N, param_divide, num_it, seedFedStart, fileNameMod, param_dim, fileOutNameMod );
end

function fFunc_stat_diff_from_file( N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = "" )
	oFileName = "fdiff";
	fFuncLst = fFunc_stat_from_file_base( fFunc_diff_stat, oFileName, N, param_divide, num_it, seedFedStart, fileNameMod, param_dim, fileOutNameMod );
end

function fFunc_stat_diff_full_from_file( N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = "" )
	oFileName = "fdiffFull";
	fFuncLst = fFunc_stat_from_file_base( fFunc_diff_full_stat, oFileName, N, param_divide, num_it, seedFedStart, fileNameMod, param_dim, fileOutNameMod );
end

function fFunc_stat_diff_full_3d_from_file( N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = "" )
	oFileName = "fdiff3dFull";
	fFuncLst = fFunc_stat_from_file_base( fFunc_diff_full_stat_3d, oFileName, N, param_divide, num_it, seedFedStart, fileNameMod, param_dim, fileOutNameMod );
end
