using Utils
using LinearAlgebra
using DegLocatorDiv
using JLD
using NPZ

struct FfObj
	N::Int64;
	param_divide::Vector{Int64};	
	fFuncLn::Int64;
	
	idPolEnd::Int64;
	
	itNum::Int64;
	
	locLst::Vector{Vector{Vector{Array{Int64}}}};
	
	fFuncLst::Array{Float64};
	
	tmp::Vector{Array{Int64}};
	flip::Vector{Array{Int64}};
	flip2::Vector{Array{Int64}};
	avg::Vector{Array{Float64}};
	
	function FfObj( N, param_dim, param_divide, itNum, posLocLst, negLocLst )
		param_divide = Int64.(param_divide);
		numThread = Threads.nthreads();
		
		if param_divide isa Array{Int64}
			fFuncLn = param_divide[1];
		else
			fFuncLn = param_divide;
		end
		
		idPolEnd = 2;
		if param_dim == 2
			idPolEnd = 1;
		end
		
		locLst = [posLocLst, negLocLst];
		
		fFuncLst = zeros( N, idPolEnd^2, fFuncLn, fFuncLn );	
		
		tmp = [zeros( Int64, fFuncLn, fFuncLn ) for it = 1 : numThread];
		flip = [zeros( Int64, fFuncLn, fFuncLn ) for it = 1 : numThread];
		flip2 = [zeros( Int64, fFuncLn, fFuncLn ) for it = 1 : numThread];
		avg = [zeros( Float64, fFuncLn, fFuncLn ) for it = 1 : numThread];
		
		new( N, param_divide, fFuncLn, idPolEnd, itNum, locLst, fFuncLst, tmp, flip, flip2, avg );
	end
end

function getShLst( fFunc, dimSh, sh )
	dim = length( size( fFunc ) );
	shLst = zeros( Int64, dim );
	shLst[dimSh] = sh;
	return shLst;
end

function shDim!( fFunc, fFuncSh, dimSh, sh )
	circshift!( fFuncSh, fFunc, getShLst( fFunc, dimSh, sh ) );
end

function reverseFfunc!( fFunc, flip, dimSh )
	shDim!( fFunc, flip, dimSh, 1 );
	myreverse!( flip; dims=dimSh );
end

function piShFfunc!( fFunc, fFuncSh, dimSh )
	ln = size( fFunc, dimSh );
	halfDivide = Int64( ln / 2 );
	shDim!( fFunc, fFuncSh, dimSh, halfDivide );
end

function piShReverseFfunc!( fFunc, dimSh )
	fFuncReverse = reverseFfunc(fFunc, dimSh);
	return piShFunc( fFuncreverse, dimsh );
end

function flipFfuncSimul!( fObj, dummyAvg )
	fObj.avg[thrId()] .= fObj.tmp[thrId()]; 
	fObj.flip[thrId()] .= fObj.tmp[thrId()];
	
	symLn = 4;
	
	reverseFfunc!( fObj.tmp[thrId()], fObj.flip2[thrId()], 1 );
	reverseFfunc!( fObj.flip2[thrId()], fObj.flip[thrId()], 2 );
	
	fObj.avg[thrId()] .+= fObj.flip[thrId()];
	
	piShFfunc!( fObj.flip[thrId()], fObj.flip2[thrId()], 1 );
	piShFfunc!( fObj.flip2[thrId()], fObj.flip[thrId()], 2 );
	
	fObj.avg[thrId()] .+= fObj.flip[thrId()];
	
	piShFfunc!( fObj.tmp[thrId()], fObj.flip2[thrId()], 1 );
	piShFfunc!( fObj.flip2[thrId()], fObj.flip[thrId()], 2 );
	
	fObj.avg[thrId()] .+= fObj.flip[thrId()];
	
	fObj.avg[thrId()] ./= symLn;
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

function flipFfuncNone!( fObj, dummyAvg )
	# return [fFunc], 1;
end

# function flipFfunc( fFunc, param_divide_num )
	# flip = reverse( fFunc; dims=1 );
	# flip = reverse( flip; dims=2 );
	# flip = circshift( flip, (1,1) );
	
	# halfDivide = Int64( param_divide_num );
	# flip2 = circshift( fFunc, (halfDivide, halfDivide) );
	# flip3 = circshift( flip, (halfDivide, halfDivide) );
	
	# return [flip, flip2, flip3];
# end

function locs12ToFFunc!( locs1, locs2, tmp )
	tmp .= 0;
	locs1z = view( locs1,:,1 );
	locs2z = view( locs2,:,1 );
	for ix in locs1z, iy in locs2z
		tmp[ix,iy] += 1;
	end
end

function locs12ToFDiffAbs!( locs1, locs2, tmp )
	tmp .= 0;
	ln = size( tmp, 1 );
	for i1 in 1:size(locs1,1), i2 in 1:size(locs2,1)
		ix = locs1[i1,:];
		iy = locs2[i2,:];
		iDiff = Int64.( abs.( ix .- iy ) ) .+ 1;
		# iDiff = mod.( Int64.( ( ix .- iy ) ) , ln ) .+ 1;
		tmp[ Tuple(iDiff)... ] += 1;
	end
end

function locs12ToNorm!( locs1, locs2, param_divide )
	tmp .= 0;
	ln = size( tmp, 1 );
	for i1 in 1:size(locs1,1), i2 in 1:size(locs2,1)
		ix = locs1[i1,:];
		iy = locs2[i2,:];
		iNorm = norm( ( ix .- iy ) ./ param_divide );
		# iDiff = mod.( Int64.( ( ix .- iy ) ) , ln ) .+ 1;
		tmp[ Tuple(iDiff)... ] += 1;
	end
end

function locs12ToFDiffMod!( locs1, locs2, tmp )
	tmp .= 0;
	ln = size( tmp, 1 );
	for i1 in 1:size(locs1,1), i2 in 1:size(locs2,1)
		ix = locs1[i1,:];
		iy = locs2[i2,:];
		# iDiff = Int64.( abs.( ix .- iy ) ) .+ 1;
		iDiff = mod.( Int64.( ( ix .- iy ) ) , ln ) .+ 1;
		tmp[ Tuple(iDiff)... ] += 1;
	end
end

function locs12ToFDiffFull!( locs1, locs2, tmp )
	tmp .= 0;
	ln = size( tmp, 1 );
	param_divide_num = Int64( (ln+1)/2 );
	for i1 in 1:size(locs1,1), i2 in 1:size(locs2,1)
		ix = locs1[i1,:];
		iy = locs2[i2,:];
		iDiff = Int64.( ( ix .- iy ) ) .+ param_divide_num;
		tmp[ Tuple(iDiff)... ] += 1;
	end
end

function locs12ToFDiff!( locs1, locs2, tmp )
	tmp .= 0;
	ln = size( tmp, 1 );
	for i1 in 1:size(locs1,1), i2 in 1:size(locs2,1)
		ix = locs1[i1,:];
		iy = locs2[i2,:];
		# iDiff = Int64.( abs.( ix .- iy ) ) .+ 1;
		iDiff = mod.( Int64.( ( ix .- iy ) ) , ln ) .+ 1;
		# println( typeof(iDiff), ' ', size(iDiff), ' ', iDiff );
		tmp[ Tuple(iDiff)... ] += 1;
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
					fFuncLst[ n, idPosNeg, :, : ] += tmp;
					idPosNeg += 1;
				end
			end
		end
	end
	
	return fFuncLst;
end

function fFunc_stat_norm!( fObj )
	Threads.@threads for n = 1:fObj.N
		for it = 1 : fObj.itNum
			idPosNeg = 1;	
			for idPol1 = 1:fObj.idPolEnd
				locs1 = fObj.locLst[idPol1][it][n];
				for idPol2 = 1:fObj.idPolEnd
					locs2 = fObj.locLst[idPol2][it][n];
					fObj.tmp[thrId()] .= 0;
					# @info( "fFuncFunc" )
					fFuncFunc!( locs1, locs2, fObj.tmp[thrId()] );
					# fFuncLst[ n, idPosNeg, :, : ] += tmp;
					# @info( "flip" )
					fObj.avg[thrId()] .= fObj.tmp[thrId()];
					flipFfunc!( fObj, fObj.avg[thrId()] );
					# @info( "add to total" )
					view( fObj.fFuncLst, n, idPosNeg, :, : ) .+= fObj.avg[thrId()];
					idPosNeg += 1;
					# readline();
				end
			end
		end
		# symFact = 4;
		# fFuncLst[ n, :, :, : ] = fFuncLst[ n, :, :, : ] ./ num_it; #  ./ symFact;
	end
	
	fObj.fFuncLst ./= fObj.itNum;
end

function fFunc_stat_base!( fObj, fFuncFunc!, flipFfunc! )
	Threads.@threads for n = 1:fObj.N
		for it = 1 : fObj.itNum
			idPosNeg = 1;	
			for idPol1 = 1:fObj.idPolEnd
				locs1 = fObj.locLst[idPol1][it][n];
				for idPol2 = 1:fObj.idPolEnd
					locs2 = fObj.locLst[idPol2][it][n];
					fObj.tmp[thrId()] .= 0;
					# @info( "fFuncFunc" )
					fFuncFunc!( locs1, locs2, fObj.tmp[thrId()] );
					# fFuncLst[ n, idPosNeg, :, : ] += tmp;
					# @info( "flip" )
					fObj.avg[thrId()] .= fObj.tmp[thrId()];
					flipFfunc!( fObj, fObj.avg[thrId()] );
					# @info( "add to total" )
					view( fObj.fFuncLst, n, idPosNeg, :, : ) .+= fObj.avg[thrId()];
					idPosNeg += 1;
					# readline();
				end
			end
		end
		# symFact = 4;
		# fFuncLst[ n, :, :, : ] = fFuncLst[ n, :, :, : ] ./ num_it; #  ./ symFact;
	end
	
	fObj.fFuncLst ./= fObj.itNum;
end

function fFunc_stat_base_nd( fFuncFunc!, flipFfunc, posLocLst, negLocLst, fFuncLnLst, N, num_it; idPolEnd = 2 )
	locLst = [ posLocLst, negLocLst ];
	
	idPosNegNum = idPolEnd^2;
	fFuncLst = zeros( N, idPosNegNum, fFuncLnLst... );
	tmp = zeros( fFuncLnLst... );
	
	for n = 1:N
		for it = 1 : num_it
			idPosNeg = 1;	
			for idPol1 = 1:idPolEnd
				locs1 = locLst[idPol1][it][n];
				for idPol2 = 1:idPolEnd
					locs2 = locLst[idPol2][it][n];
					tmp .= 0;
					fFuncFunc!( locs1, locs2, tmp );
					# fFuncLst[ n, idPosNeg, :, : ] += tmp;
					flipLst, symFact = flipFfunc( tmp );
					# for flip in flipLst
						# fFuncLst[ n, idPosNeg, :, : ] += flip;
					# end
					fFuncLst[ n, idPosNeg, .. ] += sum(flipLst) ./ symFact;
					idPosNeg += 1;
				end
			end
		end
		# symFact = 4;
		fFuncLst[ n, .. ] = fFuncLst[ n, .. ] ./ num_it; #  ./ symFact;
	end
	
	return fFuncLst;
end

function fFunc_stat( fObj )
	fFuncLst = fFunc_stat_base!( fObj, locs12ToFFunc!, flipFfuncSimul! );
	# flipFfuncSimul!
	# flipFfuncNone!
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

function fFunc_stat_from_file_base( fFunc_statFunc!, oFileNameMain, N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = "" )
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
	@info("load data");
	@time data = load(varFileName);
	
	fFuncSolver = FfObj( N, param_dim, param_divide, num_it,data["posLocLst"], data["negLocLst"] );
	
	@info("construct corr func");
	@time fFunc_statFunc!( fFuncSolver );
	
	# oFileNameMain = "ffunc";
	if fileOutNameMod != ""
		fileNameAttr = string( fileNameAttr, "_", fileOutNameMod );
	end
	fileNameNpy = string( oFileNameMain, "_",  fileNameAttr, npyType );
	npzwrite( fileNameNpy, fFuncSolver.fFuncLst );
	
	# return fFuncLst;
end

function fFuncStatNormFromFile( N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = "" )
	oFileName = "ffunc";
	fFuncLst = fFunc_stat_from_file_base( fFunc_stat, oFileName, N, param_divide, num_it, seedFedStart, fileNameMod, param_dim, fileOutNameMod );
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
