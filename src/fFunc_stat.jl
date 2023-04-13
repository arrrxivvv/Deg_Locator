using Utils
using LoopVectorization
using LinearAlgebra
using DegLocatorDiv
using JLD2
using NPZ
using EllipsisNotation

struct FfObj
	N::Int64;
	param_divide::Vector{Int64};	
	fFuncLn::Int64;
	idPolEnd::Int64;
	itNum::Int64;
	
	locLst::Vector{Vector{Vector{Array{Int64}}}};
	
	fFuncLst::Array{Float64};
	fFuncLstCross::Array{Float64};
	
	tmp::Vector{Array{Int64}};
	flip::Vector{Array{Int64}};
	flip2::Vector{Array{Int64}};
	avg::Vector{Array{Float64}};
	
	tmpMsc::Vector{Vector{Array}};
	
	function FfObj( N, param_dim, param_divide, itNum, posLocLst, negLocLst, fObjOpt )
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
		
		if fObjOpt == "1d"
			fFuncLst = zeros( N, idPolEnd^2, fFuncLn, fFuncLn );
			fFuncLstCross = zeros( N-1, idPolEnd^2, fFuncLn, fFuncLn );
			
			tmp = [zeros( Int64, fFuncLn, fFuncLn ) for it = 1 : numThread];
			flip = [zeros( Int64, fFuncLn, fFuncLn ) for it = 1 : numThread];
			flip2 = [zeros( Int64, fFuncLn, fFuncLn ) for it = 1 : numThread];
			avg = [zeros( Float64, fFuncLn, fFuncLn ) for it = 1 : numThread];
			tmpMsc = [ [] for it = 1 : numThread ];
		elseif fObjOpt == "3dDiff"
			fFuncLst = zeros( N, idPolEnd^2, param_divide... );	
			tmp = [zeros( Int64, param_divide... ) for it = 1 : numThread];
			flip = [zeros( Int64, param_divide... ) for it = 1 : numThread];
			flip2 = [zeros( Int64, param_divide... ) for it = 1 : numThread];
			avg = [zeros( Float64, param_divide... ) for it = 1 : numThread];
			# tmpMsc = [ [ zeros( Int64, length(param_divide) ) ] for it = 1 : numThread ];
			tmpMsc = [ [ zeros( Int64, param_divide... ), zeros( Int64, param_divide... ) ] for it = 1 : numThread ];
		elseif fObjOpt == "3dDiffFull"
			dimsUnpacked = param_divide .*2 .- 1;
			fFuncLst = zeros( N, idPolEnd^2, dimsUnpacked... );	
			tmp = [zeros( Int64, dimsUnpacked... ) for it = 1 : numThread];
			flip = [zeros( Int64, dimsUnpacked... ) for it = 1 : numThread];
			flip2 = [zeros( Int64, dimsUnpacked... ) for it = 1 : numThread];
			avg = [zeros( Float64, dimsUnpacked... ) for it = 1 : numThread];
			tmpMsc = [ [ zeros( Int64, dimsUnpacked... ), zeros( Int64, dimsUnpacked... ) ] for it = 1 : numThread ];
		# elseif fObjOpt == "3dFull"
			# dimsUnpacked = vcat( param_divide, param_divide );
			# fFuncLst = zeros( N, idPolEnd^2, dimsUnpacked... );	
			# tmp = [zeros( Int64, dimsUnpacked... ) for it = 1 : numThread];
			# flip = [zeros( Int64, dimsUnpacked... ) for it = 1 : numThread];
			# flip2 = [zeros( Int64, dimsUnpacked... ) for it = 1 : numThread];
			# avg = [zeros( Float64, dimsUnpacked... ) for it = 1 : numThread];
			# tmpMsc = [ [ zeros( Int64, length(dimsUnpacked) ) ] for it = 1 : numThread ];
		end
		
		new( N, param_divide, fFuncLn, idPolEnd, itNum, locLst, fFuncLst, fFuncLstCross, tmp, flip, flip2, avg, tmpMsc );
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

function locs12ToFFunc!( locs1, locs2, tmp, tmpMsc )
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

function locs12ToFDiffModCirc!( locs1, locs2, tmp, tmpMsc )
	tmp .= 0;
	tmpPreSh = tmpMsc[1];
	tmpSh = tmpMsc[2];
	# param_divide = [ size( tmp )... ];
	# iDiff = tmpMsc[1];	
	# readline();
	sh = zeros( Int64, ndims(tmp) );
	for i1 in 1:size(locs1,1)
		ix = view( locs1, i1, : )
		# @time begin
		for i2 in 1:size(locs2,1)
			iy = view( locs2, i2, : );			
			tmpPreSh[ iy... ] += 1;
		end
		# end
		sh .= -( ix .- 1 );
		circshift!( tmpSh, tmpPreSh, sh );
		tmp .+= tmpSh;
		# readline()
	end
	# readline();
end

function locs12ToFDiffMod!( locs1, locs2, tmp, tmpMsc )
	tmp .= 0;
	param_divide = [ size( tmp )... ];
	# iDiff = tmpMsc[1];	
	# readline();
	iDiff = zeros( Int64, ndims(tmp) );
	for i1 in 1:size(locs1,1), i2 in 1:size(locs2,1)
		ix = view( locs1, i1, : )
		iy = view( locs2, i2, : )
		# iDiff = Int64.( abs.( ix .- iy ) ) .+ 1;
		# @time iDiff .= ix .- iy;
		# @info(ix);
		# @time iDiff .= mod.(iDiff, param_divide);
		# @info("plus");
		# @time iDiff .+= 1;
		# iDiff .= mod.( Int64.( ( ix .- iy ) ) , param_divide ) .+ 1;
		iDiff .= mod.( ( ix .- iy ) , param_divide ) .+ 1;
		# tmp[ iDiff... ] += 1;
		tmp[ iDiff[1],iDiff[2], iDiff[3] ] += 1;
		# readline();
	end
end

function locs12ToFDiffFull!( locs1, locs2, tmp, tmpMsc )
	tmp .= 0;
	# ln = size( tmp, 1 );
	param_divide = ([ size( tmp )... ].+1)/2;
	iDiff = zeros( Int64, ndims(tmp) );
	for i1 in 1:size(locs1,1), i2 in 1:size(locs2,1)
		ix = view( locs1, i1, : )
		iy = view( locs2, i2, : )
		iDiff .= Int64.( ( ix .- iy ) ) .+ param_divide;
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

function fFunc_stat_base!( fObj, fFuncFunc!, flipFfunc!; isAcross = false )
	nSh = isAcross ? 1 : 0;
	fResult = isAcross ? fObj.fFuncLstCross : fObj.fFuncLst;
	Threads.@threads for n = (1+nSh):fObj.N
	# for n = 1:fObj.N
		for it = 1 : fObj.itNum
			idPosNeg = 1;	
			for idPol1 = 1:fObj.idPolEnd
				locs1 = fObj.locLst[idPol1][it][n-nSh];
				for idPol2 = 1:fObj.idPolEnd
					locs2 = fObj.locLst[idPol2][it][n];
					fObj.tmp[thrId()] .= 0;
					fFuncFunc!( locs1, locs2, fObj.tmp[thrId()], fObj.tmpMsc[thrId()] );
					fObj.avg[thrId()] .= fObj.tmp[thrId()];
					flipFfunc!( fObj, fObj.avg[thrId()] );
					# @time fObj.tmp[thrId()] .= fObj.avg[thrId()] .+ view( fObj.fFuncLst, n, idPosNeg, .. );
					view( fResult, n-nSh, idPosNeg, .. ) .+= fObj.avg[thrId()];
					idPosNeg += 1;
				end
			end
		end
		# symFact = 4;
		# fFuncLst[ n, :, :, : ] = fFuncLst[ n, :, :, : ] ./ num_it; #  ./ symFact;
	end
	fResult ./= fObj.itNum;
end

function fFunc_stat_base_nd( fObj, FfuncFunc!, flipFfunc! )
	Threads.@threads for n = 1:N
		for it = 1 : num_it
			idPosNeg = 1;	
			for idPol1 = 1:idPolEnd
				locs1 = locLst[idPol1][it][n];
				for idPol2 = 1:idPolEnd
					locs2 = locLst[idPol2][it][n];
					fObj.tmp[thrId()] .= 0;
					fFuncFunc!( locs1, locs2, fObj.tmp[thrId()] );
					# fFuncLst[ n, idPosNeg, :, : ] += tmp;
					fObj.avg[thrId()] .= fObj.tmp[thrId()];
					flipFfunc!( fObj, fObj.avg[thrId()] );
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

function fFunc_stat_across( fObj )
	fFuncLst = fFunc_stat_base!( fObj, locs12ToFFunc!, flipFfuncSimul!, isAcross = true );
	# flipFfuncSimul!
	# flipFfuncNone!
end

function fFunc_diff_stat( fObj )
	fFuncLst = fFunc_stat_base!( fObj, locs12ToFDiff! );
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

function fFunc_diff_mod_stat_3d( fObj )
	fFunc_stat_base!( fObj, locs12ToFDiffMod!, flipFfuncNone! );
end

function fFunc_diff_full_stat_3d( fObj )
	fFunc_stat_base!( fObj, locs12ToFDiffFull!, flipFfuncNone! );
end

function fFunc_stat_from_file_base( fFunc_statFunc!, oFileNameMain, N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = "",fObjOpt = "1d"; isDistilled = false, isAcross = false )
	fileNameAttr = fileNameAttrFunc( N, param_divide, num_it, seedFedStart );
	if fileNameMod != ""
		fileNameAttr = string( fileNameAttr, "_", fileNameMod );
	end
	if !isDistilled
		fileNameMain = "deg";
	else
		fileNameMain = "locDistilled";
	end
	if param_dim == 2 || param_dim == 3
		fileNameAttr = string( "dim_", param_dim, "_", fileNameAttr );
	end
	oFileName = string( fileNameMain, "_", fileNameAttr );
	varFileName = string( oFileName, jldType );
	if isDistilled
		locLstPol = load( varFileName, "locDistilledLst" );
		oFileNameMain = string( oFileNameMain, "_distilled" );
		NforFobj = N-1;
	else
		posLocLst = load(varFileName, "posLocLst");
		negLocLst = load(varFileName, "negLocLst");
		locLstPol = [posLocLst, negLocLst];
		NforFobj = N;
	end
	@info("load data");
	# @time data = load(varFileName);
	
	fFuncSolver = FfObj( NforFobj, param_dim, param_divide, num_it,locLstPol[1], locLstPol[2], fObjOpt );
	
	@info("construct corr func");
	@time fFunc_statFunc!( fFuncSolver );
	
	# oFileNameMain = "ffunc";
	if fileOutNameMod != ""
		fileNameAttr = string( fileNameAttr, "_", fileOutNameMod );
	end
	fileNameNpy = string( oFileNameMain, "_",  fileNameAttr, npyType );
	fResult = isAcross ? fFuncSolver.fFuncLstCross : fFuncSolver.fFuncLst;
	npzwrite( fileNameNpy, fResult );
	
	return fFuncSolver.fFuncLst;
end

function fFuncStatNormFromFile( N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = "" )
	oFileName = "ffunc";
	fFuncLst = fFunc_stat_from_file_base( fFunc_stat, oFileName, N, param_divide, num_it, seedFedStart, fileNameMod, param_dim, fileOutNameMod );
end

function fFunc_stat_from_file( N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = ""; isDistilled = false )
	oFileName = "ffunc";
	fFuncLst = fFunc_stat_from_file_base( fFunc_stat, oFileName, N, param_divide, num_it, seedFedStart, fileNameMod, param_dim, fileOutNameMod, isDistilled = isDistilled );
end

function fFunc_stat_across_from_file( N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = ""; isDistilled = false )
	oFileName = "ffunc_across";
	fFuncLst = fFunc_stat_from_file_base( fFunc_stat_across, oFileName, N, param_divide, num_it, seedFedStart, fileNameMod, param_dim, fileOutNameMod, isDistilled = isDistilled, isAcross = true );
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
	fFuncLst = fFunc_stat_from_file_base( fFunc_diff_full_stat_3d, oFileName, N, param_divide, num_it, seedFedStart, fileNameMod, param_dim, fileOutNameMod, "3dDiffFull" );
	return fFuncLst;
end

function fFunc_stat_diff_mod_3d_from_file( N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = "" )
	oFileName = "fdiff3dmod";
	fFuncLst = fFunc_stat_from_file_base( fFunc_diff_mod_stat_3d, oFileName, N, param_divide, num_it, seedFedStart, fileNameMod, param_dim, fileOutNameMod, "3dDiff" );
	return fFuncLst;
end

function zak_corr_GOE_from_file( N, param_divide, itNum, seed; fMod = "", dim = nothing, fExt = jld2Type )
	divNum = param_divide[1];
	fMain = "deg_GOE_3d";
	attrLst = deepcopy(attrLstBase);
	valLst = [N, param_divide, itNum, seed];
	if !isnothing(dim)
		attrLst = insert!( attrLst, 1, "dim" );
		valLst = insert!( valLst, 1, dim );
	end
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod );
	@info(fName);
	tFull = @timed begin
		zakLstLst = load( fName, "zakLstLst" ) ;
	end
	@info("read file: ");
	@info( timeMemStr( tFull.time, tFull.bytes ) );
	parityCorrArr = [[ zeros( N ) for x = 1:divNum, y = 1:divNum ] for d = 1 : dim ];
	tmpArr1 = [[ zeros(N) for y = 1:divNum, x = 1:divNum ] for it = 1 : itNum];
	tmpArr2 = [[ zeros(N) for y = 1:divNum, x = 1:divNum ] for it = 1 : itNum];
	# @infiltrate
	itNumLess = Int64( floor( itNum / 10 ) );
	for d = 1 : dim
		println("d = ", string(d));
		for it = 1 : itNumLess
			( arr -> arr.= 0 ).(tmpArr2[it]);
			@time begin
			for x = 1 : divNum
				for y = 1 : divNum
					# print("d = ", string(d), ", it = ", string(it), ", x,y=", string([x,y]), "\r");
					@time begin
					circshift!( tmpArr1[it], parity2dLst[it][d], (-x+1,-y+1) );
					end
					# circshift!.( tmpArr2, tempArr1, -y+1 );
					@time begin
					( (vara1,varb1)->(vara1 .= vara1 .+ varb1 .* parity2dLst[it][d][x,y] ) ).( tmpArr2[it], tmpArr1[it] );
					end
					# @infiltrate
					# parityCorrArr[d] .= parityCorrArr[d] .+ ( var2->( var->var .* parity2dLst[it][d][x][y] ).(var2) ).( tmpArr ) ;
					# parityCorrArr[d] .= parityCorrArr[d] .+ ( var2->( var->var .* parity2dLst[it][d][x][y] ).(var2) ).( circshift.( circshift( parity2dLst[it][d], -x+1 ), -y+1 ) ) ;
					# parity2dLst[it][d][x][y][n] * circshift.( circshift( parity2dLst[it][d], -x+1 ), -y+1 );
				end
			end
			end
			# @infiltrate
		end
		for it = 1 : itNum
			parityCorrArr[d] .= parityCorrArr[d] .+ tmpArr2[it];
		end
	end
	# @infiltrate
	oFmain = "parityCorr_GOE";
	oFName = fNameFunc( oFmain, attrLst, valLst, jld2Type; fMod );
	save( oFName, "parityCorr", parityCorrArr );
end

# function parity_corr_GOE_from_file( N, param_divide, itNum, seed; fMod = "", dim = nothing, fExt = jld2Type )
	# divNum = param_divide[1];
	# fMain = "deg_GOE_3d_full";
	# # fMain = "parityArr";
	# attrLst = deepcopy(attrLstBase);
	# valLst = [N, param_divide, itNum, seed];
	# if !isnothing(dim)
		# attrLst = insert!( attrLst, 1, "dim" );
		# valLst = insert!( valLst, 1, dim );
	# end
	# # @infiltrate
	# fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod );
	# @info(fName);
	# tFull = @timed begin
		# parity2dLst = load( fName, "parity2dLst" ) ;
	# end
	# @info("read file: ");
	# @info( timeMemStr( tFull.time, tFull.bytes ) );
	# # @infiltrate
	# # tFull = @timed begin
		# # parity2dArr = [ parity2dLst[it][d][x,y][n] for d=1:dim, it=itNum, x=1:divNum, y=1:divNum, n=1:N ];
	# # end
	# # @info("lst to arr: ");
	# # @info( timeMemStr( tFull.time, tFull.bytes ) )
	# # @infiltrate
	# parityCorrArr = [[ zeros( N ) for x = 1:divNum, y = 1:divNum ] for d = 1 : dim ];
	# tmpArr1 = [[ zeros(N) for y = 1:divNum, x = 1:divNum ] for it = 1 : itNum];
	# tmpArr2 = [[ zeros(N) for y = 1:divNum, x = 1:divNum ] for it = 1 : itNum];
	# # @infiltrate
	# itNumLess = Int64( floor( itNum / 10 ) );
	# for d = 1 : dim
		# println("d = ", string(d));
		# for it = 1 : itNumLess
			# ( arr -> arr.= 0 ).(tmpArr2[it]);
			# @time begin
			# for x = 1 : divNum
				# for y = 1 : divNum
					# # print("d = ", string(d), ", it = ", string(it), ", x,y=", string([x,y]), "\r");
					# @time begin
					# circshift!( tmpArr1[it], parity2dLst[it][d], (-x+1,-y+1) );
					# end
					# # circshift!.( tmpArr2, tempArr1, -y+1 );
					# @time begin
					# ( (vara1,varb1)->(vara1 .= vara1 .+ varb1 .* parity2dLst[it][d][x,y] ) ).( tmpArr2[it], tmpArr1[it] );
					# end
					# # @infiltrate
					# # parityCorrArr[d] .= parityCorrArr[d] .+ ( var2->( var->var .* parity2dLst[it][d][x][y] ).(var2) ).( tmpArr ) ;
					# # parityCorrArr[d] .= parityCorrArr[d] .+ ( var2->( var->var .* parity2dLst[it][d][x][y] ).(var2) ).( circshift.( circshift( parity2dLst[it][d], -x+1 ), -y+1 ) ) ;
					# # parity2dLst[it][d][x][y][n] * circshift.( circshift( parity2dLst[it][d], -x+1 ), -y+1 );
				# end
			# end
			# end
			# # @infiltrate
		# end
		# for it = 1 : itNum
			# parityCorrArr[d] .= parityCorrArr[d] .+ tmpArr2[it];
		# end
	# end
	# # @infiltrate
	# oFmain = "parityCorr_GOE";
	# oFName = fNameFunc( oFmain, attrLst, valLst, jld2Type; fMod );
	# save( oFName, "parityCorr", parityCorrArr );
# end

function parity_corr_GOE_arr_from_file( N, param_divide, itNum, seed; fMod = "", dim = nothing, fExt = jld2Type )
	divNum = param_divide[1];
	fMain = "parityArr";
	attrLst = deepcopy(attrLstBase);
	valLst = [N, param_divide, itNum, seed];
	if !isnothing(dim)
		attrLst = insert!( attrLst, 1, "dim" );
		valLst = insert!( valLst, 1, dim );
	end
	# @infiltrate
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod );
	@info(fName);
	tFull = @timed begin
		parity2dArr = load( fName, "parity2dArr" ) ;
	end
	
	# @infiltrate 
	parityCorrArr = zeros( size(parity2dArr) );
	parityArrSh = zeros( size(parity2dArr) );
	parityArrSummed = zeros( itNum, dim, N );
	# parityArrSummedSq = zeros( itNum, dim, N );
	parityArrTmp = zeros( size(parity2dArr) );
	dimsSummed = (3,4);
	# @infiltrate
	itNumLess = Int64( floor( itNum / 10 ) );
	
	for x = 1 : divNum, y = 1 : divNum
		print("x, y = ", string([x,y]), "\r")
		tFull = @timed begin
		circshift!(parityArrSh, parity2dArr, (0,0,x-1,y-1,0));
		parityArrSummed .= dropdims( sum( parity2dArr .* parityArrSh; dims = dimsSummed ); dims=dimsSummed);
		@view(parityCorrArr[:,:,x,y,:]) .= @view(parityCorrArr[:,:,x,y,:]) .+ parityArrSummed;
		end
		@info("read file: ");
		@info( timeMemStr( tFull.time, tFull.bytes ) );
		# @infiltrate
	end
	
	# @infiltrate
	oFmain = fNameZacCorr;
	oFName = fNameFunc( oFmain, attrLst, valLst, jld2Type; fMod );
	save( oFName, varNameZacCorr, parityCorrArr );
end

function parityCorr_GOE_AvgStd_fromFile( N, param_divide, itNum, seed; fMod = "", alpha = 0, dim = nothing, fExt = jld2Type )
	attrLst, valLst = fAttrValFunc( N, param_divide, itNum, seed, dim; alpha = alpha );
	fName = fNameFunc( fNameZacCorr, attrLst, valLst, fExt; fMod = fMod );
	parityCorrLst = load( fName, varNameZacCorr );
	
	parityCorrAvg = dropdims( mean( parityCorrLst; dims = 1 ); dims = 1 );
	parityCorrStd = dropdims( std( parityCorrLst; dims = 1 ); dims = 1 );
	
	fAvgName = fNameFunc( fNameZacCorrAvg, attrLst, valLst, fExt; fMod );
	fStdName = fNameFunc( fNameZacCorrStd, attrLst, valLst, fExt; fMod );
	
	save( fAvgName, varNameZacCorrAvg, parityCorrAvg );
	save( fStdName, varNameZacCorrStd, parityCorrStd );
end

function locLst_to_locDensity_fromfile( N, param_divide, itNum, seed; fMod = "", dim = nothing, fExt = jld2Type )
	attrLst = deepcopy(isnothing(dim) ? attrLstBase : attrLstDim);
	valLst = [N, param_divide, itNum, seed];
	if !isnothing(dim)
		insert!(valLst, 1, dim);
	end
	fName = fNameFunc( fNameDegGOE, attrLst, valLst, fExt; fMod );
	locLst = load( fName, varNamePosLoc );
	# @infiltrate
	
	thNum = Threads.nthreads();
	dimLstSh = [ deleteat!([1:dim;],d) for d = 1 : dim ];
	dimLst = [1:dim;];
	divProd = prod( param_divide );
	divLstSh = [ deleteat!(deepcopy(param_divide),d) for d = 1 : dim ];
	divProdShLst = prod.(divLstSh);
	shLst = [zeros(Int64,dim+1) for th = 1:thNum];
	xyLst = [ones(Int64,dim-1) for th = 1:thNum];
	
	locDenTmp = zeros( param_divide...,N );
	locDenTmpSh = [zeros( param_divide..., N ) for th = 1:thNum];
	locDenProd = [zeros( param_divide...,N ) for th = 1:thNum];
	locCorrTmp = zeros( param_divide..., N );
	locCorrLstAvg = [ zeros( divLstSh[d]..., N ) for d = 1 : dim ];
	locCorrLstSqAvg = [ zeros( divLstSh[d]..., N ) for d = 1 : dim ];
	locCorrLstStd = [ zeros( divLstSh[d]..., N ) for d = 1 : dim ];
	# @infiltrate
	
	for d=1:dim
		for it=1:itNum
			for n = 1 : N
				locDenTmp .= 0;
				for il = 1 : param_divide[dim]
					for id = 1 : size(locLst[itNum][dim][il,n],1)
						locDenTmp[@view(locLst[itNum][dim][il,n][id,:])...,n] += 1;
					end
				end
			end
			dLstTmp = dimLstSh[d];
			divLstTmp = divLstSh[d];
			@time begin
			for x = 1 : divLstTmp[1]
				print("it = ", string(it), ", d = ", d, ", x = ", x, "\r" )
				@time begin
				for y = 1 : divLstTmp[2]
					xyLst[Threads.threadid()][1] = x;
					xyLst[Threads.threadid()][2] = y;
					shLst[Threads.threadid()] .= 0;
					shLst[Threads.threadid()][dLstTmp] .= xyLst[Threads.threadid()] .- 1;
					circshift!(locDenTmpSh[Threads.threadid()],locDenTmp,shLst[Threads.threadid()]);
					locDenProd[Threads.threadid()] .= locDenTmp .* locDenTmpSh[Threads.threadid()];
					# @infiltrate
					selectDimLst(locCorrTmp, dLstTmp, xyLst[Threads.threadid()]) .+= dropdims( sum( locDenProd[Threads.threadid()]; dims = dLstTmp ); dims = Tuple(dLstTmp) ) ./ divProdShLst[d];
					# @infiltrate
				end
				end
			end
			end
			# @infiltrate
			locCorrLstAvg[d] .+= dropdims( mean( locCorrTmp; dims = d ); dims = d );
			locCorrLstSqAvg[d] .+= dropdims( mean( locCorrTmp.^2; dims = d ); dims = d );
		end
		locCorrLstAvg[d] ./= itNum;
		locCorrLstSqAvg[d] ./= itNum;
		locCorrLstStd[d] .= sqrt.( locCorrLstSqAvg[d] .- locCorrLstAvg[d].^2 );
	end
	
	fMainOut = fNameLocCorr;
	fNameOut = fNameFunc( fMainOut, attrLst, valLst, jld2Type; fMod );
	save( fNameOut, varNameLocCorrAvg, locCorrLstAvg, varNameLocCorrStd, locCorrLstStd );
end
