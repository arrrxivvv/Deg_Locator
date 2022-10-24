using Utils
using DegLocatorDiv
using JLD
using NPZ

function flipFfunc( fFunc, param_divide_num )
	fFuncFlip = reverse( fFunc; dims=1 );
	fFuncFlip = reverse( fFuncFlip; dims=2 );
	fFuncFlip = circshift( fFuncFlip, (1,1) );
	
	halfDivide = Int64( param_divide_num );
	fFuncFlip2 = circshift( fFunc, (halfDivide, halfDivide) );
	fFuncFlip3 = circshift( fFuncFlip, (halfDivide, halfDivide) );
	
	return [fFuncFlip, fFuncFlip2, fFuncFlip3];
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

function fFunc_stat( posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = 2 )
	locLst = [ posLocLst, negLocLst ];
	
	idPosNegNum = 4;
	fFuncLst = zeros( N, idPosNegNum, param_divide_num, param_divide_num );	
	
	fFunc1 = zeros( param_divide_num );
	fFunc2 = zeros( param_divide_num );
	
	fFuncTemp = zeros( param_divide_num, param_divide_num );
	
	for n = 1:N
		for it = 1 : num_it
			# @time begin
				idPosNeg = 1;	
				for idPol1 = 1:idPolEnd
					locs1 = locLst[idPol1][it][n][:,1];
					for idPol2 = 1:idPolEnd
						locs2 = locLst[idPol2][it][n][:,1];
						fFuncTemp .= 0;
						for ix in locs1, iy in locs2
							fFuncTemp[ix,iy] += 1;
						end
						fFuncLst[ n, idPosNeg, :, : ] += fFuncTemp;
						fFuncFlipLst = flipFfunc( fFuncTemp, param_divide_num );
						for fFuncFlip in fFuncFlipLst
							fFuncLst[ n, idPosNeg, :, : ] += fFuncFlip;
						end
						idPosNeg += 1;
					end
				end
			# end
		end
		symFact = 4;
		fFuncLst[ n, :, :, : ] = fFuncLst[ n, :, :, : ] ./ num_it ./ symFact;
	end
	
	return fFuncLst;
end

function fFunc_diff_stat( posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = 2 )
	locLst = [ posLocLst, negLocLst ];
	
	idPosNegNum = 4;
	fFuncLst = zeros( N, idPosNegNum, param_divide_num, param_divide_num );
	fFuncTemp = zeros( param_divide_num, param_divide_num );
	
	for n = 1:N
		for it = 1 : num_it
			idPosNeg = 1;
			for idPol1 = 1:idPolEnd
				locs1 = locLst[idPol1][it][n][:,1];
				for idPol2 = 1:idPolEnd
					fFuncTemp .= 0; 
					locs2 = locLst[idPol2][it][n][:,1];
					for ix in locs1, iy in locs2
						iDiff = Int64.( abs.( ix .- iy ) ) .+ 1;
						fFuncTemp[ Tuple(iDiff)... ] += 1;
					end
					fFuncLst[ n, idPosNeg, :, : ] += fFuncTemp;
					idPosNeg += 1;
				end
			end
		end
		fFuncLst[ n, :, :, : ] = fFuncLst[ n, :, :, : ] ./ num_it;
	end
	
	return fFuncLst;
end

function fFunc_stat_from_file( N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3, fileOutNameMod = "" )
	fileNameAttr = fileNameAttrFunc( N, param_divide, num_it, seedFedStart );
	if fileNameMod != ""
		fileNameAttr = string( fileNameAttr, "_", fileNameMod );
	end
	fileNameMain = "deg";
	if param_dim == 2
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
	
	@time fFuncLst = fFunc_stat( data["posLocLst"], data["negLocLst"], param_divide_num, N, num_it; idPolEnd = idPolEnd );
	
	oFileNameMain = "ffunc";
	if fileOutNameMod != ""
		fileNameAttr = string( fileNameAttr, "_", fileOutNameMod );
	end
	fileNameNpy = string( oFileNameMain, "_",  fileNameAttr, npyType );
	npzwrite( fileNameNpy, fFuncLst );
	
	return fFuncLst;
end
