using Utils
using DegLocatorDiv
using JLD
using NPZ

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
		fFuncLst = fFuncLst ./ num_it;
	end
	
	return fFuncLst;
end

function fFunc_stat( posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = 2 )
	locLst = [ posLocLst, negLocLst ];
	
	idPosNegNum = 4;
	fFuncLst = zeros( N, idPosNegNum, param_divide_num, param_divide_num );	
	
	fFunc1 = zeros( param_divide_num );
	fFunc2 = zeros( param_divide_num );
	
	for n = 1:N
		for it = 1 : num_it
			@time begin
				idPosNeg = 1;	
				for idPol1 = 1:idPolEnd
					locs1 = locLst[idPol1][it][n][:,1];
					fFunc1 .= 0;
					for ix in locs1
						fFunc1[ix] += 1;
					end
					for idPol2 = 1:idPolEnd
						locs2 = locLst[idPol2][it][n][:,1];			
						fFunc2 .= 0;
						for iy in locs2
							fFunc2[iy] += 1;
						end
						fFuncTemp = fFunc1 .* fFunc2';
						fFuncLst[ n, idPosNeg, :, : ] += fFuncTemp;
						idPosNeg += 1;
					end
				end
			end
		end
		@time fFuncLst[ n, :, :, : ] = fFuncLst[ n, :, :, : ] ./ num_it;
	end
	
	return fFuncLst;
end

function fFunc_diff_stat( posLocLst, negLocLst, param_divide_num, N, num_it; idPolEnd = 2 )
	locLst = [ posLocLst, negLocLst ];
	
	idPosNegNum = 4;
	fFuncLst = zeros( N, idPosNegNum, param_divide_num, param_divide_num );	
	
	for n = 1:N
		for it = 1 : num_it
			idPosNeg = 1;
			for idPol1 = 1:idPolEnd
				locs1 = sort( locLst[idPol1][it][n][:,1] );
				fFunc1 = zeros( param_divide_num );
				for ix in locs1
					fFunc1[ix] += 1;
				end
				for idPol2 = 1:idPolEnd
					locs2 = sort( locLst[idPol2][it][n][:,1] );
					fFunc2 = zeros( param_divide_num );
					for iy in locs2
						fFunc2[iy] += 1;
					end
					fFuncTemp = fFunc1 .* fFunc2';
					fFuncLst[ n, idPosNeg, :, : ] += fFuncTemp;
					idPosNeg += 1;
				end
			end
		end
		fFuncLst[ n, :, :, : ] = fFuncLst[ n, :, :, : ] ./ num_it;
	end
	
	return fFuncLst;
end

function fFunc_stat_from_file( N, param_divide, num_it, seedFedStart, fileNameMod = "", param_dim = 3 )
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
	
	@time fFuncLst = fFunc_stat( data["posLocLst"], data["negLocLst"], param_divide_num, N, num_it );
	
	oFileNameMain = "ffunc";
	# if fileNameMod != ""
		# fileNameAttr = string( fileNameAttr, "_", fileNameMod );
	# end
	fileNameNpy = string( oFileNameMain, "_",  fileNameAttr, npyType );
	npzwrite( fileNameNpy, fFuncLst );
	
	return fFuncLst;
end
