using Utils
using DegLocatorDiv
using JLD
using NPZ

function fFunc_stat( posLocLst, negLocLst, param_divide_num, N, num_it )
	locLst = [ posLocLst, negLocLst ];
	
	idPosNegNum = 4;
	fFuncLst = zeros( N, idPosNegNum, param_divide_num, param_divide_num );	
	
	for n = 1:N
		for it = 1 : num_it
			idPosNeg = 1;
			for idPol1 = 1:2
				println("sort locs1");
				@time locs1 = sort( locLst[idPol1][it][n][:,1] );
				for idPol2 = 1:2
					println("sort locs2");
					@time locs2 = sort( locLst[idPol2][it][n][:,1] );
					fFuncTemp = zeros( param_divide_num, param_divide_num );
					println("place onto grid");
					@time begin
						for ix = 1 : size(locs1)[1]
							for iy = 1 : size(locs2)[1]
								fFuncTemp[ locs1[ix], locs2[iy] ] += 1;
							end
						end
					end
					fFuncLst[ n, idPosNeg, :, : ] += fFuncTemp;
				end
				idPosNeg += 1;
			end
		end
		fFuncLst[ n, :, :, : ] = fFuncLst[ n, :, :, : ] ./ num_it;
	end
	
	return fFuncLst;
end

function fFunc_stat_from_file( N, param_divide, num_it, seedFedStart )
	fileNameAttr = fileNameAttrFunc( N, param_divide, num_it, seedFedStart );
	fileNameMain = "deg"; 
	oFileName = string( fileNameMain, "_", fileNameAttr );
	varFileName = string( oFileName, jldType );
	data = load(varFileName);
	
	param_divide = Int64.(param_divide);
	if param_divide isa Array{Int64}
		param_divide_num = param_divide[1];
	else
		param_divide_num = param_divide;
	end
	
	fFuncLst = fFunc_stat( data["posLocLst"], data["negLocLst"], param_divide_num, N, num_it );
	
	oFileNameMain = "ffunc";
	fileNameNpy = string( oFileNameMain, "_",  fileNameAttr, npyType );
	npzwrite( fileNameNpy, fFuncLst );
	
	return fFuncLst;
end
