using LinearAlgebra
# using DegLocatorDiv

function mapDeg( avgNum, N, param_divide, seed; fileNameMod = "", dim = nothing, pow = 2 )
	fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed; dim = dim );
	fileNameAttr = string( fileNameAttr, "_", fileNameMod );
	degFileName = string( "deg", "_", fileNameAttr );
	varFileName = string( degFileName, jldType );
	@info(varFileName);
	posLocLst = load(varFileName, "posLocLst");
	negLocLst = load(varFileName, "negLocLst");
	posNlst = load(varFileName, "posNlst");
	
	H_GUE_lst = load(varFileName, "H_GUE_lst");
	
	n=2;
	it=1;
	
	divide=80;
	param_step = 2*pi/80;
	
	vecParam = zeros(2*dim);
	locSh = zeros( 2*dim );
	vecDim = 2*dim;
	nearGapMat = zeros(ones(Int64,vecDim)*3);
	
	Hsh = zeros( Complex{Float64}, N, N );
	
	for ii= 1 : size( posLocLst[it][n],1 )
		locP = posLocLst[it][n][ii,:];
		for d = 1 : dim
			vecParam[2*dim-1] = cos( locP[d] );
			vecParam[2*dim] = sin( locP[d] );
		end
		
		vecParam /= norm(vecParam);
		for iSh in CartesianIndices( nearGapMat )
			iArr = [iSh[d] for d = 1:vecDim];
			iArr -= 2;
			locSh = vecParam + iArr * param_step;
			Hsh .= 0;
			for ih in eachindex(H_GUE_lst[it])
				Hsh += locSh[ih] * H_GUE_lst[it][ih];
			end
		end
	end
end

function locLstPurify( posLocLst, negLocLst, posNlst )
	itNum = size(posLocLst, 1);
	N = size(posLocLst, 2);
	nLst = ones(Int64, itNum, N);
	pureNlst = posNlst[:,2:end] - posNlst[:,1:end-1];
	
	posLocLstPure = [ [ zeros( Int64, pureNlst[it,n], 3 ) for n = 1:N-1 ] for it = 1:itNum ];
	negLocLstPure = [ [ zeros( Int64, pureNlst[it,n], 3 ) for n = 1:N-1 ] for it = 1:itNum ];
	
	locLstPols = [posLocLst, negLocLst];
	locLstPolsPure = [posLocLstPure, negLocLstPure];
	polInvLst = [1,2];
	
	Threads.@threads for it = 1 : itNum
		for n = 2 : N
			for pol = 1:2
				if n == 1
					locLstPolsPure[pol][it][1,:] = locLstPols[pol][it][1,:];
					continue ;
				end
				polInv = polInvLst[pol];
				
				lnThis = posNlst[itNum,n];
				lnPrev = pureNlst[itNum,n-1];
				idPrev = 1;
				idPure = 1;
				for id = 1 : lnThis
					locTested = locLstPos[pol][it][n][id,:];
					if idPrev > lnPrev
						locLstPolsPure[pol][it][n][idPure,:] = locTested;
						idPure++;
					elseif locTested != locLstPolsPure[polInv][it][n][idPrev,:]
						locLstPolsPure[pol][it][n][idPure,:] = locTested;
						idPure++;
					else
						idPrev++;
					end
				end
			end
		end
	end
	
	return locLstPolsPure[1], locLstPolsPure[2], pureNlst;
end
