N=10;
num_it = 5;

param_dim = 3;
param_min = fill(0.0,param_dim);
param_max = fill(2*pi,param_dim);
param_divide = [80,16,16];
ratio = [1,1,1];

degObj = DegObj( N, param_dim, param_divide, param_min, param_max );

H_GUElst = DegLocatorDiv.HlstFunc( H_GUE, degObj.param_dim, degObj.N, false );

HmatFun = (H,xlst)-> Hmat_3comb_ratio!(H,xlst,H_GUElst, ratio);

@time begin
	Threads.@threads for ind in degObj.posLst
		HmatFun( degObj.HmatLst[ind], degObj.param_mesh[ind] );
	end
end

@time begin
	for it = 1 : num_it 
	   Threads.@threads for ii in degObj.posLst		 
		   degObj.Elst[ii], degObj.vecLst[ii] = eigen( degObj.HmatLst[ii] );
	   end
	end
end

@time begin
   Threads.@threads for ii in degObj.posLst		 
	   degObj.Elst[ii], degObj.vecLst[ii] = eigen( degObj.HmatLst[ii] );
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst		 
	   callZheevr!( degObj.HmatLst[ii], degObj.Elst[ii], degObj.vecLst[ii] );
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst		 
		matTmp = copy( degObj.HmatLst[ii] );
	   degObj.Elst[ii], degObj.vecLst[ii] = LAPACK.geevx!('B', 'N', 'V', 'N', matTmp)[[2,4]];
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst		 
		matTmp = copy( degObj.HmatLst[ii] );
	   LAPACK.geevx!('B', 'N', 'V', 'N', matTmp)[[2,4]];
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst
	   degObj.Elst[ii], degObj.vecLst[ii] = LAPACK.geevx!('B', 'N', 'V', 'N', degObj.HmatLst[ii])[[2,4]];
   end
end

ElstCx = [ zeros( Complex{Float64}, N ) for i1 in degObj.posLst ];

@time begin
   Threads.@threads for ii in degObj.posLst		 
	   degObj.Elst[ii], degObj.vecLst[ii] = eigen( degObj.HmatLst[ii] );
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst		 
	   ElstCx[ii], degObj.vecLst[ii] = LAPACK.geevx!('B', 'N', 'V', 'N', degObj.HmatLst[ii])[[2,4]];
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst		 
	   ElstCx[ii], degObj.vecLst[ii] = LAPACK.syevr!('B', 'N', 'V', 'N', degObj.HmatLst[ii])[[2,4]];
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst	
		matTmp = Hermitian( degObj.HmatLst[ii] );
	   ElstCx[ii], degObj.vecLst[ii] = LAPACK.syevr!('V', 'A', matTmp.uplo, matTmp.data, 0.0, 0.0, 0, 0, -1.0);
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst		 
	   degObj.Elst[ii], degObj.vecLst[ii] = eigen( degObj.HmatLst[ii] );
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst		 
	   eLstTmp, vecLstTmp = eigen!( degObj.HmatLst[ii] );
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst		 
	   factTmp = eigen( degObj.HmatLst[ii] );
   end
end

testEigArr = Array{Eigen{Complex{Float64},Float64,Array{Complex{Float64},2},Array{Float64,1}},3}(undef,param_divide...);

@time begin
   Threads.@threads for ii in degObj.posLst		 
	   factTmp = eigen( degObj.HmatLst[ii] );
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst		 
	   testEigArr[ii] = eigen( degObj.HmatLst[ii] );
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst		 
	   testEigArr[ii] = eigen!( degObj.HmatLst[ii] );
   end
end

@time begin
   Threads.@threads for ii in itLnNum
	   eigen!( mat40 );
   end
end

@time begin
   Threads.@threads for ii in itLnNum
	   LAPACK.geevx!(permute ? (scale ? 'B' : 'P') : (scale ? 'S' : 'N'), 'N', 'V', 'N', mat40)[[2,4]];
   end
end

@time begin
   Threads.@threads for ii in itLnNum
		matTmp = copy(mat40);
	   val, vec = LAPACK.geevx!(permute ? (scale ? 'B' : 'P') : (scale ? 'S' : 'N'), 'N', 'V', 'N', matTmp)[[2,4]];
   end
end

matTmpArr = [similar(mat40) for ii = 1 : itLnNum];

@time begin
   Threads.@threads for ii in 1:itLnNum
		matTmpArr[ii] .= mat40;
		# matTmpH = Hermitian(matTmp);
	   val, vec = LAPACK.syevr!('V', 'A', Hermitian(matTmpArr[ii]).uplo, Hermitian(matTmpArr[ii]).data, 0.0, 0.0, 0, 0, -1.0);
   end
end

@time begin
   Threads.@threads for ii in degObj.posLst
		matTmp .= mat40;
		matTmpH = Hermitian(matTmp);
	   LAPACK.syevr!('V', 'A', matTmpH.uplo, matTmpH.data, 0.0, 0.0, 0, 0, -1.0);
   end
end

matTmp = [ similar(mat40) for 1 : itLnNum ];

@time begin
   Threads.@threads for ii in 1:itLnNum
		matTmp .= mat40;
		matTmpH = Hermitian(matTmp);
	   eigen(matTmpH);
   end
end

@time begin
   Threads.@threads for ii in 1:itLnNum
		matTmp .= mat40;
		matTmpH = Hermitian(matTmp);
	   eigen!(matTmpH);
   end
end



  # val, degObj.vecLst[ii] = eigen( degObj.HmatLst[ii] );
		   # degObj.Elst[ii] = Real.(val);
		   
		   
 isuppz = similar(A, BlasInt, 2*n)
work   = Vector{$elty}(undef, 1)
lwork  = BlasInt(-1)
rwork  = Vector{$relty}(undef, 1)
lrwork = BlasInt(-1)
iwork  = Vector{BlasInt}(undef, 1)
liwork = BlasInt(-1)
info   = Ref{BlasInt}()

w = similar(mat40,Float64,40);
Z = similar(mat40);

mat40 = rand()

itLnNum = 32000;

matTmpArr = [ similar(mat40) for it = 1 : itLnNum ];
ZLst = [ similar(mat40) for it = 1 : itLnNum ];
wLst = [ similar(mat40,Float64,40) for it = 1 : itLnNum ];

@time begin
   Threads.@threads for ii in 1:itLnNum
		matTmpArr[ii] .= mat40;
		callZheevr!( matTmpArr[ii], wLst[ii], ZLst[ii] );
   end
end

@time begin
   Threads.@threads for ii in 1:itLnNum
		matTmpArr[ii] .= mat40;
		Fobj = eigen!( matTmpArr[ii] );
		wLst[ii] .= Fobj.values;
		ZLst[ii] .= Fobj.vectors;
   end
end


function getWorkLn( N )
	mat = rand(ComplexF64, N, N);
	mat = mat + mat';
	wN = similar( mat, Float64, N );
	ZN = similar( mat );
	
	@time eigenZheevr!( mat, wN, ZN );
end

function testZheevTime( N, itLnNum; preWork = false )
	mat = rand(ComplexF64, N, N);
	mat = mat + mat';
	
	matTmpArr = [ similar(mat) for it = 1 : itLnNum ];
	ZLst = [ similar(mat) for it = 1 : itLnNum ];
	wLst = [ similar(mat,Float64,N) for it = 1 : itLnNum ];
	
	if preWork
		lwork, lrwork, liwork = eigenPrework!( matTmpArr[1], wLst[1], Zlst[1] );
		work = [ Vector{ComplexF64}(undef, lwork) for ii = 1 : itLnNum ];
		rwork = [ Vector{Float64}(undef, lrwork) for ii = 1 : itLnNum ];
		iwork = [ Vector{BlasInt}(undef, liwork) for ii = 1 : itLnNum ];
	end
	
	@time begin
		Threads.@threads for ii in 1:itLnNum
			matTmpArr[ii] .= mat;
			if preWork
				callZheevr!( matTmpArr[ii], wLst[ii], ZLst[ii], lwork, lrwork, liwork; work=work, rwork=rwork, iwork=iwork );
			else
				callZheevr!( matTmpArr[ii], wLst[ii], ZLst[ii] );
			end
		end
	end
	
	# @time eigenZheevr!( mat, wN, ZN );
end
