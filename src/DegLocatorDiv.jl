module DegLocatorDiv

using MKL
using UMat
using Utils
using DivB
using EigCustom
using ThreadedArrays
import LinearAlgebra: BlasFloat, BlasInt

# @enum EnumDegMethod degMethodFlux=1 degMethodRootFind=2
# @enum EnumSaveMem memNone=1 memEig=2 memEigLink=3
# @enum EnumRootEig eigFull=1 eigLanczos=2
@enum EnumSaveMem memNone=1 memEig=2 memEigLink=3 rootFind=4 rootFindLanczos=5
export EnumSaveMem, memNone, memEig, memEigLink, rootFind, rootFindLanczos

include("DegParams.jl");
export DegParams, degParamsNonPeriodic, makeArrOverGrid, linIdFromIdVec, wrapIdVec!

export DegObj
struct DegObj
	param_divide::Vector{Int64};
	param_dim::Int64;
	N::Int64;
	Bfield_ln::Int64;
	posLst::AbstractArray{CartesianIndex{N}, N} where N;
	posLstSlice;
	nonPeriodic::Bool;
	isStrongGC::Bool;
	enumSaveMem::EnumSaveMem;
	
	param_min::Vector{Float64};
	param_max::Vector{Float64};
	param_step::Vector{Float64};
	param_grids::Vector{ Vector{Float64} };
	param_mesh::Array{ Vector{Float64} };
	
	HmatLst::Array{Array{Complex{Float64}}};
	linkLst::Vector{ Array{ Vector{Complex{Float64}} } };
	BfieldLst::Vector{ Array{ Vector{Complex{Float64}} } };
	divBLst::Array{ Vector{Complex{Float64}} };
	Elst::Array{ Vector{Float64} };
	vecLst::Array{ Array{Complex{Float64}} };
	vecLstPrev::Array{ Array{Complex{Float64}} };
	vecLst1st::Array{ Array{Complex{Float64}} };
	non0LstRe::Array{ Vector{Float64} };
	divBposArr::Array{ Bool };
	divBnegArr::Array{ Bool };
	
	lwork::BlasInt;
	lrwork::BlasInt;
	liwork::BlasInt;
	workLst::Array{ComplexF64};
	rworkLst::Array{Float64};
	iworkLst::Array{BlasInt};
	
	linkRatio::Vector{ Array{ Vector{Complex{Float64}} } };
	linkLstSh::Array{ Vector{Complex{Float64}} };
	BfieldSh::Array{ Vector{Complex{Float64}} };
	
	lnSimp::Int64;
	lnExtra::Int64;
	lnSimpAll::Int64;
	ptLst::ThrArray{Float64};
	ptLstEnd::ThrArray{Float64};
	ptLstTmp::ThrArray{Float64};
	ptIntLst::ThrArray{Int64};
	shSimpCartLst::Array{AbstractArray{CartesianIndex{N}, N},2} where{N} ;# where {N,M};
	shMeshLst::Array{Array{Vector{Float64}},2};
	valLst::ThrArray{Float64};
	valLstTmp::ThrArray{Float64};
	ixLst::ThrArray{Int64};
	shPt1::Vector{Int64};
	shPtOther::Matrix{Int64};
	Htmp::ThrArray{Complex{Float64}};
	vLstTmp::ThrArray{Complex{Float64}};
	eLstTmp::ThrArray{Float64};
	auxLst;
	auxLstTmp;
	eDiffLst::Array{Float64};
	locLstFromSimp::Array{Array{Float64}};
	locLstHash::Array{Vector{Vector{Float64}}};
	# locLstNoDup::Vector{Vector{Vector{Float64}}};
	thresNM::Float64;
	thresEDeg::Float64;
	thresDegCollision::Float64;
	
	function DegObj(N,param_dim,param_divide, param_min, param_max; nonPeriodic = false, isStrongGC = false, enumSaveMem = memNone, thresNM = 1e-10, thresEDeg = 1e-7, thresDegCollision = 1e-7)
		posLst = CartesianIndices(Tuple(param_divide));
		posLstSlice = CartesianIndices( selectdim( posLst, param_dim, 1 ) );
	
		Bfield_ln = Int64( round( param_dim * (param_dim-1) / 2 ) );
		param_step = ( param_max .- param_min ) ./ param_divide;
		param_grids = [ collect( range( param_min[iDim], param_max[iDim] - param_step[iDim], length = param_divide[iDim] ) ) for iDim = 1:param_dim ];
		param_mesh = [ [ param_grids[j][ind[j]] for j = 1:param_dim ] for ind in posLst ];
		@debug "param_mesh" param_mesh[1,1,2]
		
		if enumSaveMem == rootFind || enumSaveMem == rootFindLanczos
			Elst = [ zeros( Float64, N ) for i1 in posLst ];
			vecLst = [];
			vecLstPrev = [];
			vecLst1st = [];
		elseif enumSaveMem == memNone || enumSaveMem == memEig
			Elst = [ zeros( Float64, N ) for i1 in posLst ];
			vecLst = [ zeros( Complex{Float64}, N, N ) for i1 in posLst ];
			vecLstPrev = [];
			vecLst1st = [];
		else
			Elst = [ zeros( Float64, N ) for i1 in posLstSlice ];
			vecLst = [ zeros( Complex{Float64}, N, N ) for i1 in posLstSlice ];
			vecLstPrev = [ zeros( Complex{Float64}, N, N ) for idSlice in posLstSlice ];
			vecLst1st = [ zeros( Complex{Float64}, N, N ) for idSlice in posLstSlice ];
		end
		
		if enumSaveMem == memNone
			HmatLst = [ zeros( Complex{Float64}, N, N ) for i1 in posLst ];
		else
			HmatLst = [ zeros( Complex{Float64}, N, N ) for i1 in posLstSlice ];
		end
		
		lwork, lrwork, liwork, workLst, rworkLst, iworkLst = eigenWorkThrdInit!( HmatLst[1], Elst[1], HmatLst[1] );
		
		if enumSaveMem != rootFind && enumSaveMem != rootFindLanczos
			linkLst = [ [ zeros( Complex{Float64}, N ) for pos in posLst ] for dim = 1:param_dim ];
			BfieldLst = [ [ zeros( Complex{Float64}, N ) for pos in posLst ] for l = 1 : Bfield_ln ];
			divBLst = [ zeros(Float64,N) for it in posLst ];
			non0LstRe = [ zeros(Float64,N) for it in posLst ];
			divBposArr = fill( false, param_divide...,N );
			divBnegArr = fill( false, param_divide...,N );
			
			linkRatio = [ [ zeros( Complex{Float64}, N ) for pos in posLst ] for id = 1:2 ];
			linkLstSh = [ zeros(Complex{Float64}, N) for pos in posLst ];
			BfieldSh = [ zeros(Complex{Float64},N) for pos in posLst ];
			lnSimp = 1;
			lnExtra = 0;
			lnSimpAll = lnSimp + lnExtra;
			ptLst = thrArr_empty();
			ptLstEnd = thrArr_empty();
			ptLstTmp = thrArr_empty();
			ptIntLst = thrArr_empty(Int64);
			shSimpCartLst = Matrix{AbstractArray{CartesianIndex{param_dim},param_dim}}(undef,0,0);
			shMeshLst = Array{Array{Vector{Float64}},2}(undef,0,0);
			valLst = thrArr_empty();
			valLstTmp = thrArr_empty();
			ixLst = thrArr_empty(Int64);
			shPt1 = [];
			shPtOther = Matrix{Int64}(undef,0,0);
			Htmp = thrArr_empty(ComplexF64);
			eLstTmp = thrArr_empty();
			vLstTmp = thrArr_empty(ComplexF64);
			auxLst = nothing;
			auxLstTmp = nothing;
			eDiffLst = [];
			locLstFromSimp = [];
			locLstHash = [];
			# locLstNoDup = [];
		else
			linkLst = [];
			BfieldLst = [];
			divBLst = [];
			non0LstRe = [];
			divBposArr = [];
			divBnegArr = [];
			linkRatio = [];
			linkLstSh = [];
			BfieldSh = [];
			lnSimp = 2^( param_dim - 1 );
			lnExtra = ( param_dim == 3 ? 1 : 0 );
			lnSimpAll = lnSimp + lnExtra;
			ptIntLst = threaded_zeros( Int64, param_dim, param_dim+1 );
			ptLst = threaded_zeros( param_dim, param_dim+1 );
			ptLstEnd = similar( ptLst );
			ptLstTmp = similar( ptLst );
			valLst = threaded_zeros(param_dim+1);
			valLstTmp = similar( valLst );
			ixLst = threaded_zeros( Int64, param_dim+1 );
			shPt1 = append!(zeros(Int64, lnSimp),ones(Int64,lnExtra));
			shPtOther = ones( Int64, param_dim, lnSimp + lnExtra );
			Htmp = threaded_zeros( Complex{Float64}, N, N );
			if enumSaveMem == rootFindLanczos
				vLstTmp = threaded_zeros( Complex{Float64}, N, N );
				auxLst = threaded_fill( auxLanczos, (0, zeros(ComplexF64,N)), (param_dim+1,) );
				auxLstTmp = deepcopy(auxLst);
			else
				vLstTmp = thrArr_empty(ComplexF64);
				auxLst = thrArr_empty();
				auxLstTmp = thrArr_empty();
			end
			eLstTmp = threaded_zeros( N );
			id1 = 1; id2 = id1 + 1;
			for it = 2 : lnSimp
				shPtOther[id1,it] = -1;
				shPtOther[id2,it] = -1;
				id2 += 1;
				if id2 > param_dim
					id1 += 1;
					id2 = id1 + 1;
				end
			end
			shSimpCartLst = Array{AbstractArray{CartesianIndex{param_dim},param_dim}, 2}(undef, param_dim+1, lnSimpAll);
			shMeshLst = [ deepcopy(param_mesh) for id = 1:(param_dim+1), iSimp = 1:lnSimpAll ];
			# CircShiftedArray{CartesianIndex{param_dim},param_dim,typeof(posLst)}
			shArr = zeros(Int64, param_dim);
			for iSimp = 1 : lnSimpAll
				shArr .= shPt1[iSimp];
				shSimpCartLst[1,iSimp] = ShiftedArrays.circshift( posLst, Tuple(shArr) );
				for iCart in posLst
					shMeshLst[1,iSimp][iCart] .-= shArr.*param_step;
				end
				for id = 1 : param_dim
					shArr .= 0;
					shArr[id] = shPtOther[id,iSimp];
					shSimpCartLst[id+1,iSimp] = ShiftedArrays.circshift( posLst, shArr );
					for iCart in posLst
						shMeshLst[id+1,iSimp][iCart] .-=  shArr.*param_step;
					end
				end
			end
			eDiffLst = zeros( param_divide..., lnSimpAll );
			locLstFromSimp = [ zeros( param_dim, lnSimpAll ) for loc in posLst ];
			locLstHash = [ [] for ii in posLst ];
			# locLstNoDup = [[] for ii = 1 : N-1 ];
		end
		new(param_divide, param_dim, N, Bfield_ln, posLst, posLstSlice, nonPeriodic, isStrongGC, enumSaveMem, param_min, param_max, param_step, param_grids, param_mesh, HmatLst, linkLst, BfieldLst, divBLst, Elst, vecLst, vecLstPrev, vecLst1st, non0LstRe, divBposArr, divBnegArr, lwork, lrwork, liwork, workLst, rworkLst, iworkLst, linkRatio, linkLstSh, BfieldSh, lnSimp, lnExtra, lnSimpAll, ptLst, ptLstEnd, ptLstTmp, ptIntLst, shSimpCartLst, shMeshLst, valLst, valLstTmp, ixLst, shPt1, shPtOther, Htmp, vLstTmp, eLstTmp, auxLst, auxLstTmp, eDiffLst, locLstFromSimp, locLstHash, thresNM, thresEDeg, thresDegCollision );
	end
end

include("divBIOFuncs.jl")
export fileNameAttrFunc, fNameAttrLstFunc, fNameFunc, fAttrOptLstFunc

include("degLocator_funcs.jl")
export locator_div, locator_div_GUE, locator_div_sin3, locLstDistill, locator_div_GUE_scale, locator_div_GUE_ratio, locator_div_GOE_ratio

include("degLocator_funcs_rootfind.jl")

# include()

include("distillLocs_func.jl")
export distillLocsFromFile, distillLocN, whichLocsFromFile, distillLocsFromWhich, locLstPurify, locLstPurify_detailedOutput, locNvarFromFile, parityGOE_resave_fromFile, parityAvg_fromFile, collisionPurifyPerLevel!, collisionPurifyFromFile

include("divBProfile_funcs.jl")
export divB_profile, divB_profile_GOE_3d
degOptAttrLst = ["scale", "ratio", "alpha"];
degOptDefaultLst = [1, nothing, 0];
rtFndAttrLst = ["thresNM", "thresEDeg"];
rtFndDefaultLst = [1e-7, 1e-5];

include("deltaN_funcs.jl")
export deltaN

include("deltaNAvg_funcs.jl")
export deltaN_avg, deltaN_avg_fromFile, deltaN_avg_lst, deltaN_avg_lst_fromFile, varG_fromFile, varG_lst_fromFile, deltaN_var_lst_fromFile, deltaN_var_fromFile, deltaNCum_avg_fromFile, deltaNCum_avg_lst_fromFile, deltaN_var_fromFileDirect

include("fFunc_stat.jl")
export FfuncObj, fFunc_stat, fFunc_stat_from_file,  fFunc_stat_across_from_file, fFunc_stat_diff_from_file, fFunc_stat_diff_full_from_file, fFunc_stat_diff_full_3d_from_file, fFunc_stat_diff_mod_3d_from_file, parity_corr_GOE_from_file, parity_corr_GOE_arr_from_file, locLst_to_locDensity_fromfile, parityCorr_GOE_AvgStd_fromFile

end
