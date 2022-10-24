using LinearAlgebra
using ArgParse
using JLD2, FileIO
include("Umat.jl")
include("Chern.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-o"
            help = "output file names"
			nargs = *
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

N = 3;
V1 = 10;
V2 = 10;
tau = 1;
ang_divide = 20;
ang_ln = ang_divide+1;

u=3.4;

V_step = 1;
V_min = 1;
V_max = 10;

V1_lst = [V_min:V_step:V_max];
V2_lst = [V_min:V_step:V_max];

chernLst = zeros( Complex(Float64), ( V1_lst.shape, V2_lst.shape, N );
quasiELst = zeros( Complex(Float64), ( V1_lst.shape, V2_lst.shape, ang_ln, ang_ln, N ) );
eVecLst = zeros( Complex(Float64), ( V1_lst.shape, V2_lst.shape, ang_ln, ang_ln, N, N ) );

for indV1 = 1:length(V1_lst)
	V1 = V1_lst[indV1];
	for indV2 = 1:length(V2_lst)
		V2 = V2_lst[indV2];
		print(V1, V2);
		U_mat_flq = ( ang1, ang2) -> ( U_mat_flq_raw( N, V1, V2, tau, ang1, ang2 );
		(chernLstTmp, quasiELstTmp, eVecLstTmp) = Chern( N, U_mat_flq, ang_divide, isFlq = true, isBandMatch = True );
		chernLst[indV1,indV2,:] .= chernLstTmp;
		quasiELst[indV1,indV2,:,:,:] .= quasiELstTmp;
		eVecLst[indV1,indV2,:,:,:,:] .= eVecLstTmp;
		printArrPres( chernLst[indV1,indV2,:] );
	end	
end

if parseargs( "o" )
	filenames = parseargs( "o" );
	@save filename[1] [V1_lst;V2_lst];
	@save filename[2] round.( real.(chernLst); digits=2 );
	if length(filename) >= 4
		@save filename[3] quasiELst;
		@save filename[3] eVecLst;
	end
end
