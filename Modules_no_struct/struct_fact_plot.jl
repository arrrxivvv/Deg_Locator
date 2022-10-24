using Plots
using JLD
using NPZ
include("structure_factors.jl")

# const fileType = ".pdf";
# const jldType = ".jld";
# const npyType = ".npy";

function plotNsave_abs( N, structFact, sliceLst, filenameMain, fileType )
	for n = 1:N
		for z in sliceLst
			figabs = contour( abs.( structFactLst[idStruct][:,:,z,n] ), fill=true );
			filename = string( filenameMain, "_abs", "_N_", N, "_n_", n, "_z_", z, fileType );
			savefig( figabs, filename );
		end
	end
end

function struct_plots()
	Nlst = [3, 5, 10, 15];
	# Nlst = [3, 5];
	param_dim = 3;
	param_divide = 50;
	param_step = 2*pi/param_divide;
	sliceLst = [1,10,20,25,30,40,50];
	# sliceLst = [20,30];
	
	filenameLst = [ "structFactFFT", "structFactTrue" ];
	
	
	for N in Nlst
		posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div_GUE( N, param_divide );
		
		structFactFFT = structure_factor_FFT( N, param_dim, param_divide, param_step, BfieldLst );
		structFactTrue = structure_factor_true( N, param_dim, param_divide, param_step, posLocLst, negLocLst );
		
		structFactLst = [structFactFFT, structFactTrue];
		
		for idStruct = 1:2
			varName = string( filenameLst[idStruct], "_N_", N );
			filename = string( varName, jldType );
			filename_npy = string( varName, npyType );
			save( filename, varName, structFactLst[idStruct] );
			npzwrite( filename_npy, structFactLst[idStruct] );
			# plotNsave_abs( N, structFactLst[idStruct], sliceLst, filenameLst[idStruct], fileType );
		end
		
		GC.gc();
	end	
end

function struct_plots_divide()
	Nlst = [3, 5];
	param_divide_lst = [25,50,100]
	param_dim = 3;
	# param_divide = 50;
	
	# sliceLst = [1,10,20,25,30,40,50];
	sliceLst = [20,30];
	
	filenameLst = [ "structFactFFT", "structFactTrue" ];
	
	for N in Nlst
		for param_divide in param_divide_lst
			sliceLst = convert.( Int64, round.([1,0.2*param_divide, 0.5*param_divide, 0.8*param_divide, param_divide]) );
			param_step = 2*pi/param_divide;
			posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div_GUE( N, param_divide );
			
			structFactFFT = structure_factor_FFT( N, param_dim, param_divide, param_step, BfieldLst );
			structFactTrue = structure_factor_true( N, param_dim, param_divide, param_step, posLocLst, negLocLst );
			
			structFactLst = [structFactFFT, structFactTrue];
			
			varName = "";
			for idStruct = 1:2
				varName = string( filenameLst[idStruct], "_N_", N, "_divide_", param_divide );
				filename = string( varName, jldType );
				filename_npy = string( varName, npyType );
				save( filename, varName, structFactLst[idStruct] );
				npzwrite( filename_npy, structFactLst[idStruct] );
				# plotNsave_abs( N, structFactLst[idStruct], sliceLst, filenameLst[idStruct], fileType );
			end
			GC.gc();
		end
	end	
end

function struct_plots_pol()
	Nlst = [3, 5, 10, 15];
	# Nlst = [3, 5];
	param_dim = 3;
	param_divide = 50;
	param_step = 2*pi/param_divide;
	sliceLst = [1,10,20,25,30,40,50];
	# sliceLst = [20,30];
	
	filenameMain = "structFact";
	filenamePolLst = ["pos", "neg"];
	
	for N in Nlst
		posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div_GUE( N, param_divide );
		# sizeLst = push!( param_divide * ones(param_dim), N );
		
		locLstPols = [posLocLst, negLocLst];
		
		# structFactLst = [ zeros( sizeLst... ) for idPols = 1:2 ];
		for idPol = 1:2
			polarity = (-1)^idPol;
			structFactPol = structure_factor_pol( N, param_dim, param_divide, param_step, locLstPols[idPol], polarity );
			varName = string( filenameMain, "_", filenamePolLst[idPol], "_N_", N );
			filename = string( varName, jldType );
			filename_npy = string( varName, npyType );
			save( filename, varName, structFactPol );
			npzwrite( filename_npy, structFactPol );
		end
		GC.gc();
	end	
end