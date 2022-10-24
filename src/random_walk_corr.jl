using Utils
using JLD
using NPZ

function monopole_wrap( locLst, Nparam )
	locLstWrapped = mod.( locLst .- 1, Nparam ) .+ 1;
	return locLstWrapped;
end

function scatter_monopole_no_wrap( Npos , dPosNeg, Nparam )
	posLocLst = rand( 1:Nparam, Npos );
	negLocDiffLst =rand( -dPosNeg:dPosNeg );
	negLocLst = posLocLst .+ negLocDiffLst;

	return posLocLst, negLocLst;
end

function scatter_monopole( Npos , dPosNeg, Nparam )
	posLocLst = rand( 1:Nparam, Npos );
	negLocDiffLst =rand( -dPosNeg:dPosNeg );
	negLocLst = mod.( posLocLst .+ negLocDiffLst .- 1, Nparam ) .+ 1;

	return posLocLst, negLocLst;
end

function scatter_monopole_sym_half( Npos, dPosNeg, Nparam )
	NparamHalf = Int64( Nparam/2 );
	NposHalf = Int64( Npos / 2 );
	posLocLstHalf, negLocLstHalf = scatter_monopole( NposHalf, dPosNeg, NparamHalf );
	posLocLst = vcat( posLocLstHalf, posLocLstHalf .+ NparamHalf );
	negLocLst = vcat( negLocLstHalf, negLocLstHalf .+ NparamHalf );

	return posLocLst, negLocLst;
end

function scatter_monopole_sym_half_switch( Npos, dPosNeg, Nparam )
	NparamHalf = Int64( Nparam/2 );
	NposHalf = Int64( Npos /2 );
	# NnegHalf = Npos - NposHalf;
	posLocLstHalf, negLocLstHalf = scatter_monopole_no_wrap( NposHalf, dPosNeg, NparamHalf );
	posLocLst = monopole_wrap( vcat( posLocLstHalf, negLocLstHalf .+ NparamHalf ), Nparam );
	negLocLst = monopole_wrap( vcat( negLocLstHalf, posLocLstHalf .+ NparamHalf ), Nparam );

	return posLocLst, negLocLst;
end

function scatter_monopole_part_sym_half( Npos, dPosNeg, Nparam )
	NposSym = Int64( floor( Npos/2 ) );
	NposNonSym = Npos - NposSym;
	posLocLstSym, negLocLstSym = scatter_monopole_sym_half( NposSym, dPosNeg, Nparam );
	posLocLstNonSym, negLocLstNonSym = scatter_monopole( NposNonSym, dPosNeg, Nparam );
	posLocLst = vcat( posLocLstSym, posLocLstNonSym );
	negLocLst = vcat( negLocLstSym, negLocLstNonSym );
	
	return posLocLst, negLocLst;
end

function scatter_monopole_part_sym_half_switch( Npos, dPosNeg, Nparam )
	NposSym = Int64( floor( Npos/2 ) );
	NposNonSym = Npos - NposSym;
	posLocLstSym, negLocLstSym = scatter_monopole_sym_half_switch( NposSym, dPosNeg, Nparam );
	posLocLstNonSym, negLocLstNonSym = scatter_monopole( NposNonSym, dPosNeg, Nparam );
	posLocLst = vcat( posLocLstSym, posLocLstNonSym );
	negLocLst = vcat( negLocLstSym, negLocLstNonSym );
	
	return posLocLst, negLocLst;
end

function scatter_monopole_free( Npos , dPosNeg, Nparam )
	posLocLst = rand( 1:Nparam, Npos );
	negLocLst = rand( 1:Nparam, Npos );

	return posLocLst, negLocLst;
end

function deltaC_scatter( Npos, Nparam, posLocLst, negLocLst )
	posChargeLst = ones(Float64,Npos);
	negChargeLst = -ones(Float64,Npos);
	
	chargeLst = vcat(posChargeLst, negChargeLst);
	posLst = vcat(posLocLst, negLocLst);
	
	posLstSorted, chargeLstSorted = SortABbyA( posLst, chargeLst );
	
	Clst = zeros( Nparam );
	pos = 1;
	Cnow = 0;
	for id = 1 : length( posLstSorted )
		Npos = posLstSorted[id];
		Clst[pos : Npos-1] .= Cnow;
		Cnow += chargeLstSorted[id];
		pos = Npos;
	end
	Clst[pos : end] .= Cnow;
	# Ccumulated = accumulate( +, chargeLstSorted );
	
	return Clst;
end

function deltaC_scatter_avg( Npos, Nparam, dPosNeg, avgNum, posFun = scatter_monopole, fileNameMod = "" )
	CsqAvg = zeros( Nparam );
	for it = 1 : avgNum
		posLocLst, negLocLst = posFun( Npos, dPosNeg, Nparam );
		Clst = deltaC_scatter( Npos, Nparam, posLocLst, negLocLst );
		CsqAvg += Clst.^2;
	end
	CsqAvg = CsqAvg ./ avgNum;
	oFileName = string( "Clst_RndWalk_Npos_", Npos, "_Nparam_", Nparam, "_dPosNeg_", dPosNeg, "_avgNum_", avgNum, fileNameMod, npyType );
	npzwrite( oFileName, CsqAvg );
	
	return CsqAvg;
end
