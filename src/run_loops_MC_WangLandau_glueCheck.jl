using Loops_MC
using Utils
using FilenameManip
using JLD2

using DelimitedFiles

dosArrZonedSampleLst = load( fNameGotten, "dosArrReplicaSample" );

nDim = 3;
D_hist = 2;
divNum = 6;

numZones = 8;
EMinRatio = -2.0;
EMaxRatio = 2.0;
EOverlapRatio = 0.7;

histDosType = Loops_MC.WLHistDosZonedInE{nDim,D_hist};

idMinLst, idMaxLst = Loops_MC.genWLZone_idMinMaxLst( histDosType, divNum, EMinRatio, EMaxRatio, EOverlapRatio, numZones );

histDosLst = [ histDosType( divNum, idMinLst[ii], idMaxLst[ii] ) for ii = 1 : numZones ];

for ii = 1 : numZones
	histDosLst[ii].dosArr .= @view( dosArrZonedSampleLst[ii,1][:,:,end] );
end

iMatchLst, dosShLst = Loops_MC.findGluePtsDosArrReplica( histDosLst );

# dosArrFull = glueDosArrReplicaFromGluePts(  );

# dosArr2dLst = Loops_MC.getDosArr.( histDosLst );
# dosMinLst = minimum.(a -> a !=0 ? a : Inf, dosArr2dLst);
# for iHist = 1 : length( dosArr2dLst )
	# for i2d = 1 : length( dosArr2dLst[iHist] )
		# if dosArr2dLst[iHist][i2d] != 0
			# dosArr2dLst[iHist][i2d] -= dosMinLst[iHist];
			# dosArr2dLst[iHist][i2d] += log(2);
		# end
	# end
# end

# dosMaxLst = maximum.(dosArr2dLst);

# dosArrExpShMaxLst = [ similar( dosArr2dLst[ii] ) for ii = 1 : length(histDosLst) ];
# for iHist = 1 : length( dosArr2dLst )
	# for i2d = 1 : length( dosArr2dLst[iHist] )
		# if dosArr2dLst[iHist][i2d] == 0
			# dosArrExpShMaxLst[iHist][i2d] = 0;
		# else
			# dosArrExpShMaxLst[iHist][i2d] = exp(dosArr2dLst[iHist][i2d] - dosMaxLst[iHist]);
		# end
	# end
# end

# # dosArr1dLst = ( d -> dropdims( log.( sum( exp.( d ); dims = 1 ) ); dims = 1 ) ).(dosArr2dLst);
# dosArr1dLst = ( (d,mx) -> dropdims( log.( sum( d; dims = 1 ) ); dims = 1 ) .+ mx ).(dosArrExpShMaxLst, dosMaxLst);

# iMatchLst, dosShLst = Loops_MC.findGluePtsDosArrReplica( dosArr1dLst, idMinLst, idMaxLst );

dosArrFull = Loops_MC.glueDosArrReplicaFromGluePts( histDosLst, iMatchLst, dosShLst );

fAllAttr = FilenameManip.extractFNameAttrs( fNameGotten );

fMainDosArrFull = "dosArrFullTest";

fNameDosArrFull = fNameFunc( fMainDosArrFull, [], [], jld2Type; fMod = fAllAttr );

save( fNameDosArrFull, "dosArrFull", dosArrFull );

open(Loops_MC.fNameFileLstWL, "w") do io
	println(io, fNameDosArrFull);
end
