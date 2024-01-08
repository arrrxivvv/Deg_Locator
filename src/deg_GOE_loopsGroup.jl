struct LoopGroup
	ptLst::Vector{Vector{Int64}};
	
	xMin::Int64;
	xMax::Int64;
	yMin::Int64;
	yMax::Int64;
end

function deg_GOE_loopsGroup( mSz, divLst, itNum, seedFed; nDim = 3, fMod = "", fExt = jld2Type, dThres = sqrt(3) )
	fMain = "deg_GOE3";
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seedFed; dim = nDim );
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	locLst = load( fName, "locLst" );
	
	locGrpLst = [ Vector{Vector{Vector{Int64}}}(undef,0) for iM = 1 : mSz ];
	
	for it = 1 : itNum
		for iM = 1 : mSz, dim = 1 : nDim, i3 = 1 : divLst[end]
			numPts = size(locLst[it, i3][iM],2);
			for iPt = 1 : numPts
				if true
					# push!( locGrpLst[iM],  );
				else
				end
			end
		end
	end
end
