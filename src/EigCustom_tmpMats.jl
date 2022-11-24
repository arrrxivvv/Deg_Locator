struct EigTmpMats
	Hmat::Matrix{ComplexF64};
	vLst::Matrix{ComplexF64};
	Elst::Vector{Float64};
	
	eigWork::EigWork;
end

function eigTmpMatsInit( mSz::Int64 )
	Hmat = zeros(ComplexF64,mSz,mSz);
	vLst = zeros(ComplexF64,mSz,mSz);
	Elst = zeros(Float64,mSz);
	
	eigWork = eigenPreworkStruct!(Hmat, Elst, vLst);
	
	return EigTmpMats( Hmat, vLst, Elst, eigWork );
end

function eigOnTmpMats!( eigTmpMats::EigTmpMats )
	eigenZheevrStruct!( eigTmpMats.Hmat, eigTmpMats.Elst, eigTmpMats.vLst, eigTmpMats.eigWork );
end
