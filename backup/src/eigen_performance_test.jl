using LinearAlgebra
using Random

include("utils.jl")
include("Umat.jl")

function testEigen()
	num_it = 50^3;
	N = 20;
	for it = 1:num_it
		# Hmat = H_GUE( N );
		HmatRaw = rand( Complex{Float64}, N, N );
		Hmat = HmatRaw + adjoint(HmatRaw);
		eigenObj = eigen(Hmat);
	end
end
