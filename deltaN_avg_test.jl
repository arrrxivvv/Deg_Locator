include("deltaN_avg.jl")

function deltaN_avg_test()
	avgNumLst = [500,1000,1500];
	
	for avgNum in avgNumLst
		deltaN_avg_lst( avgNum );
	end
end
