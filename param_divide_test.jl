include("divB_profile.jl")

function testDivide()
	param_divide_lst = [26, 40, 50];
	Nlst = [10];
	num_it = 5;
	for param_divide in param_divide_lst
		divB_profile( Nlst, num_it, param_divide );
	end
end
