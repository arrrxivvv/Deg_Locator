include("divB_profile.jl")

function div_profile_res( Nlst = [20], param_divide_lst = [25,50,75] )
	for param_divide in param_divide_lst
		num_it = 2;
		divB_profile( Nlst, num_it, param_divide );
	end
end
