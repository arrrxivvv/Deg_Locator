push!(LOAD_PATH,"./");

using Utils
using DegLocatorDiv

num_it = 1000;
sym_fact = 4;
num_dim = 3;
it_sym = 4^num_dim;
num_it_full = num_it * it_sym;

Nlst = [10];
param_divide = 26;
seed = -1;

@time divB_profile( Nlst, num_it, param_divide, seed );
print("instance number: ", num_it);
deltaN_avg_lst_fromFile( Nlst, num_it_full, param_divide, seed );
