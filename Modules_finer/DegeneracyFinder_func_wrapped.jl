using LinearAlgebra

function Hamiltonian(x,y,z)
	(sin(x)+cos(x))*[0 1; 1 0] + (sin(y)+cos(y))*[0 -1im;1im 0] + (sin(z)+cos(z))*[1 0;0 -1];
end

function FieldCompute(Vec0,VecX,VecY,VecXY)

	tmp = dot(Vec0,VecX);
	tmp1 = dot(VecX,VecXY);
	tmp2 = dot(VecXY,VecY);
	tmp3 = dot(VecY,Vec0);
                                               
	ntmp = tmp/norm(tmp);
	ntmp1 = tmp1/norm(tmp1);
	ntmp2 = tmp2/norm(tmp2);
	ntmp3 = tmp3/norm(tmp3);

	log(ntmp*ntmp1*ntmp2*ntmp3)/1im;
end

function PointCompute(x,y,z)
	dx = 0.001;
	dy = 0.001;
	dz = 0.001;

	Do,Vo = eigen(Hamiltonian(x,y,z));
	Dx,Vx = eigen(Hamiltonian(x+dx,y,z));
	Dy,Vy = eigen(Hamiltonian(x,y+dy,z));
	Dz,Vz = eigen(Hamiltonian(x,y,z+dz));
	Dxy,Vxy = eigen(Hamiltonian(x+dx,y+dy,z));
	Dxz,Vxz = eigen(Hamiltonian(x+dx,y,z+dz));
	Dyz,Vyz = eigen(Hamiltonian(x,y+dy,z+dz));



	nlVo = Vo[:,1]/(norm(Vo[:,1]));
	nlVx = Vx[:,1]/(norm(Vx[:,1]));
	nlVy = Vy[:,1]/(norm(Vy[:,1]));
	nlVz = Vz[:,1]/(norm(Vz[:,1]));
	nlVxy = Vxy[:,1]/(norm(Vxy[:,1]));
	nlVxz = Vxz[:,1]/(norm(Vxz[:,1]));
	nlVyz = Vyz[:,1]/(norm(Vyz[:,1]));

	Bx = FieldCompute(nlVo,nlVy,nlVz,nlVyz);
	By = FieldCompute(nlVo,nlVz,nlVx,nlVxz);
	Bz = FieldCompute(nlVo,nlVx,nlVy,nlVxy);

	real([Bx By Bz]);
end

using Distributions

function main()

	ini = rand(Uniform(0,pi),20,3);
	depoint = zeros(20,3);


	for j = 1:20
		tmp = ini[j,:];
		for i = 1:1000
			dtmp = PointCompute(tmp[1],tmp[2],tmp[3]);
			println(norm(dtmp));
			scale = 0.001/norm(dtmp);
			ntmp = zeros(3);
			# ntmp[1] = mod(tmp[1] - scale*dtmp[1],2pi);
			# ntmp[2] = mod(tmp[2] - scale*dtmp[2],2pi);
			# ntmp[3] = mod(tmp[3] - scale*dtmp[3],2pi);
			ntmp[1] = mod(tmp[1] + scale*dtmp[1],2pi);
			ntmp[2] = mod(tmp[2] + scale*dtmp[2],2pi);
			ntmp[3] = mod(tmp[3] + scale*dtmp[3],2pi);
			tmp = ntmp;
		end

		if j == 1
			depoint[j,:] = round.(tmp;digits=2);
		else
			# for i

			# end
		end
	end
	depoint
end
