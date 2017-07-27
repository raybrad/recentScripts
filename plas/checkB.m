function checkB
clear;

mu = 4*pi*1e-7;
epsilon = 8.854e-12;
q = 1.6021892e-19;
lambda=1e-9;
lunitx = [0.0:1.0:15.0];
lunity = [0.0:1.0:15.0];
lunitz = [0.0:1.0:15.0];
lunitx=lunitx*1e-9/lambda;
lunity=lunity*1e-9/lambda;
lunitz=lunitz*1e-9/lambda;
Vt=2.5852e-2;
dt=5e-17;
s_D = 1;
tao = lambda^2/s_D;
s_A = tao*Vt/lambda;
s_B = tao*Vt/lambda^2;
s_E = Vt/lambda;
ni = epsilon*Vt/q/lambda^2;
s_J = q*ni*s_D/lambda;

s_sigma = epsilon/tao; % A/V/m
s_Curr = s_sigma * Vt * lambda; % A
K = epsilon*mu*(lambda/tao)^2; % dimensionless

no_of_nodes_x=16;
no_of_nodes_y=16;
no_of_nodes_z=16;

fpx = fopen('Bxline.dat','w');
fpy = fopen('Byline.dat','w');
fpz = fopen('Bzline.dat','w');

for m=1:150
fname=strcat('variables_',num2str(m),'.mat');
load(fname);

%read dA/dt
VecPotx(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
VecPoty(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
VecPotz(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
Bx(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
By(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
Bz(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;

tmpcounter = 0;
for k = 1:no_of_nodes_z
    for j = 1:no_of_nodes_y
        for i = 1:no_of_nodes_x
                if (i < no_of_nodes_x) 
	        tmpcounter = 1 + tmpcounter;
		VecPotx(i,j,k)=A(tmpcounter);
		else
		VecPotx(i,j,k)=VecPotx(i-1,j,k);
		end
		   if (i>1 && i<no_of_nodes_x)
		   VecPotx(i,j,k)=0.5*(VecPotx(i,j,k)+VecPotx(i-1,j,k));
		   end

           if (j < no_of_nodes_y) 
           tmpcounter = 1 + tmpcounter;
		VecPoty(i,j,k)=A(tmpcounter);
		else
		VecPoty(i,j,k)=VecPoty(i,j-1,k);
		end
		   if (j>1 && j<no_of_nodes_y)
		   VecPoty(i,j,k)=0.5*(VecPoty(i,j,k)+VecPoty(i,j-1,k));
		   end

                if (k < no_of_nodes_z) 
                tmpcounter = 1 + tmpcounter;
		VecPotz(i,j,k)=A(tmpcounter);
		else
		VecPotz(i,j,k)=VecPotz(i,j,k-1);
		end
		   if (k>1 && k<no_of_nodes_z)
		   VecPotz(i,j,k)=0.5*(VecPotz(i,j,k)+VecPotz(i,j,k-1));
		   end
        end
    end
 end
[py,px,pz]=gradient(VecPotx,lunity,lunitx,lunitz);
[qy,qx,qz]=gradient(VecPoty,lunity,lunitx,lunitz);
[ry,rx,rz]=gradient(VecPotz,lunity,lunitx,lunitz);

Bx = ry-qz;
By = pz-rx;
Bz = qx-py;

Bx=Bx*s_B;
By=By*s_B;
Bz=Bz*s_B;

fprintf(fpx,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15,Bx(1,8,6),Bx(3,8,6),Bx(5,8,6),Bx(7,8,6),Bx(9,8,6),Bx(11,8,6),Bx(13,8,6),Bx(16,8,6));

fprintf(fpy,'%12.5e  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15,Bx(8,1,6),Bx(8,3,6),Bx(8,5,6),Bx(8,7,6),Bx(8,9,6),Bx(8,11,6),Bx(8,13,6),Bx(8,16,6));

fprintf(fpz,'%12.5e  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15,Bx(8,8,1),Bx(8,8,3),Bx(8,8,5),Bx(8,8,7),Bx(8,8,9),Bx(8,8,11),Bx(8,8,13),Bx(8,8,16));
end

fclose(fpx);
fclose(fpy);
fclose(fpz);

clear;
