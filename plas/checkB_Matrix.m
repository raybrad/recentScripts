function checkB_Matrix(m)

fname=strcat('variables_',num2str(m),'.mat');
load(fname);

mu = 4*pi*1e-7;
epsilon = 8.854e-12;
q = 1.6021892e-19;
lambda=1e-9;
lunit=1e-9/lambda;
Vt=2.5852e-2;
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
[py,px,pz]=gradient(VecPotx,lunit,lunit,lunit);
[qy,qx,qz]=gradient(VecPoty,lunit,lunit,lunit);
[ry,rx,rz]=gradient(VecPotz,lunit,lunit,lunit);

Bx = ry-qz;
By = pz-rx;
Bz = qx-py;

Bx=Bx*s_B;
By=By*s_B;
Bz=Bz*s_B;

BB(:,:)=Bx(:,8,:);
outputname=strcat('yplaneBx',num2str(m),'.dat');
fp = fopen(outputname,'w');

for j=1:16
	for i=1:16
	fprintf(fp,'%12.5e ',BB(i,j));
	end
	fprintf(fp,'\n');
end
fclose(fp);

clear;
