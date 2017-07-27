function checkE
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

fpx = fopen('E2xline.dat','w');
fpy = fopen('E2yline.dat','w');
fpz = fopen('E2zline.dat','w');

for m=1:350
fname=strcat('variables_',num2str(m),'.mat');
load(fname);

%read dA/dt
DiVecPotx(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
DiVecPoty(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
DiVecPotz(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;

tmpcounter = 0;
for k = 1:no_of_nodes_z
    for j = 1:no_of_nodes_y
        for i = 1:no_of_nodes_x
                if (i < no_of_nodes_x) 
	        tmpcounter = 1 + tmpcounter;
		DiVecPotx(i,j,k)=H(tmpcounter);
		else
		DiVecPotx(i,j,k)=DiVecPotx(i-1,j,k);
		end
		   if (i>1 && i<no_of_nodes_x)
		   DiVecPotx(i,j,k)=0.5*(DiVecPotx(i,j,k)+DiVecPotx(i-1,j,k));
		   end

                if (j < no_of_nodes_y) 
                tmpcounter = 1 + tmpcounter;
		DiVecPoty(i,j,k)=H(tmpcounter);
		else
		DiVecPoty(i,j,k)=DiVecPoty(i,j-1,k);
		end
		   if (j>1 && j<no_of_nodes_y)
		   DiVecPoty(i,j,k)=0.5*(DiVecPoty(i,j,k)+DiVecPoty(i,j-1,k));
		   end

                if (k < no_of_nodes_z) 
                tmpcounter = 1 + tmpcounter;
		DiVecPotz(i,j,k)=H(tmpcounter);
		else
		DiVecPotz(i,j,k)=DiVecPotz(i,j,k-1);
		end
		   if (k>1 && k<no_of_nodes_z)
		   DiVecPotz(i,j,k)=0.5*(DiVecPotz(i,j,k)+DiVecPotz(i,j,k-1));
		   end
        end
    end
end
% get V and and div V
tmpcounter = 0;
for k = 1:no_of_nodes_z
    for j = 1:no_of_nodes_y
        for i = 1:no_of_nodes_x
            tmpcounter = (no_of_nodes_x*(j-1))+(no_of_nodes_y*no_of_nodes_x*(k-1))+i;
            Voltage(i,j,k)=V(tmpcounter);
        end
    end
end
% column first , then rows
[dVy,dVx,dVz]=gradient(Voltage,lunity,lunitx,lunitz);

Ex(:,:,:)=-dVx(:,:,:)*s_E-DiVecPotx(:,:,:)*s_E;
Ey(:,:,:)=-dVy(:,:,:)*s_E-DiVecPoty(:,:,:)*s_E;
Ez(:,:,:)=-dVz(:,:,:)*s_E-DiVecPotz(:,:,:)*s_E;
for k = 1:no_of_nodes_z
    for j = 1:no_of_nodes_y
        for i = 1:no_of_nodes_x
		E2(i,j,k)=Ex(i,j,k)^2+Ey(i,j,k)^2+Ez(i,j,k)^2;
        end
    end
end

fprintf(fpx,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15,E2(1,8,8),E2(3,8,8),E2(5,8,8),E2(7,8,8),E2(8,8,8),E2(12,8,8),E2(14,8,8),E2(16,8,8));

fprintf(fpy,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15,E2(8,1,8),E2(8,3,8),E2(8,5,8),E2(8,7,8),E2(8,8,8),E2(8,12,8),E2(8,14,8),E2(8,16,8));

fprintf(fpz,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15,E2(8,8,1),E2(8,8,3),E2(8,8,5),E2(8,8,7),E2(8,8,8),E2(8,8,12),E2(8,8,14),E2(8,8,16));
end

fclose(fpx);
fclose(fpy);
fclose(fpz);

clear;
