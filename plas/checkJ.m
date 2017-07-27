function checkJ
clear;
mu = 4*pi*1e-7;
epsilon = 8.854e-12;
q = 1.6021892e-19;
lambda=1e-9;
Vt=2.5852e-2;
dt=1e-18;
s_D = 1; %m^2/s
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

fpx = fopen('Jxline.dat','w');
fpy = fopen('Jyline.dat','w');
fpz = fopen('Jzline.dat','w');

for m=1:4
fname=strcat('variables_',num2str(m),'.mat');
load(fname);

%read dA/dt
Jx(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
Jy(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
Jz(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;

tmpcounter = 0;
for k = 1:no_of_nodes_z
    for j = 1:no_of_nodes_y
        for i = 1:no_of_nodes_x
                if (i < no_of_nodes_x) 
	        tmpcounter = 1 + tmpcounter;
		Jx(i,j,k)=Js(tmpcounter);
		else
		Jx(i,j,k)=Jx(i-1,j,k);
		end
		   if (i>1 && i<no_of_nodes_x)
		   Jx(i,j,k)=0.5*(Jx(i,j,k)+Jx(i-1,j,k));
		   end

                if (j < no_of_nodes_y) 
                tmpcounter = 1 + tmpcounter;
		Jy(i,j,k)=Js(tmpcounter);
		else
		Jy(i,j,k)=Jy(i,j-1,k);
		end
		   if (j>1 && j<no_of_nodes_y)
		   Jy(i,j,k)=0.5*(Jy(i,j,k)+Jy(i,j-1,k));
		   end

                if (k < no_of_nodes_z) 
                tmpcounter = 1 + tmpcounter;
		Jz(i,j,k)=Js(tmpcounter);
		else
		Jz(i,j,k)=Jz(i,j,k-1);
		end
		   if (k>1 && k<no_of_nodes_z)
		   Jz(i,j,k)=0.5*(Jz(i,j,k)+Jz(i,j,k-1));
		   end
        end
    end
end

Jx=Jx*s_J;
Jy=Jy*s_J;
Jz=Jz*s_J;
fprintf(fpx,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15/1e-15,Jx(1,8,8),Jx(3,8,8),Jx(5,8,8),Jx(7,8,8),Jx(9,8,8),Jx(11,8,8),Jx(15,8,8),Jx(16,8,8));

fprintf(fpy,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15/1e-15,Jy(8,1,8),Jy(8,3,8),Jy(8,5,8),Jy(8,7,8),Jy(8,9,8),Jy(8,11,8),Jy(8,15,8),Jy(8,16,8));

fprintf(fpz,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15/1e-15,Jy(8,8,1),Jy(8,8,3),Jy(8,8,5),Jy(8,8,8),Jy(8,8,9),Jy(8,8,11),Jy(8,8,13),Jy(8,8,16));
end

fclose(fpx);
fclose(fpy);
fclose(fpz);

clear;
