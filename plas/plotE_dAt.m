function plotE_dAt(axisxyz,m)
%clear;
fname=strcat('variables_',num2str(m),'.mat');
load(fname);
 
mu = 4*pi*1e-7;
epsilon = 8.854e-12;
q = 1.6021892e-19;
lambda=1e-9;
lunit=1e-9/lambda;
Vt=2.5852e-2;
dt=1e-18;
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
DiVecPotx(:,:,:)=-DiVecPotx(:,:,:)*s_E;
DiVecPoty(:,:,:)=-DiVecPoty(:,:,:)*s_E;
DiVecPotz(:,:,:)=-DiVecPotz(:,:,:)*s_E;

[X,Y]=ndgrid(1.0:1.0:16.0,1.0:1.0:16.0);
switch axisxyz
case 1
dtA(:,:)=sqrt(DiVecPotx(8,:,:).^2+DiVecPoty(8,:,:).^2+DiVecPotz(8,:,:).^2);
EE1(:,:)=DiVecPoty(8,:,:);
EE2(:,:)=DiVecPotz(8,:,:);
case 2
dtA(:,:)=sqrt(DiVecPotx(:,8,:).^2+DiVecPoty(:,8,:).^2+DiVecPotz(:,8,:).^2);
EE1(:,:)=DiVecPotx(:,8,:);
EE2(:,:)=DiVecPotz(:,8,:);
case 3
dtA(:,:)=sqrt(DiVecPotx(:,:,8).^2+DiVecPoty(:,:,8).^2+DiVecPotz(:,:,8).^2);
EE1(:,:)=DiVecPotx(:,:,8);
EE2(:,:)=DiVecPoty(:,:,8);
otherwise
display('error');
end

figure;
contourf(X,Y,dtA,10);
%surf(X,Y,dtA);

figure;
quiver(X,Y,EE1,EE2);

%%quiver 3D vector arrows
figure;
[X,Y,Z]=ndgrid(1.0:2.0:16.0,1.0:2.0:16.0,1.0:2.0:16.0);
Exx(:,:,:)=-DiVecPotx(1:2:16,1:2:16,1:2:16);
Eyy(:,:,:)=-DiVecPoty(1:2:16,1:2:16,1:2:16);
Ezz(:,:,:)=-DiVecPotz(1:2:16,1:2:16,1:2:16);
quiver3(X,Y,Z,Exx,Eyy,Ezz);

clear;
