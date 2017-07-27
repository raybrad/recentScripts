function plotJ(axisxyz,m)
%clear;
%fname=strcat('dump/variables_',num2str(m),'.mat');
fname=strcat('topold/variables_',num2str(m),'.mat');
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

[X,Y]=ndgrid(1.0:1.0:16.0,1.0:1.0:16.0);
%Vecx(:,:)=Jx(:,:,8);Vecy(:,:)=Jy(:,:,8);
%figure;
%quiver(X,Y,Vecx,Vecy);
switch axisxyz
case 1
aveJ(:,:)=sqrt(Jx(8,:,:).^2+Jy(8,:,:).^2+Jz(8,:,:).^2);
J1(:,:)=Jy(8,:,:);
J2(:,:)=Jz(8,:,:);
case 2
aveJ(:,:)=sqrt(Jx(:,8,:).^2+Jy(:,8,:).^2+Jz(:,8,:).^2);
J1(:,:)=Jx(:,8,:);
J2(:,:)=Jz(:,8,:);
case 3
aveJ(:,:)=sqrt(Jx(:,:,8).^2+Jy(:,:,8).^2+Jz(:,:,8).^2);
J1(:,:)=Jx(:,:,8);
J2(:,:)=Jy(:,:,8);
otherwise
display('error');
end

figure;
contourf(X,Y,aveJ,10);
%surf(X,Y,aveJ);

figure;
quiver(X,Y,J1,J2);

%%quiver 3D vector arrows
figure;
[X,Y,Z]=ndgrid(1.0:2.0:16.0,1.0:2.0:16.0,1.0:2.0:16.0);
Jxx(:,:,:)=Jx(1:2:16,1:2:16,1:2:16);
Jyy(:,:,:)=Jy(1:2:16,1:2:16,1:2:16);
Jzz(:,:,:)=Jz(1:2:16,1:2:16,1:2:16);
quiver3(X,Y,Z,Jxx,Jyy,Jzz);

clear;
