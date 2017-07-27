function plotB(axisxyz,m)
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

Bx = ry-qz;%qz large
By = pz-rx;
Bz = qx-py;%qx ??

Bx=Bx*s_B;
By=By*s_B;
Bz=Bz*s_B;


[X,Y]=ndgrid(1.0:1.0:16.0,1.0:1.0:16.0);
switch axisxyz
case 1
aveB(:,:)=sqrt(Bx(8,:,:).^2+By(8,:,:).^2+Bz(8,:,:).^2);
%aveB(:,:)=Bz(8,:,:).^2;
BB1(:,:)=By(8,:,:);
BB2(:,:)=Bz(8,:,:);
case 2 
aveB(:,:)=sqrt(Bx(:,8,:).^2+By(:,8,:).^2+Bz(:,8,:).^2);
%aveB(:,:)=Bz(:,8,:).^2;
BB1(:,:)=Bx(:,8,:);
BB2(:,:)=Bz(:,8,:);
case 3 
aveB(:,:)=sqrt(Bx(:,:,8).^2+By(:,:,8).^2+Bz(:,:,8).^2);
%aveB(:,:)=Bz(:,:,8).^2;
BB1(:,:)=Bx(:,:,8);
BB2(:,:)=By(:,:,8);
otherwise
display('error');
end

%for i=1:16
%	disp(Bx(8,8,i));
%end

%for i=1:16
%AAA(i)=Bx(8,8,i);
%end
%figure;
%plot(AAA);

figure;
contourf(X,Y,aveB,10);
%surf(X,Y,dtA);

figure;
quiver(X,Y,BB1,BB2);
%%quiver 3D vector arrows
figure;
[X,Y,Z]=ndgrid(1.0:2.0:16.0,1.0:2.0:16.0,1.0:2.0:16.0);
Bxx(:,:,:)=Bx(1:2:16,1:2:16,1:2:16);
Byy(:,:,:)=By(1:2:16,1:2:16,1:2:16);
Bzz(:,:,:)=Bz(1:2:16,1:2:16,1:2:16);
quiver3(X,Y,Z,Bxx,Byy,Bzz);
