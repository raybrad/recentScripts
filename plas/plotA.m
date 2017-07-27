function plotA(axisxyz,m)
%clear;
fname=strcat('dump/variables_',num2str(m),'.mat');
load(fname);
epsilon = 8.854e-12;
q = 1.6021892e-19;
lambda=1e-9;
dt=1e-18;
Vt=2.5852e-2;
s_D = 1;
tao = lambda^2/s_D;
s_A = tao*Vt/lambda;
s_B = tao*Vt/lambda^2;
s_E = Vt/lambda;
ni = epsilon*Vt/q/lambda^2;
s_J = q*ni*s_D/lambda;
 
no_of_nodes_x=16;
no_of_nodes_y=16;
no_of_nodes_z=16;


VecPotx(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
VecPoty(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
VecPotz(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;

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
VecPotx=VecPotx*s_A;
VecPoty=VecPoty*s_A;
VecPotz=VecPotz*s_A;

[X,Y]=ndgrid(1.0:1.0:16.0,1.0:1.0:16.0);
%Vecx(:,:)=VecPotx(:,:,8);Vecy(:,:)=VecPoty(:,:,8);
%figure;
%quiver(X,Y,Vecx,Vecy);
switch axisxyz
case 1
aveA(:,:)=sqrt(VecPotx(8,:,:).^2+VecPoty(8,:,:).^2+VecPotz(8,:,:).^2);
AA1(:,:)=VecPoty(8,:,:);
AA2(:,:)=VecPotz(8,:,:);
case 2
aveA(:,:)=sqrt(VecPotx(:,8,:).^2+VecPoty(:,8,:).^2+VecPotz(:,8,:).^2);
AA1(:,:)=VecPotx(:,8,:);
AA2(:,:)=VecPotz(:,8,:);
case 3
aveA(:,:)=sqrt(VecPotx(:,:,8).^2+VecPoty(:,:,8).^2+VecPotz(:,:,8).^2);
AA1(:,:)=VecPotx(:,:,8);
AA2(:,:)=VecPoty(:,:,8);
otherwise
display('error');
end

figure;
contourf(X,Y,aveA,10);
%surf(X,Y,aveA);

figure;
quiver(X,Y,AA1,AA2);

%%quiver 3D vector arrows
figure;
[X,Y,Z]=ndgrid(1.0:2.0:16.0,1.0:2.0:16.0,1.0:2.0:16.0);
Axx(:,:,:)=VecPotx(1:2:16,1:2:16,1:2:16);
Ayy(:,:,:)=VecPoty(1:2:16,1:2:16,1:2:16);
Azz(:,:,:)=VecPotz(1:2:16,1:2:16,1:2:16);
quiver3(X,Y,Z,Axx,Ayy,Azz);

clear;
