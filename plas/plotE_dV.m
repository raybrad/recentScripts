function plotE_dV(axisxyz,m)
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
[dVy,dVx,dVz]=gradient(Voltage,lunit,lunit,lunit);

Ex(:,:,:)=-dVx(:,:,:)*s_E;
Ey(:,:,:)=-dVy(:,:,:)*s_E;
Ez(:,:,:)=-dVz(:,:,:)*s_E;



%
[X,Y]=ndgrid(1.0:1.0:16.0,1.0:1.0:16.0);

switch axisxyz
case 1
aveE(:,:)=sqrt(Ex(8,:,:).^2+Ey(8,:,:).^2+Ez(8,:,:).^2);
%aveE(:,:)=Ez(8,:,:).^2;
EE1(:,:)=Ey(8,:,:);
EE2(:,:)=Ez(8,:,:);
case 2 
aveE(:,:)=sqrt(Ex(:,8,:).^2+Ey(:,8,:).^2+Ez(:,8,:).^2);
%aveE(:,:)=Ez(:,8,:).^2;
EE1(:,:)=Ex(:,8,:);
EE2(:,:)=Ez(:,8,:);
case 3 
aveE(:,:)=sqrt(Ex(:,:,8).^2+Ey(:,:,8).^2+Ez(:,:,8).^2);
%aveE(:,:)=Ez(:,:,8).^2;
EE1(:,:)=Ex(:,:,8);
EE2(:,:)=Ey(:,:,8);
otherwise
display('error');
end

figure;

contourf(X,Y,aveE,10);
%surf(X,Y,dtA);

%%quiver 3D vector arrows
figure;
quiver(X,Y,EE1,EE2);

%[X,Y,Z]=meshgrid(1.0:1.0:10.0,1.0:1.0:16.0,1.0:1.0:10.0);
%quiver3(X,Y,Z,DiVecPotx,DiVecPoty,DiVecPotz);

figure;
[X,Y,Z]=ndgrid(1.0:2.0:16.0,1.0:2.0:16.0,1.0:2.0:16.0);
Exx(:,:,:)=Ex(1:2:16,1:2:16,1:2:16);
Eyy(:,:,:)=Ey(1:2:16,1:2:16,1:2:16);
Ezz(:,:,:)=Ez(1:2:16,1:2:16,1:2:16);
quiver3(X,Y,Z,Exx,Eyy,Ezz,'color',[0,0,1]);
clear;
