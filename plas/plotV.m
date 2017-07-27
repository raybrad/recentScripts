function plotV(axisxyz,m)
fname=strcat('variables_',num2str(m),'.mat');
load(fname);
 
Vt=2.5852e-2;

no_of_nodes_x=16;
no_of_nodes_y=16;
no_of_nodes_z=16;

% get V 
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

Voltage=Voltage*Vt;
%
[X,Y]=ndgrid(1.0:1.0:16.0,1.0:1.0:16.0);

switch axisxyz
case 1
V_plane(:,:)=Voltage(8,:,:);
case 2 
V_plane(:,:)=Voltage(:,8,:);
case 3 
V_plane(:,:)=Voltage(:,:,8);
otherwise
display('error');
end

figure;

contourf(X,Y,V_plane,10);
%surf(X,Y,dtA);

%%quiver 3D vector arrows

%[X,Y,Z]=meshgrid(1.0:1.0:10.0,1.0:1.0:16.0,1.0:1.0:10.0);
%quiver3(X,Y,Z,DiVecPotx,DiVecPoty,DiVecPotz);

clear;
