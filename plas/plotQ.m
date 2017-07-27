function plotQ(axisxyz)

charge = dlmread('charge.all');
 

Nx=129;
Ny=65;
Nz=65;


tmpcounter = 0;
for k = 1:Nz
    for j = 1:Ny
        for i = 1:Nx
	        tmpcounter = 1 + tmpcounter;
		chargeXYZ(i,j,k)=charge(tmpcounter);
        end
    end
end
%

switch axisxyz
case 1
[Grid1,Grid2]=ndgrid(1.0:1.0:Ny,1.0:1.0:Nz);
plane(:,:)=chargeXYZ(round(Nx/2),:,:);
case 2 
[Grid1,Grid2]=ndgrid(1.0:1.0:Nx,1.0:1.0:Nz);
plane(:,:)=chargeXYZ(:,round(Ny/2),:);
case 3 
[Grid1,Grid2]=ndgrid(1.0:1.0:Nx,1.0:1.0:Ny);
plane(:,:)=chargeXYZ(:,:,round(Nz/2));
otherwise
display('error');
end

figure;

contour(Grid1,Grid2,plane,100);

clear;
