function checkB_movie(axisxyz)

figure('Renderer','zbuffer');
%axis tight manual;
set(gca,'NextPlot','replaceChildren');
set(gca,'zlim',[-3e-18 3e-18]);
% Preallocate the struct array for the struct returned by getframe
num=80
Fa(num+1) = struct('cdata',[],'colormap',[]);
Fb(num+1) = struct('cdata',[],'colormap',[]);

[X,Y]=ndgrid(1.0:1.0:16.0,1.0:1.0:16.0);

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


for m = 1:4:num
    fname=strcat('variables_',num2str(m),'.mat');
    display(['load matrix:',fname]);
    load(fname);
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

switch axisxyz
case 1
BB(:,:)=Bx(8,:,:);
case 2
BB(:,:)=Bx(:,8,:);
case 3
BB(:,:)=Bx(:,:,8);
otherwise
display('error');
end

surf(X,Y,BB);
Fb(m) = getframe;

end
clear;
