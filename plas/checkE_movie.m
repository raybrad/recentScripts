function checkE_movie(axisxyz)

hfig=figure('Renderer','zbuffer');
%axis tight manual;
set(gca,'NextPlot','replaceChildren');
set(gca,'zlim',[-5e-10 5e-10]);
% Preallocate the struct array for the struct returned by getframe
num=150;
Fb(1:num) = struct('cdata',[],'colormap',[]);

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


for m = 1:1:num
    fname=strcat('variables_',num2str(m),'.mat');
%    display(['load matrix:',fname]);
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
[dVy,dVx,dVz]=gradient(Voltage,lunit,lunit,lunit);

Ex(:,:,:)=-dVx(:,:,:)*s_E-DiVecPotx(:,:,:)*s_E;
Ey(:,:,:)=-dVy(:,:,:)*s_E-DiVecPoty(:,:,:)*s_E;
Ez(:,:,:)=-dVz(:,:,:)*s_E-DiVecPotz(:,:,:)*s_E;

%Ex=Ex/(sqrt(mu/epsilon)/2);
%Ey=Ey/(sqrt(mu/epsilon)/2);
%Ez=Ez/(sqrt(mu/epsilon)/2);
switch axisxyz
case 1
EE(:,:)=Ey(8,:,:);
case 2
EE(:,:)=Ey(:,8,:);
case 3
EE(:,:)=Ey(:,:,8);
otherwise
display('error');
end

surf(X,Y,EE);
Fb(m) = getframe(hfig);

end
%figure;
%movie(gcf,Fb);
%movie2avi(Fb,'E_Matrix_x8yzplane.avi');
clear;
