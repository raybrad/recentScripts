function checkA
clear;

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

fpx = fopen('Axline.dat','w');
fpy = fopen('Ayline.dat','w');
fpz = fopen('Azline.dat','w');

for m=1:150
fname=strcat('variables_',num2str(m),'.mat');
load(fname);

%read dA/dt
Ax(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
Ay(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;
Az(1:no_of_nodes_x,1:no_of_nodes_y,1:no_of_nodes_z) =0.0;

tmpcounter = 0;
for k = 1:no_of_nodes_z
    for j = 1:no_of_nodes_y
        for i = 1:no_of_nodes_x
                if (i < no_of_nodes_x) 
	        tmpcounter = 1 + tmpcounter;
		Ax(i,j,k)=A(tmpcounter);
		else
		Ax(i,j,k)=Ax(i-1,j,k);
		end
		   if (i>1 && i<no_of_nodes_x)
		   Ax(i,j,k)=0.5*(Ax(i,j,k)+Ax(i-1,j,k));
		   end

                if (j < no_of_nodes_y) 
                tmpcounter = 1 + tmpcounter;
		Ay(i,j,k)=A(tmpcounter);
		else
		Ay(i,j,k)=Ay(i,j-1,k);
		end
		   if (j>1 && j<no_of_nodes_y)
		   Ay(i,j,k)=0.5*(Ay(i,j,k)+Ay(i,j-1,k));
		   end

                if (k < no_of_nodes_z) 
                tmpcounter = 1 + tmpcounter;
		Az(i,j,k)=A(tmpcounter);
		else
		Az(i,j,k)=Az(i,j,k-1);
		end
		   if (k>1 && k<no_of_nodes_z)
		   Az(i,j,k)=0.5*(Az(i,j,k)+Az(i,j,k-1));
		   end
        end
    end
 end
Ax=Ax*s_A;
Ay=Ay*s_A;
Az=Az*s_A;


fprintf(fpx,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15,Ay(1,8,6),Ay(3,8,6),Ay(5,8,6),Ay(7,8,6),Ay(9,8,6),Ay(11,8,6),Ay(13,8,6),Ay(16,8,6));

fprintf(fpy,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15,Ay(8,1,6),Ay(8,3,6),Ay(8,5,6),Ay(8,7,6),Ay(8,9,6),Ay(8,11,6),Ay(8,13,6),Ay(8,16,6));

fprintf(fpz,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n',dt*m/1e-15,Ay(8,8,1),Ay(8,8,3),Ay(8,8,5),Ay(8,8,7),Ay(8,8,9),Ay(8,8,11),Ay(8,8,13),Ay(8,8,16));
end

fclose(fpx);
fclose(fpy);
fclose(fpz);

clear;
