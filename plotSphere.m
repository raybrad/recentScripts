temp=dlmread('spheres.dat');
X=temp(:,1);
Y=temp(:,2);
Z=temp(:,3);
S=temp(:,4);
scatter3sph(X,Y,Z,'size',S);

