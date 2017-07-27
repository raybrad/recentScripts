function sMat=getSparse(infile,nOrbs)

%infile="Fock_AO_orig.mat";
%nOrbs=711;
temp=dlmread(infile);
Idx=temp(:,1);
Jdx=temp(:,2);
Vdx=temp(:,3);
sMat=sparse(Idx,Jdx,Vdx,nOrbs,nOrbs);

