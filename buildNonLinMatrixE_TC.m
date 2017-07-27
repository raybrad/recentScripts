function [J,F,I,Q2,D] = buildNonLinMatrixE_TC(X,option)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% J = [GE+JIV, JGET*D(V)+JIT  ], F = [GE*V+I]
%     [-JQV , GT+JGTT*D(T)-JQT]     [GT*T-Q]

global scl;
global Nnode eqnNodes;
% global nodeLinks linkVolumes nodeVolumes;
global linkL links Nlink;
global dirNodes isDirNodes zeroNodes;
global devices;
global s_Ids;
global TdirNodes;
global linkMaterialS;
global row0 col0;

Nn = length(eqnNodes);
V = zeros(Nnode,1);
V(eqnNodes) = X(1:Nn);
T = X(Nn+1:end);
TlkAll = mean(T(links(:,1:2)),2);
sgmAll = zeros(Nlink,2); dsgmAll = sgmAll;
[sgmAll(:,1),dsgmAll(:,1)] = getTDsigma(TlkAll,2);
[sgmAll(:,2),dsgmAll(:,2)] = getTDsigma(TlkAll,4);
kpaAll = zeros(Nlink,4); dkpaAll = kpaAll;
[kpaAll(:,1), dkpaAll(:,1)] = getTDkappa(TlkAll,1);
[kpaAll(:,2), dkpaAll(:,2)] = getTDkappa(TlkAll,2);
[kpaAll(:,3), dkpaAll(:,3)] = getTDkappa(TlkAll,3);
[kpaAll(:,4), dkpaAll(:,4)] = getTDkappa(TlkAll,4);



gradV = (V(links(:,2))-V(links(:,1)))./linkL;
valE = sum(sgmAll.*linkMaterialS(:,[2,4]),2).';
valGE = (valE./linkL.');
valGE = [valGE;-valGE;-valGE;valGE]; valGE = valGE(:);
GE = sparse(row0,col0,valGE,Nnode,Nnode); GE(dirNodes,:) = 0;
GE = GE(eqnNodes,eqnNodes);
valQV = -valE.*gradV.';
valQ = -0.5*valQV.*gradV.'.*linkL.';
valQ = [valQ;zeros(1,Nlink);zeros(1,Nlink);valQ];
Qmat = sparse(row0,col0,valQ,Nnode,Nnode);
Q1 = diag(Qmat);

valT = kpaAll.*linkMaterialS;
valT = sum(valT,2).'./linkL.';
valGT = [valT;-valT;-valT;valT];
valGT = valGT(:);
GT = sparse(row0,col0,valGT,Nnode,Nnode);
GT(TdirNodes,:) = 0;
%%% GT entries at TdirNodes should be 1's at the diagonals
% rowone = TdirNodes; colone = TdirNodes; valone = ones(length(TdirNodes),1);
% GTone = sparse(rowone,colone,valone,Nnode,Nnode);
% GT = GT+GTone;
if strcmp(option,'JF')
    valJQV = [valQV;-valQV;valQV;-valQV];
    JQV1 = sparse(row0,col0,valJQV,Nnode,Nnode);
    JQV1(TdirNodes,:) = 0;
    JQV1 = JQV1(:,eqnNodes);
    valGET = dsgmAll.*linkMaterialS(:,[2,4]);
    valGET = -sum(valGET,2).'.*gradV.';
    valJGET = [valGET;-valGET;valGET;-valGET]; % dltV for n1 row and n2 row are different
    valJGET = valJGET(:);
    JGET = sparse(row0,col0,valJGET,Nnode,Nnode); JGET(dirNodes,:) = 0;
    JGET = JGET(eqnNodes,:);
    valQT = 0.5*valGET.*gradV.';
    valJQT = [valQT;valQT;valQT;valQT];
    valJQT = valJQT(:);
    JQT1 = sparse(row0,col0,valJQT,Nnode,Nnode);
    JQT1(TdirNodes,:) = 0;
    
    gradT = (T(links(:,2))-T(links(:,1)))./linkL;
    valGTT = dkpaAll.*linkMaterialS;
    valGTT = -sum(valGTT,2).'.*gradT.';
    valJGTT = [valGTT;valGTT;-valGTT;-valGTT]; % dltV for n1 row and n2 row are different
    valJGTT = valJGTT(:);
    JGTT = sparse(row0,col0,valJGTT,Nnode,Nnode);
    JGTT(TdirNodes,:) = 0;
    JGT = GT+JGTT;
end

%%% fill matrices related to semiconductor devices
Ndevice = size(devices,1);
% JF = sparse(Nnode,Nnode);
I = zeros(Nnode,1); Q2 = zeros(Nnode,1); J = 0; 
row = []; col = []; valJIV = []; valJIT = []; valJQV = []; valJQT = [];
for i = 1:Ndevice   
%     ndS = devices{i}.source(:,1); ndG = devices{i}.gate(:,1); ndD = devices{i}.drain(:,1);
%     wtS = devices{i}.source(:,2); wtG = devices{i}.gate(:,2); wtD = devices{i}.drain(:,2);
%     Vs = sum(V(ndS).*wtS); Vg = sum(V(ndG).*wtG); Vd = sum(V(ndD).*wtD);
%     Vgs = Vg-Vs; Vds = Vd-Vs;
%     Tdevice = mean(T([ndS;ndG]));
%     [Ids,dIdV,dIdT] = Fitting(Vds,Vgs,Tdevice);
%     Ids = Ids*s_Ids; dIdV = dIdV*s_Ids*1; dIdT = dIdT*s_Ids;
%     dIdT = [dIdT,-dIdT,0];
%     I(ndD) = I(ndD) + Ids*wtD; I(ndS) = I(ndS) - Ids*wtS;
%     Q2(ndD) = Q2(ndD)+0.5*Ids*Vds*wtD; Q2(ndS) = Q2(ndS)+0.5*Ids*Vds*wtS;
%     if strcmp(option,'F'), continue; end % compute F only
%     NndS = length(ndS); NndG = length(ndG); NndD = length(ndD); % # of nodes of each terminal
%     Nnd = NndS+NndG+NndD;
%     row = [row;ndD*ones(Nnd,1);ndS*ones(Nnd,1)];
%     col = [col;[ndD;ndS;ndG];[ndD;ndS;ndG]];
%     valJIV = [valJIV;dIdV.';-dIdV.']; % dI/dV
%     valJIT = [valJIT;dIdT.';-dIdT.']; % dI/dT
%     dQdV = 0.5*(dIdV.'*Vds+Ids*[1;-1;0]);
%     valJQV = [valJQV;dQdV;dQdV];
%     dQdT = 0.5*dIdT.'*Vds;
%     valJQT = [valJQT;dQdT;dQdT];
    
    %%%% single-node terminal %%%%%%%%%%%%
    ndS = devices(i,1); ndG = devices(i,2); ndD = devices(i,3);
    Vgs = V(ndG)-V(ndS); Vds = V(ndD)-V(ndS);
    Tdevice = (T(ndS)+T(ndD))/2;
    [Ids,dIdV,dIdT] = Fitting(Vds,Vgs,Tdevice);
    Ids = Ids*s_Ids; dIdV = dIdV*s_Ids*1; dIdT = dIdT*s_Ids;
    dIdT = [dIdT,-dIdT,0];
    I(ndD) = I(ndD) + Ids; I(ndS) = I(ndS) - Ids;
    Q2(ndD) = Q2(ndD)+0.5*Ids*Vds; Q2(ndS) = Q2(ndS)+0.5*Ids*Vds;
    if strcmp(option,'F'), continue; end % compute F only
    row = [row;ndD*ones(3,1);ndS*ones(3,1)];
    col = [col;[ndD;ndS;ndG];[ndD;ndS;ndG]];
    valJIV = [valJIV;dIdV.';-dIdV.']; % dI/dV
    valJIT = [valJIT;dIdT.';-dIdT.']; % dI/dT
    dQdV = 0.5*(dIdV.'*Vds+Ids*[1;-1;0]);
    valJQV = [valJQV;dQdV;dQdV];
    dQdT = 0.5*dIdT.'*Vds;
    valJQT = [valJQT;dQdT;dQdT];
end
I = I(eqnNodes);
if strcmp(option,'JF')
    JIV = sparse(row,col,valJIV,Nnode,Nnode);
    JIV = JIV(eqnNodes,eqnNodes);
    JIT = sparse(row,col,valJIT,Nnode,Nnode);
    JIT = JIT(eqnNodes,:);
    
    JQV2 = sparse(row,col,valJQV,Nnode,Nnode); JQV2 = JQV2(:,eqnNodes);
    JQT2 = sparse(row,col,valJQT,Nnode,Nnode); %JQT2 = JQT2(TeqnNodes,TeqnNodes);
    JQV = JQV1+JQV2;
    JQT = JQT1+JQT2;
    J = [GE+JIV,(JGET+JIT);-JQV,JGT-JQT];
end
FE = GE*V(eqnNodes)+I;
Q = Q1+Q2; Q(TdirNodes) = 0;
FT = GT*T-Q;
F = [FE;FT];
1;

% [kpa_sdAll, dkpaT_sdAll] = getTDkappa(TlkAll,1);
% [kpa_mtAll, dkpaT_mtAll] = getTDkappa(TlkAll,2);
% [kpa_inAll,dkpaT_inAll] = getTDkappa(TlkAll,3);
% [kpa_mt1All, dkpaT_mt1All] = getTDkappa(TlkAll,4);
% kpaAll = [kpa_sdAll,kpa_mtAll,kpa_inAll,kpa_mt1All]; dkpaAll = [dkpaT_sdAll,dkpaT_mtAll,dkpaT_inAll,dkpaT_mt1All];
% row = []; col = []; valGT = []; valJGTT = [];
% FT = zeros(Nnode,1);
% allNodes = (1:Nnode);
% for n1 = allNodes      % to solve T
%     if isTDirNodes(n1)
%         row = [row;n1]; col = [col;n1]; valGT = [valGT;1];
%         continue;
%     end
%     ajlk_n1 = nodeLinks{n1}(1,:); % adjacent links of n1
%     ajnd_n1 = nodeLinks{n1}(2,:); % adjacent nodes of n1
%     ajvol_n1 = nodeVolumes{n1}(1,:); % adjacent volumes of n1
%     ajvolM_n1 = volumeM(ajvol_n1,1); % material of adjacent volumes of n1
%     nz = length(ajnd_n1)+1;
%     row_n1 = n1*ones(nz,1); col_n1 = [n1;ajnd_n1']; val_n1JGT = zeros(nz,1);
%     if nnz(diff(ajvolM_n1)) == 0 % all surrounding volumes are of the same material
%         switch ajvolM_n1(1)
%             case 1
%                 kpa_vec = kpa_sdAll(ajlk_n1);dkpaT_vec = dkpaT_sdAll(ajlk_n1);
%             case 2
%                 kpa_vec = kpa_mtAll(ajlk_n1);dkpaT_vec = dkpaT_mtAll(ajlk_n1);
%             case 3
%                 kpa_vec = kpa_inAll(ajlk_n1);dkpaT_vec = dkpaT_inAll(ajlk_n1);
%             case 4
%                 kpa_vec = kpa_mt1All(ajlk_n1);dkpaT_vec = dkpaT_mt1All(ajlk_n1);
%             otherwise
%                 error('undefined material');
%         end
%         L_vec = linkL(ajlk_n1); S_vec = linkS(ajlk_n1);
%         gradT_vec = (T(ajnd_n1)-T(n1))./L_vec;
%         Phi = -(kpa_vec.*gradT_vec); % dPhiT = dPhiT1+dPhiT2 = kpa*d(-gradT)/dT+dkpa/dT*(-gradT)
%         dPhiT1 = kpa_vec./L_vec;
%         dPhiT2 = -dkpaT_vec.*gradT_vec;
%         val_n1JGT(1) = val_n1JGT(1)+sum((dPhiT1+dPhiT2).*S_vec);
%         val_n1JGT(2:end) = val_n1JGT(2:end)+(-dPhiT1+dPhiT2).*S_vec;
%         FT(n1) = FT(n1)+sum(Phi.*S_vec);
%     else % volumen materials may be different. Assemble the matrix by looping over all links
%         L_vec = linkL(ajlk_n1);
%         dPhiT_vec = [sum(kpaAll(ajlk_n1,:).*linkMaterialS(ajlk_n1,:),2),sum(dkpaAll(ajlk_n1,:).*linkMaterialS(ajlk_n1,:),2)];
%         gradT_vec = (T(ajnd_n1)-T(n1))./L_vec;
%         Phi_vec = -dPhiT_vec(:,1).*gradT_vec;
%         FT(n1) = sum(Phi_vec);
%         dPhiT1_vec = dPhiT_vec(:,1)./L_vec;
%         dPhiT2_vec = -dPhiT_vec(:,2).*gradT_vec;
%         val_n1JGT(1) = sum(dPhiT1_vec+dPhiT2_vec);
%         val_n1JGT(2:nz) = -dPhiT1_vec+dPhiT2_vec;
%         1;
%     end
%     row = [row;row_n1]; col = [col;col_n1]; valGT = [valGT;val_n1JGT];
% end
% JGT = sparse(row,col,valGT,Nnode,Nnode);
% FT = FT-Q;
% J = [GE+JIV,(JGET+JIT);-JQV*1,JGT-JQT];
% F = [FE;FT];

1;
% Ids = Ids(eqnNodes);
% F = zeros(Nnode,1);

% Qlk = 0; Qlk1 = 0;
% for i = 1:length(ajlk_n1)
%     n2 = ajnd_n1(i); lk = ajlk_n1(i);
%     if isMetalLinks(lk,1)
%         L = linkL(lk);
%         gradV = (V(n2)-V(n1))/L;
%         sgm = sgmAll(lk);
%         Qlk = Qlk+sgm*gradV^2;
%     end
%     if isMetalLinks(lk,2)
%         L = linkL(lk);
%         gradV = (V(n2)-V(n1))/L;
%         sgm1 = sgm1All(lk);
%         Qlk1 = Qlk1+sgm1*gradV^2;
%     end
% end
% Q(n1) = (metalV(n1,1)*Qlk+metalV(n1,2)*Qlk1)*0.5;

%         dJV_vec0 = zeros(nz-1,2);
%         for i = 1:length(ajlk_n1)
%             n2 = ajnd_n1(i); lk = ajlk_n1(i);
%             ajvol_lk = linkVolumes{lk}(1,:); ajvolS_lk = linkVolumes{lk}(2,:);
%             ajvolM_lk = volumeM(ajvol_lk,2);
%             L = linkL(lk);
%             sgm = sgmAll(lk); dsgm = dsgmAll(lk);
%             metalS = sum(ajvolS_lk(ajvolM_lk == 2));
%             dJV1 = [sgm,dsgm]*(metalS/L);
%             sgm1 = sgm1All(lk); dsgm1 = dsgm1All(lk);
%             metalS1 = sum(ajvolS_lk(ajvolM_lk == 4));
%             dJV2 = [sgm1,dsgm1]*(metalS1/L);
%             dJV = dJV1+dJV2;
%             dJV_vec0(i,:) = dJV;
%             dV = V(n2)-V(n1); gradV = dV/L;
%             val_n1GE(1) = val_n1GE(1)+dJV(1);val_n1GE(i+1) = val_n1GE(i+1)-dJV(1);
%             val_n1JGET(1) = val_n1JGET(1)+dJV(2)*(-dV);val_n1JGET(i+1) = val_n1JGET(i+1)-dJV(2)*(-dV);
%             dQV = -dJV(1)*gradV;
%             val_n1JQV(1) = val_n1JQV(1)+dQV;val_n1JQV(i+1) = val_n1JQV(i+1)-dQV;
%             dQT = (0.5/L)*dJV(2)*dV*dV;
%             val_n1JQT(1) = val_n1JQT(1)+dQT;val_n1JQT(i+1) = val_n1JQT(i+1)+dQT;
%         end
%         dJV1 = [sgm,dsgm]*(metalS/L);

%         for i = 1:length(ajlk_n1)
%             n2 = ajnd_n1(i); lk = ajlk_n1(i);
%             ajvol_lk = linkVolumes{lk}(1,:);
%             ajvolS_lk = linkVolumes{lk}(2,:);
%             ajvolM_lk = volumeM(ajvol_lk,2);
%             L = linkL(lk);
%             gradT = (T(n2)-T(n1))/L;
%             kpa_sd = kpa_sdAll(lk);dkpaT_sd = dkpaT_sdAll(lk);
%             kpa_mt = kpa_mtAll(lk);dkpaT_mt = dkpaT_mtAll(lk);
%             kpa_in = kpa_inAll(lk);dkpaT_in = dkpaT_inAll(lk);
%             kpa_mt1 = kpa_mt1All(lk);dkpaT_mt1 = dkpaT_mt1All(lk);
%             % vectorize calculations for better efficiency
%             kpa_vec = [kpa_sd,kpa_mt,kpa_in,kpa_mt1];
%             dkpaT_vec = [dkpaT_sd,dkpaT_mt,dkpaT_in,dkpaT_mt1];
%             tmp = ajvolS_lk.*kpa_vec(ajvolM_lk);
%             Phi = -sum(tmp)*gradT;
%             dPhiT1 = sum(tmp)/L;
%             dPhiT2 = -sum(ajvolS_lk.*dkpaT_vec(ajvolM_lk))*gradT;
%             val_n1(1) = val_n1(1)+(dPhiT1+dPhiT2);
%             val_n1(i+1) = val_n1(i+1)+(-dPhiT1+dPhiT2);
%             FT(n1) = FT(n1)+Phi;
%         end

% %%% fill matrices not related to semiconductor devices
% row = []; col = []; valGE = []; valJGET = []; valJQV = []; valJQT = [];
% rowT = []; colT = []; valGT = []; valJGTT = [];
% FT = zeros(Nnode,1);
% allNodes = (1:Nnode);
% for n1 = allNodes
%     if isTDirNodes(n1)
%         rowT = [rowT;n1]; colT = [colT;n1]; valGT = [valGT;1];
%         continue;
%     end
%     ajlk_n1 = nodeLinks{n1}(1,:); % adjacent links of n1
%     ajnd_n1 = nodeLinks{n1}(2,:); % adjacent nodes of n1
%     ajvol_n1 = nodeVolumes{n1}(1,:); % adjacent volumes of n1
%     ajvolM_n1 = volumeM(ajvol_n1,1); % material of adjacent volumes of n1
%     nz = length(ajnd_n1)+1;
%     row_n1 = repmat(n1,nz,1); col_n1 = [n1;ajnd_n1']; val_n1GE = zeros(nz,1); val_n1JGET = zeros(nz,1);
%     val_n1JQV = zeros(nz,1); val_n1JQT = zeros(nz,1);val_n1JGT = zeros(nz,1);
%     L_vec = linkL(ajlk_n1); S_vec = linkS(ajlk_n1);
%     dV_vec = V(ajnd_n1)-V(n1);
%     gradT_vec = (T(ajnd_n1)-T(n1))./L_vec;
%     %-----------------------------------------------------------------------------------------------
%     if nnz(diff(ajvolM_n1)) == 0 % all surrounding volumes are of the same material
%         switch ajvolM_n1(1)
%             case 1
%                 kpa_vec = kpa_sdAll(ajlk_n1);dkpaT_vec = dkpaT_sdAll(ajlk_n1);
%             case 2
%                 sgm_vec = sgmAll(ajlk_n1); dsgm_vec = dsgmAll(ajlk_n1);
%                 kpa_vec = kpa_mtAll(ajlk_n1);dkpaT_vec = dkpaT_mtAll(ajlk_n1);
%             case 3
%                 kpa_vec = kpa_inAll(ajlk_n1);dkpaT_vec = dkpaT_inAll(ajlk_n1);
%             case 4
%                 sgm_vec = sgm1All(ajlk_n1); dsgm_vec = dsgm1All(ajlk_n1);
%                 kpa_vec = kpa_mt1All(ajlk_n1);dkpaT_vec = dkpaT_mt1All(ajlk_n1);
%             otherwise
%                 error('undefined material');
%         end
%         %%% fill E matrices %%%%%%%%%%%
%         if isMetalNodes(n1) && ~isDirNodes(n1)
%             dJV_vec = [sgm_vec.*S_vec./L_vec,dsgm_vec.*S_vec./L_vec];
%             dGET_vec = -dJV_vec(:,2).*dV_vec;
%             dQV_vec = -dJV_vec(:,1).*dV_vec./L_vec;
%             dQT_vec = 0.5*dJV_vec(:,2)./L_vec.*(dV_vec).^2;
%             val_n1GE(1) = sum(dJV_vec(:,1)); val_n1GE(2:end) = -dJV_vec(:,1);
%             val_n1JGET(1) = sum(dGET_vec(:,1)); val_n1JGET(2:end) = -dGET_vec;
%             val_n1JQV(1) = sum(dQV_vec); val_n1JQV(2:end) = -dQV_vec;
%             val_n1JQT(1) = sum(dQT_vec);val_n1JQT(2:end) = dQT_vec;
%             row(end+1:end+nz) = row_n1; col(end+1:end+nz) = col_n1; valGE(end+1:end+nz) = val_n1GE;
%             valJGET(end+1:end+nz) = val_n1JGET; valJQV(end+1:end+nz) = val_n1JQV; valJQT(end+1:end+nz) = val_n1JQT;
%         end
%         %%% fill T matrices %%%%%%%%%%%
%         Phi_vec = -(kpa_vec.*gradT_vec).*S_vec; % dPhiT = dPhiT1+dPhiT2 = kpa*d(-gradT)/dT+dkpa/dT*(-gradT)
%         dPhiT1_vec = kpa_vec./L_vec.*S_vec;
%         dPhiT2_vec = -dkpaT_vec.*gradT_vec.*S_vec;
%     else
%         %%% fill E matrices %%%%%%%%%%%%%%%%%%%
%         if isMetalNodes(n1) && ~isDirNodes(n1)
%             dJV_vec = [(sgmAll(ajlk_n1).*linkMaterialS(ajlk_n1,2)+sgm1All(ajlk_n1).*linkMaterialS(ajlk_n1,4))./L_vec,...
%                 (dsgmAll(ajlk_n1).*linkMaterialS(ajlk_n1,2)+dsgm1All(ajlk_n1).*linkMaterialS(ajlk_n1,4))./L_vec];
%             dGET_vec = -dJV_vec(:,2).*dV_vec;
%             dQV_vec = -dJV_vec(:,1).*dV_vec./L_vec;
%             dQT_vec = 0.5*dJV_vec(:,2)./L_vec.*(dV_vec).^2;
%             val_n1GE(1) = sum(dJV_vec(:,1)); val_n1GE(2:end) = -dJV_vec(:,1);
%             val_n1JGET(1) = sum(dGET_vec(:,1)); val_n1JGET(2:end) = -dGET_vec;
%             val_n1JQV(1) = sum(dQV_vec); val_n1JQV(2:end) = -dQV_vec;
%             val_n1JQT(1) = sum(dQT_vec);val_n1JQT(2:end) = dQT_vec;
%             row(end+1:end+nz) = row_n1; col(end+1:end+nz) = col_n1; valGE(end+1:end+nz) = val_n1GE;
%             valJGET(end+1:end+nz) = val_n1JGET; valJQV(end+1:end+nz) = val_n1JQV; valJQT(end+1:end+nz) = val_n1JQT;
%         end
%         %%% fill T matrices %%%%%%%%%%%%%%%
%         dPhiT_vec = [sum(kpaAll(ajlk_n1,:).*linkMaterialS(ajlk_n1,:),2),sum(dkpaAll(ajlk_n1,:).*linkMaterialS(ajlk_n1,:),2)];
%         Phi_vec = -dPhiT_vec(:,1).*gradT_vec;
%         dPhiT1_vec = dPhiT_vec(:,1)./L_vec;
%         dPhiT2_vec = -dPhiT_vec(:,2).*gradT_vec;
%     end
%     val_n1JGT(1) = sum(dPhiT1_vec+dPhiT2_vec);
%     val_n1JGT(2:end) = -dPhiT1_vec+dPhiT2_vec;
% %     rowT = [rowT;row_n1]; colT = [colT;col_n1]; valGT = [valGT;val_n1JGT];
%     rowT(end+1:end+nz) = row_n1; colT(end+1:end+nz) = col_n1; valGT(end+1:end+nz) = val_n1JGT;
%     FT(n1) = sum(Phi_vec);
% end
% GE = sparse(row,col,valGE,Nnode,Nnode); GE = GE(eqnNodes,eqnNodes);
% JGET = sparse(row,col,valJGET,Nnode,Nnode);
% JGET = JGET(eqnNodes,:);
% JQV1 = sparse(row,col,valJQV,Nnode,Nnode); JQV1 = JQV1(:,eqnNodes);
% JQT1 = sparse(row,col,valJQT,Nnode,Nnode); %JQT1 = JQT1(TeqnNodes,TeqnNodes);
% JGT = sparse(rowT,colT,valGT,Nnode,Nnode);
