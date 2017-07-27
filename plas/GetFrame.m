load('AO177.dat');
load('AO178.dat');
load('OccNo.dat');

nOrbs=350;
nTime=61;
AO177R=reshape(AO177(:,1),nOrbs,nTime);
AO177I=reshape(AO177(:,2),nOrbs,nTime);
AO178R=reshape(AO178(:,1),nOrbs,nTime);
AO178I=reshape(AO178(:,2),nOrbs,nTime);
OccDis=reshape(OccNo(:,2),nOrbs,nTime);
XNO=OccNo(1:nOrbs,1);
clear AO177 AO178 OccNo;

%axis tight manual;
hfig=figure;
for ii=1:nTime
  h=plot(XNO, AO177R(:,ii),'r-+', XNO, AO178R(:,ii), 'b-x', XNO, OccDis(:,ii), 'g-');
  axis ( [0.0 350.0 -1.5 1.5] );
  set(gca, 'FontSize', 14, 'LineWidth', 1.0);
  xlabel('\bf{Natural Orbital No.}', 'FontSize', 18);
  ylabel('\bf{Coeffecient}', 'FontSize', 18);
  legend([h(1), h(2), h(3)],'AO H1', 'AO H1', 'Occ No.');
  legend('boxoff');
  F(ii) = getframe(hfig);
end
figure;
%movie(gcf,F,1,1);
movie2avi(F,'test-7.avi','fps',1);
