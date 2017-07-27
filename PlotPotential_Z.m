periodx=257;
periody=129;
periodt = 0.05;
totalt=0.70;
nn=floor(totalt/periodt)+1;
format = '%7.3f\n';
load('Xgrid.dat');
load('Ygrid.dat');
X=reshape(Xgrid,periody,periodx);
Y=reshape(Ygrid,periody,periodx);
hfig=figure('Renderer','zbuffer');
set(gca,'NextPlot','replaceChildren');
set(gca,'zlim',[-6e-1 1e-1]);
F(1:nn) = struct('cdata',[],'colormap',[]);
for ii=1:nn
   time = periodt*(ii-1);
   FileName  = ['potential' strtrim(num2str(time,format)) '.dat'];
   Potential = load(FileName);
   Pot = reshape(Potential,periody,periodx);
  % h=plot(X(1,:), Pot(65,:));
  % axis ( [10.0 32.0 -9.0 0.5] );
  % set(gca, 'FontSize', 14, 'LineWidth', 1.0);
  % xlabel('\bf{X Direction (A)}', 'FontSize', 18);
  % ylabel('\bf{Potential (eV)}', 'FontSize', 18);
  % text(19,-0.5,['Time=' num2str(time) 'fs'], 'FontSize', 18)
  % line([17.9 17.9],[-9 0.5],'Color','r','LineWidth',3);
  % line([24.4 24.4],[-9 0.5],'Color','r','LineWidth',3);
  %F(ii) = getframe(gcf, [0 0 540 410]);
  surfc(X,Y,Pot);
  F=getframe(hfig);
end

%movie2avi(F,'H2_Potential_-7Z.avi','fps',1);
