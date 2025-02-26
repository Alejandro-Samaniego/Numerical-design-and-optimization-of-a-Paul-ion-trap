function [posmat1,posmat2]=dos_iones(p,q,m,deltat,N,V,X,Y,Z)
% DIFERENCIAS CENTRADAS
%P posicions y vel iniciales, q carga, m massa,(tenemos dos cargas y dos masas por eso entramos un vector de cada) deltat time-step, N num de
%its deseado, V potencial, X,Y,Z puntos donde se ha calculado el potencial.
epsilon0=1;
x=p(1); y=p(2); z=p(3); vx=p(4); vy=p(5); vz=p(6);
x2=p(7); y2=p(8); z2=p(9); vx2=p(10); vy2=p(11); vz2=p(12);
%primeros 4 steps
[ex, ey, ez] = gradient(-V);
posmat1=zeros(3,N);
posmat1(:,1)=[x; y; z];
posmat1(:,[2:4])=posmat1(:,1)*[1;1;1]'+([deltat/2;deltat;deltat*3/2]*[vx;vy;vz]')';
posmat2=zeros(3,N);
posmat2(:,1)=[x2; y2; z2];
posmat2(:,[2:4])=posmat2(:,1)*[1;1;1]'+([deltat/2;deltat;deltat*3/2]*[vx2;vy2;vz2]')';
% SIguientes steps
for ii=5:N
  Erepulsion2=q(1)*q(2)/(4*pi*epsilon0)*(posmat2(:,ii-1)-posmat1(:,ii-1))/((norm(posmat2(:,ii-1)-posmat1(:,ii-1)))^3);%Campo de 1 sobre 2
  Erepulsion1=-Erepulsion2;%Campo de 2 sobre 1
  Ex = interp3(X,Y,Z,ex,posmat1(1,ii-1),posmat1(2,ii-1),posmat1(3,ii-1));
  Ey = interp3(X,Y,Z,ey,posmat1(1,ii-1),posmat1(2,ii-1),posmat1(3,ii-1));
  Ez = interp3(X,Y,Z,ez,posmat1(1,ii-1),posmat1(2,ii-1),posmat1(3,ii-1));
  acc=[q(1)*(Ex+Erepulsion1(1))/m(1);q(1)*(Ey+Erepulsion1(2))/m(1);q(1)*(Ez+Erepulsion1(3))/m(1)];
  posmat1(:,ii)=4*deltat^2*acc+2*posmat1(:,ii-2)-posmat1(:,ii-4);   
  Ex2 = interp3(X,Y,Z,ex,posmat2(1,ii-1),posmat2(2,ii-1),posmat2(3,ii-1));
  Ey2 = interp3(X,Y,Z,ey,posmat2(1,ii-1),posmat2(2,ii-1),posmat2(3,ii-1));
  Ez2 = interp3(X,Y,Z,ez,posmat2(1,ii-1),posmat2(2,ii-1),posmat2(3,ii-1));
  acc2=[q(2)*(Ex2+Erepulsion2(1))/m(2);q(2)*(Ey2+Erepulsion2(2))/m(2);q(2)*(Ez2+Erepulsion2(3))/m(2)];
  posmat2(:,ii)=4*deltat^2*acc2+2*posmat2(:,ii-2)-posmat2(:,ii-4); 
end
end