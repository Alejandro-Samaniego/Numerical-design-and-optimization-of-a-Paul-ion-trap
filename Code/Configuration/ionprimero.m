function posmat=ionprimero(p,q,m,deltat,N,V,X,Y,Z)

%  v=ionprimero(p,q,m,epsilon)
% %p(1)=x, p(2)=y, p(3)=z, p(4)=vx, p(5)=vy, p(6)=vz
% %q carga, m massa ion
% % 2*epsilon es el lado del cubo en el qual hacemos meshgrid al lado de donde
% % se encuentra el ion con tal de calcular ahi V, E.
% 
% % F=ma=mx''=qE
% % Let x'=v--> 2nd order diff eq 
% % x'=v
% % v'=qE/m
% %V(1), v(2), v(3) son las eqs para cada componente de la velocidad del ion
% V(1)=p(4);
% v(2)=p(5);
% v(3)=p(6);
% 
% x = -epsilon+p(1):epsilon:epsilon+p(1);
% y = -epsilon+p(2):epsilon:epsilon+p(2);
% z = -epsilon+p(3):epsilon:epsilon+p(3);
% [xx,yy,zz] = meshgrid(x,y,z);
% V=zeros(3*3*3,1);
% for ii=1:max(size(obj.cent))
%    RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los 
%    %puntos de carga y los del meshgrid.
%    RR=vecnorm(RR);
%    V=V+1/(4*pi*epsilon0)*q(ii)*obj.ds(ii)./RR'; %Coulomb
% end
% Vmat=zeros(3,3,3);
% for jj=1:10
%    for kk=1:10
%        Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
%    end
% end
% % ip
% [Ex, Ey, Ez]=gradient(-V);
% v(4)=q*Ex/m;
% v(5)=q*Ey/m;
% v(6)=q*Ez/m;


% DIFERENCIAS CENTRADAS
%P posicions y vel iniciales, q carga, m massa, deltat time-step, N num de
%its deseado, V potencial, X,Y,Z puntos donde se ha calculado el potencial.

x=p(1); y=p(2); z=p(3); vx=p(4); vy=p(5); vz=p(6);
%primeros 4 steps
[ex, ey, ez] = gradient(-V);
posmat=zeros(3,N);
posmat(:,1)=[x; y; z];
posmat(:,[2:4])=posmat(:,1)*[1;1;1]'+([deltat/2;deltat;deltat*3/2]*[vx;vy;vz]')';
% SIguientes steps
for ii=5:N
  Ex = interp3(X,Y,Z,ex,posmat(1,ii-1),posmat(2,ii-1),posmat(3,ii-1));
  Ey = interp3(X,Y,Z,ey,posmat(1,ii-1),posmat(2,ii-1),posmat(3,ii-1));
  Ez = interp3(X,Y,Z,ez,posmat(1,ii-1),posmat(2,ii-1),posmat(3,ii-1));
  acc=[q*Ex/m;q*Ey/m;q*Ez/m];
  posmat(:,ii)=4*deltat^2*acc+2*posmat(:,ii-2)-posmat(:,ii-4);      
end

end