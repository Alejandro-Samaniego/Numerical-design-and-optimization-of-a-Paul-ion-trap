function posmat=ionprimeroAC(p,q,m,deltat,N,V,X,Y,Z,w)

x=p(1); y=p(2); z=p(3); vx=p(4); vy=p(5); vz=p(6);
%primeros 4 steps
[exs, eys, ezs] = gradient(-V); %static
posmat=zeros(3,N);
posmat(:,1)=[x; y; z];
posmat(:,[2:4])=posmat(:,1)*[1;1;1]'+([deltat/2;deltat;deltat*3/2]*[vx;vy;vz]')';
% SIguientes steps
for ii=5:N
  ex=exs*cos(w*ii*deltat);
  ey=eys*cos(w*ii*deltat);
  ez=ezs*cos(w*ii*deltat);
  Ex = interp3(X,Y,Z,ex,posmat(1,ii-1),posmat(2,ii-1),posmat(3,ii-1));
  Ey = interp3(X,Y,Z,ey,posmat(1,ii-1),posmat(2,ii-1),posmat(3,ii-1));
  Ez = interp3(X,Y,Z,ez,posmat(1,ii-1),posmat(2,ii-1),posmat(3,ii-1));
  acc=[q*Ex/m;q*Ey/m;q*Ez/m];
  posmat(:,ii)=4*deltat^2*acc+2*posmat(:,ii-2)-posmat(:,ii-4);      
end

end