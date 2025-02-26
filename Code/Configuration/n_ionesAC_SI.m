function posmat=n_ionesAC_SI(p,q,m,deltat,N,V,X,Y,Z,niones,w)
% NO PONEMOS EL GLOBAL

% DIFERENCIAS CENTRADAS
%P posicions y vel iniciales(de todos los iones, estan origanizadas como:columna 1=posiciones iniciales en las tres primeras filas y velocidades en las 3 filas siguientes, columna 2= igual pero con la segunda particula...), q carga, m massa,(tenemos dos cargas y dos masas por eso entramos un vector de cada) deltat time-step, N num de
%its deseado, V potencial, X,Y,Z puntos donde se ha calculado el potencial.
%
%w es la freq del potencial AC
%IMPORTANTE: q y m han de entrar en fila.
epsilon0=8.854187817599999e-12;
% epsilon0=1;
%primeros 4 steps
posmat=zeros(3,N,niones);
posmat(:,1,:)=p(1:3,:); %posicion inicial del los n iones
vecaux=zeros(3,1,niones);
vecaux(:,1,:)=p(4:6,:); %velocidades iniciales de los n iones
posmat(:,[2:4],:)=posmat(:,1,:)+([deltat/2;deltat;deltat*3/2]'.*vecaux);
% asumimos los 4 primeros steps que necesitamos para inicializar como si no
% hubiera fuerza/campo--> PONER UN DELTAT MAS PEQUEÑO????
Erepulsion=zeros(3,niones);
% SIguientes steps
[exs, eys, ezs] = gradient(-V); %static field
for ii=5:N
  for jj=1:niones
      %debemos computar el campo que ve un ion debido a todos los demas
     qjj=zeros(1,1,niones-1);
     qjj(1,1,:)=q([1:jj-1 jj+1:niones]); %todaS LAS CARGAS MENOS LA QUE ESTAMOS TRATANDO
     rjj=posmat(:,ii-1,[1:jj-1 jj+1:niones]);
     %todas las pos menos las del ion que estamos observando: todas filas
     %(x, y, x) columa ii-1 (posiciones anteriores), fondo todas menos la
     %del ion en cuestion
     Erepulsion(:,jj)=sum(qjj/(4*pi*epsilon0).*(posmat(:,ii-1,jj)-rjj)./((vecnorm(posmat(:,ii-1,jj)-rjj)).^3),3);%Campo de todas las cargas sobre carga jj
     %el ultimo 3 indica en que direccion sumamos (para dentro en este caso)
     %CUIDADO: no seria necesario que qjj multiplicara componenete a
     %componente y de hecho qjj no tendria que ser vec en direccion para dentro???
  end
  
  ex=exs*cos(w*ii*deltat);
  ey=eys*cos(w*ii*deltat);
  ez=ezs*cos(w*ii*deltat);
  %[ex, ey, ez] = gradient(-V*cos(w*(ii-5)*deltat));
  %OBS: en caso que funcine lo de poner un delta mas pequeño para los
  %primeros 4 its...
  Ex = interp3(X,Y,Z,ex,reshape(posmat(1,ii-1,:),[1 niones]),reshape(posmat(2,ii-1,:),[1 niones]),reshape(posmat(3,ii-1,:),[1 niones])); %Necesitamos esto en vec fila, para esto es necesario que el vector de xq,yq,zq entre en fila.
  Ey = interp3(X,Y,Z,ey,reshape(posmat(1,ii-1,:),[1 niones]),reshape(posmat(2,ii-1,:),[1 niones]),reshape(posmat(3,ii-1,:),[1 niones])); 
  Ez = interp3(X,Y,Z,ez,reshape(posmat(1,ii-1,:),[1 niones]),reshape(posmat(2,ii-1,:),[1 niones]),reshape(posmat(3,ii-1,:),[1 niones])); 
  acc=[q.*(Ex+Erepulsion(1,:))./m;q.*(Ey+Erepulsion(2,:))./m;q.*(Ez+Erepulsion(3,:))./m];
  acc=reshape(acc,[3 1 niones]);
  posmat(:,ii,:)=4*deltat^2*acc+2*posmat(:,ii-2,:)-posmat(:,ii-4,:);   
end
end