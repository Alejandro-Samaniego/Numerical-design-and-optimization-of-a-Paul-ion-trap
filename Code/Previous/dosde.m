%% 2D CIRCUMFERENCE (line charge 3D)
format long
d=2;
R=1;
N=100; %Número de puntos de carga en UNA bola en R2.
h=2*pi*R/N*ones(2*N,1);
%epsilon0=8.8541878176e-12;
epsilon0=1;
Cteorico=pi*epsilon0/acosh(d/R);

V0=4;
b=V0/2*ones(2*N,1); % Término independiente.
b(N+1:end)=-b(N+1:end);
rni=R*exp(1j*linspace(0,2*pi,N))-d; %Puntos de las circumeferencia de la izquierda
rnd=R*exp(1j*linspace(0,2*pi,N))+d; %Puntos de las circumeferencia de la derecha
rn=[rni rnd]';
Z=zeros(2*N);
for kk=1:2*N
    Z(:,kk)=-h(kk)/(2*pi*epsilon0)*log(abs(rn-rn(kk)));
    Z(kk,kk)=-h(kk)/(2*pi*epsilon0)*(log(h(kk)/2)-1);
end

q=Z\b;
Cexperimental=h(1:N)'*q(1:N)/V0;

m=20;%m de meshgrid
x = 0.2+linspace(-2*R-d,2*R+d,m);
y = 0.2+linspace(-2*R-d,2*R+d,m);
[xx,yy] = meshgrid(x,y); 
rr = xx + 1j*yy;
V=zeros(m,m);
for ii=1:2*N
   RR=abs(rr-rn(ii)); % RR la distancia entre los puntos de carga y los del meshgrid.
   V=V+1/(4*pi*epsilon0)*q(ii)*h(ii)./RR; %Coulomb
end
[Ex, Ey] = gradient(-V);
figure(1)
surf(x,y,V);
shading interp
colormap jet
colorbar
%axis([-d-R d+R -d-R d+R -3 3])

figure(2)
pcolor(x,y,V); 
shading interp
colormap jet
colorbar
hold on; 
quiver(xx,yy,Ex,Ey); 
contour(xx,yy,V,'color','k');


return

%% 2D PLANAR CAPACITOR

format long
d=1;
L=10;
N=200; %Número de puntos de carga en UNA placa.
S=(L)^2;
%epsilon0=8.8541878176*10^(-12);
epsilon0=1;
Cteorico=S*epsilon0/d;
h=L/N*ones(2*N,1);

V0=4;
b=V0/2*ones(2*N,1); % Término independiente.
b(N/2+1:end)=-b(N/2+1:end);

rni=-d/2*ones(1,N)+1j*linspace(0,L,N); %Puntos de la placa de la izquierda
rnd=d/2*ones(1,N)+1j*linspace(0,L,N);%Puntos de la placa de la derecha
rn=[rni rnd]';
Z=zeros(2*N);
for kk=1:2*N
Z(:,kk)=-h(kk)/(2*pi*epsilon0)*log(abs(rn-rn(kk)));
end
Z(isinf(Z))=0;
Zd=diag(-h/(2*pi*epsilon0).*(log(h/2)-1));
Z=Z+Zd; %MATRIZ DEL SISTEMA A SOLUCIONAR

q=Z\b;
Cexperimental=sum(h.*q)/V0;