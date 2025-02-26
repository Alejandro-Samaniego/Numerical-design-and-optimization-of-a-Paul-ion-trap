% Vamos a diseñar una tramopa de unos 10 centimetros
delta=[0.01,2*pi/30];
R=0.002; it=7; v=[0.04,sqrt(2)*0.04,-0.03,0.03]; 
%ponemos altura 0.04 y radio sqrt(2)*altura para mas info bibliografia Paul
%trap aix marseille pag 9
[~,~,~,obj]=hiperbolic(delta,R,it,v);
epsilon0=8.8541878176*10^(-12);
N=10;
Lx=0.05; Ly=0.05; Lz=0; d=0.001;
v1 = obj.Posmat(:,obj.topol(1,:));
v2 = obj.Posmat(:,obj.topol(2,:));
v3 = obj.Posmat(:,obj.topol(3,:));
obj.cent =(v1+v2+v3)/3;
c = cross(v2-v1, v3-v1);
obj.ds = sqrt(sum(c.^2))/2;
obj.un = c./repmat(2*obj.ds,3,1); 
Z=zeros(max(size(obj.cent)));
for ii=1:max(size(obj.cent))
Z(ii,:) = int_S_1divR(obj.cent(:,ii) , v1 , v2 , v3 , obj.un , obj.cent)/(4*pi*epsilon0);
end
V0=0.002;
b=V0/2*ones(max(size(obj.cent)),1);
b(end-max(size(obj.topolring))+1:end)=-b(end-max(size(obj.topolring))+1:end);
Nfb1p=numel(obj.ds);
q=Z\b;
% Computacion de V y E
x = linspace(-Lx,Lx,10);
y = linspace(-Ly,Ly,10); 
z = 0.001+linspace(Lz-d-Lx,Lz+Lx,10);
[xx,yy,zz] = meshgrid(x,y,z);
V=zeros(10*10*10,1);
for ii=1:max(size(obj.cent))
   RR=[xx(:),yy(:),zz(:)]'-obj.cent(:,ii); % RR la distancia entre los 
   %puntos de carga y los del meshgrid.
   RR=vecnorm(RR);
   V=V+1/(4*pi*epsilon0)*q(ii)*obj.ds(ii)./RR'; %Coulomb
end
Vmat=zeros(10,10,10);
for jj=1:10
   for kk=1:10
       Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
   end
end
% 1 ION de Berilio+ MOVIMIENTO MARIPOSA
n=1;
p=rand(6,n)-rand(6,n); 
qion=1.602*10^(-19)*ones(1,n);
p([1:3],:)=0.025*p([1:3],:); 
m=9.012182*1.6605*10^(-27)*ones(1,n); 
deltat=2*10^(-4); 
N=500; 
w=1000;
posmat=ionprimeroAC(p,qion,m,deltat,N,Vmat,xx,yy,zz,w);