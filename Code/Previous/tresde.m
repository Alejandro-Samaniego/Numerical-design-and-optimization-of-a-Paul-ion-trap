%% CONDENSADOR PLANO-PARALELO TRIDIMENSIONAL
% Para definir los vertices de cada triangulo cojmos una fila y todas las
% columnas de la mat de topologia, donde estan los vertices indexados. Para
% cojer sus cordenadas cojemos una columna de la MATRIZ DE VERTICES 
% (1 COLUMNA, 3 FILAS== X Y Z)
N=10;
Lx=2; Ly=2; Lz=1; d=0.05;
r=[-Lx Lx -Ly Ly Lz];
[obj]=geomtopol(r,N);
v1 = obj.vertex(:,obj.topol(1,:));
v2 = obj.vertex(:,obj.topol(2,:));
v3 = obj.vertex(:,obj.topol(3,:));
obj.cent =(v1+v2+v3)/3;
c = cross(v2-v1, v3-v1);
obj.ds = sqrt(sum(c.^2))/2;%Sum de una matriz suma columnas, luego al hacer 
% esto hacemos norma de cada vector de c (cada vect es una columna).
% /2 porque es triangulo y no paralelogramo
obj.un = c./repmat(2*obj.ds,3,1); % Normalización. Multiplicamos por 2 ya 
% que ahora LO QUE QUEREMOS ES NORMALIZAR LOS VECTORES COLUMA DE C. 
% OBS: obj.un sale una matriz de todo de vectores (0,0,1), tiene sentido ya
% que la normal a un plano paralelo a XY es z...
% OBS: ¿Porque repmat?--> pk c tiene muchos vectores "seguidos" cada uno de
%3 componentes. Cada componente tiene que ir dividida por la norma del 
% vector al cual corresponde, luego hay que dividir 3 escalares por un
% numero, los siguientes 3 por otro etc 

% Plano 2
r2=[-Lx Lx -Ly Ly Lz-d];
[obj2]=geomtopol(r2,N);
v1_2 = obj2.vertex(:,obj2.topol(1,:));
v2_2 = obj2.vertex(:,obj2.topol(2,:));
v3_2 = obj2.vertex(:,obj2.topol(3,:));
obj2.cent =(v1_2+v2_2+v3_2)/3;
c2 = cross(v2_2-v1_2, v3_2-v1_2);
obj2.ds = sqrt(sum(c2.^2))/2;
obj2.un = c2./repmat(2*obj2.ds,3,1);

figure(1)
patch('Faces',obj.topol','Vertices',obj.vertex','FaceColor','g','EdgeColor','k');
hold on
patch('Faces',obj2.topol','Vertices',obj2.vertex','FaceColor','g','EdgeColor','k');
axis equal;
hold off

% Concatenamos la info de los 2 planos
objfinal.cent=[obj.cent obj2.cent];
objfinal.un=[obj.un obj2.un];
objfinal.ds=[obj.ds obj2.ds];
objfinal.vertex=[obj.vertex obj2.vertex];
v1_final=[v1 v1_2];
v2_final=[v2 v2_2];
v3_final=[v3 v3_2];
Z=zeros(max(size(objfinal.cent)));
for ii=1:max(size(objfinal.cent))
Z(ii,:) = int_S_1divR(objfinal.cent(:,ii) , v1_final , v2_final , v3_final , objfinal.un , objfinal.cent)/(4*pi*epsilon0);
end
V0=4;
b=V0/2*ones(max(size(objfinal.cent)),1);
b(max(size(objfinal.cent))/2+1:end)=-b(max(size(objfinal.cent))/2+1:end);
%Si objfinal.cent no te un numero parell d'elements que...?

Nfb1p=numel(obj.ds);
q=Z\b;

% CAPACIDAD
S=2*Lx*2*Ly;
%epsilon0=8.8541878176*10^(-12);
epsilon0=1;
Cteorico=S*epsilon0/d
Cexp=objfinal.ds(1:Nfb1p)*q(1:Nfb1p)/V0


% Computacion de V y E
x = linspace(-Lx,Lx,10);
y = linspace(-Ly,Ly,10); 
z = 0.05+linspace(Lz-d-Lx,Lz+Lx,10);
[xx,yy,zz] = meshgrid(x,y,z);
V=zeros(10*10*10,1);
for ii=1:max(size(objfinal.cent))
   RR=[xx(:),yy(:),zz(:)]'-objfinal.cent(:,ii); % RR la distancia entre los 
   %puntos de carga y los del meshgrid.
   RR=vecnorm(RR);
   V=V+1/(4*pi*epsilon0)*q(ii)*objfinal.ds(ii)./RR'; %Coulomb
end
Vmat=zeros(10,10,10);
for jj=1:10
   for kk=1:10
       Vmat(:,kk,jj)=V((jj-1)*100+(kk-1)*10+1:(jj-1)*100+(kk)*10); % (todas la filas, columna kk i plano jj)
   end
end
% Computation of potential at xx,yy,zz
% [Ex, Ey, Ez] = gradient(-V);
% Lx=2; Ly=2; Lz=1; d=0.000001;
xs=0; ys=0; zs=[]; % zs=Lz-d/2;
figure
slice(xx,yy,zz, Vmat, xs,ys,zs);
axis equal;
%contourslice(xx,yy,zz, Vmat, xs,ys,zs);
% coneplot(xx,yy,zz, Ex,Ey,Ez, cxx,cyy,czz, absE);
colormap jet; shading interp; colorbar