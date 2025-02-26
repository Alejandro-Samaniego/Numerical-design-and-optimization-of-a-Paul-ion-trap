function [Posmattop,Posmatbottom,Posmatring,obj]=hiperbolic(delta,R,it,v)

%% OBS
% N== numero de puntos de discretizacion por circumferencia. No lo podemos
% entrar pk N es diferente para cada circumf... queremos misma densidad de
% ptos por cirucmf.

% Delta== incremento de radio de las circumferencias concentricas.
% R == radio de la 1a circumf, mas pequeña (afecta a los top bottom)
% it==numero de circumferencias concentricas que queremos
% v==vector con cordenadas x e y donde "empiezan" los casquetes hiperbolicos
% v(1) son las z de las puntas de los casquetes top/bottom
% v(2) es la distancia a la punta del casquete central desde el origen
% V(3) v(4) es el alto hasta el que queremos llegar

%% Top and bottom
R=[R:delta(1):R+it*delta(1)];
N=round(2*pi*R/delta(1));
pdc=zeros(1,sum(N));
pdc(1:N(1))=R(1)*exp(1j*linspace(0,2*pi,N(1))); %pdc==puntos de circumferencia
escalon=ones(1:it+1);
for ii=1:it
    pdc(escalon(1:ii)*N(1:ii)'+1:escalon(1:ii+1)*N(1:ii+1)')=R(ii+1)*exp(1j*linspace(0,2*pi-2*pi/N(ii+1),N(ii+1))); %pdc==puntos de circumferencia
    %-delta con tal de que no haya puntos repetidos????????????
end
z2=(2*v(1)^2+abs(pdc).^2)/2;
zmas=sqrt(z2); 
zmenos=-sqrt(z2); %rama de abajo

Posmattop=[real(pdc);imag(pdc);zmas];
obj.Posmattop=Posmattop;
obj.topoltop = delaunay(Posmattop(1,:)', Posmattop(2,:)')';

Posmatbottom=[real(pdc);imag(pdc);zmenos];
obj.Posmatbottom=Posmatbottom;
obj.topolbottom = delaunay(Posmatbottom(1,:)', Posmatbottom(2,:)')';


%% Ring
[fifi,zz]=meshgrid(0:delta(2):2*pi, v(3):delta(1):v(4));
r2=v(2)^2+2*zz.^2;
r=sqrt(r2); %radios

polar=r.*exp(1j*fifi);
polar1=polar(1:floor(end/2)+1,1:end-1);
polar2=polar(floor(end/2)+1:end,1:end-1);
zz1=zz(1:floor(end/2)+1,1:end-1);
zz2=zz(floor(end/2)+1:end,1:end-1);
Posmatring1=[real(polar1(:))';imag(polar1(:))';zz1(:)'];
Posmatring2=[real(polar2(:))';imag(polar2(:))';zz2(:)'];
Posmatring=[Posmatring1 Posmatring2];
obj.Posmatring=Posmatring;

% DT = delaunayTriangulation(x,y,C) specifies the edge constraints in a matrix C.
topol2_1 = delaunayTriangulation(Posmatring1(1,:)', Posmatring1(2,:)')';
topol2_2 = delaunayTriangulation(Posmatring2(1,:)', Posmatring2(2,:)')';
% topol2_1 = delaunay(Posmatring1(1,:)', Posmatring1(2,:)')';
% topol2_2 = delaunay(Posmatring2(1,:)', Posmatring2(2,:)')';


%% Eliminar triangulos no deseados con cross product
topol2_1=topol2_1.ConnectivityList;
topol2_2=topol2_2.ConnectivityList;
topol2_2=topol2_2+length(Posmatring1);
obj.topolring=[topol2_1' topol2_2'];
v1 = obj.Posmatring(:,obj.topolring(1,:));
v2 = obj.Posmatring(:,obj.topolring(2,:));
v3 = obj.Posmatring(:,obj.topolring(3,:));
normales = cross(v2-v1, v3-v1);
zz=[0 0 1]'; %% Quizas se debe generalizar o ¿poner esta z como entrada a la función...?
for jj=1:length(obj.topolring)
     if norm(cross(normales(:,jj),zz))==0
          obj.topolring(:,jj)=[0 0 0]';
     else

     end
end
[~,col]=find(obj.topolring==0);
obj.topolring(:,col)=[];

%% CONCATENAMOS TODOS LOS OBJS
obj.Posmat=[Posmattop Posmatbottom Posmatring];
topoltop=obj.topoltop;
topolbottom=obj.topolbottom;
topolring=obj.topolring;
topolbottom=topolbottom+length(Posmattop);
topolring=topolring+length(Posmattop)+length(Posmatbottom);
obj.topol=[topoltop topolbottom topolring];

% obj.Posmat=[Posmattop Posmatring Posmatbottom];
% topoltop=obj.topoltop;
% topolbottom=obj.topolbottom;
% topolring=obj.topolring;
% topolring=topolring+length(Posmattop);
% topolbottom=topolbottom+length(Posmattop)+length(Posmatring);
% obj.topol=[topoltop topolring topolbottom];
end