clear
%% 1.- Parámetros, ecuaciones, condiciones iniciales & rango temporal
% m = 1;                             
% g = 9.8;                            
% theta0 = 0.174;                    
% omega0 = 0;                        
tf = 10;                          
FR=50;
p1=1280;
p2=0.9*720;
e1=0;
e2=0;
% lvec = [];  
% cvec = []; 
% Evec = [];
%Programamos un bucle para obtener los valores de los 8 péndulos
%diferentes y los almacenamos en los vectores vacíos definidos previamente
%para almacenar los datos de forma compacta.
% for n=0:7
% T=60/(20+n);
% lvec=[lvec,g/(4*pi^2)*T^2];
% cvec = [cvec,g/lvec(end)];                     
% Evec = [Evec,m*lvec(end)^2*omega0^2/2 + m*g*lvec(end)*(1-cos(theta0))];                              
% end

%% 2.- Características del plot (ventana, colores, tamaño de las fuentes y los puntos, grosor de la líneas).
figh=figure('WindowState','fullscreen');
axis equal
axis manual
AR = 16/9;                  %'Screen aspect ratio'
%Se reserva un espacio a la derecha de la ventana con el objetivo de
%insertar texto más adelante.
%El objetivo de la siguiente línea es encajar las componentes que
%aparecerán en el vídeo de una forma organizada y estética.
axis([-0.07 0.07 -0.07 0.07 -0.07 0.07])
ax = gca
ax.GridColor = [0.25 0.25 0.25];
ax.GridAlpha = 0.5;
box on
bgColor = 1.5*[0.4 0.4 0.4]; %Color del backgroud (en el modelo de colores RGB).
fgColor = 1.5*[0.1 0.1 0.1]; %Color del Foreground.
bggColor = [1 1 1]; %Color de fondo del Figure.
mColor = [0.635 0.0780 0.184];        %Color de la masa (granate).
% sColor = [0.15 0.15 0.15];        %Color del soporte.
% pColor = [1 1 0.75];        %Color del péndulo.
% ma= 17;                     %Tamaño de la bola de soporte.
ms = 11;                    %Tamaño de las masas.
% fs = 15;                    %Tamaño de la fuente.
% lw = 2;                     %Grosor de las líneas.
set(gca,'Color',bgColor)
set(gcf,'color',bggColor);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[],'tickDir','none')
set(figh, 'ToolBar', 'none');
hold on
grid on
%% 3.- Características del vídeo
nameFile = 'IonPaulTrap';                 
vwo = VideoWriter(nameFile,'MPEG-4');         
FR = 50;                                               
prescribedTimes = 0:1/FR:tf;  
vwo.Quality = 95;             
vwo.FrameRate = FR;         
open(vwo);
%% 4.- Computación numérica del movimiento del péndulo con el comando ode45
%El siguiente bucle permite definir la ecuación diferencial que seguirá
%cada uno de entre los ocho péndulos que hemos definido. Los vectores
%obtenidos tras la resolución de las respectivas edos se almacenan en una
%celda para llamarlos más adelante.
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
p=[-0.0048;-0.0039;0.0058;-0.1698;-0.0276;-0.0391];
%p=rand(6,n)-rand(6,n); 
qion=1.602*10^(-19)*ones(1,n);
%p([1:3],:)=0.025*p([1:3],:); 
m=9.012182*1.6605*10^(-27)*ones(1,n); 
deltat=2.5*10^(-4); 
N=500; 
w=1000;
posmat=ionprimeroAC(p,qion,m,deltat,N,Vmat,xx,yy,zz,w);
%% 5.- Frame inicial
%Creamos un bucle para dibujar cada péndulo en el instante t=0 y adaptamos
%las coordenadas añadiendo la componente 3D.
camera=cross(p(4:6),[0;0;1]);
electrodo=patch('Faces',obj.topol','Vertices',obj.Posmat','CData',q,'FaceColor','flat','EdgeColor',[0.2 0.2 0.2]);
alpha(0.2)
colormap jet; 
Ion=plot3(posmat(1,1),posmat(2,1),posmat(3,1),'o','MarkerFaceColor',mColor,'MarkerEdgeColor',mColor);
view(45,25)
%A continuación añadimos el texto que aparecerá en el vídeo. Los parámetros
%que varían con el tiempo habrá que volver a definirlos más adelante para
%el resto de frames.
%IMPORTANTE: Usamos el interpretador de texto LaTeX para los símbolos matemáticos.
frame = getframe(figh,[e1 e2 p1 p2]);
writeVideo(vwo,frame);
%% 6.- Loop para el resto de frames del vídeo
%Definimos todo lo que el bucle necesitará borrar.
for i=1:N-1
  delete(Ion);
  Ion=plot3(posmat(1,1:i+1),posmat(2,1:i+1),posmat(3,1:i+1),'-','Color',mColor);
  view(45,25)
   frame = getframe(figh,[e1 e2 p1 p2]);
   writeVideo(vwo,frame);
end

% %7.- Loop para obtener un segundo de un 'frame' estático con la figura ya dibujada
% nFrames = FR*1;  %Número de 'frames' en 5 segundos
% for i=1:nFrames
%    frame = getframe;
%    writeVideo(vwo,frame);
% end
close(vwo);
