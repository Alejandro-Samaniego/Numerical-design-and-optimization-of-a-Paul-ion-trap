function []=hiperbolic(N,v,xcap, xring)
% v vector con x0==r0 y z0, parametros de la hiperbola
% v(1)= x0, v(2)= z0
% N numero de puntos de discretizacion de los 2pi radianes
% ATENCION: x representa la r. X tiene que ser un vect. Xcap y xring son las 
% regiones donde definimos las ramas de las hipèrbolas correspondientes a
% los caps y el anillo respectivamente. xcap y xring entra en fila.

% Rotacion the un angulo theta. La aplicaremos n VECES para trotar n*theta
% grados.
theta=2*pi/N;
Rot=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

%End cap 1
z2=(2*v(2)^2+xcap.^2)/2;
zmas=sqrt(z2); %rama de arriba
Pmat=[xcap; zeros(1,length(xcap));zmas];


%End cap 2
zmenos=-sqrt(z2); %rama de abajo

%Ring
zr=(-v(1)^2+xring.^2)/2;
zring=sqrt(zr);


end