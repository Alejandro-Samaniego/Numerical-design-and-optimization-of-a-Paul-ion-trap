function [obj]=geomtopol(r,N)
%Entramos un vector de posiciones con la x inicial, x final, y inicial, y
%final, z, respectivamente.
x = linspace(r(1),r(2),N);
y = linspace(r(3),r(4),N);
[xx,yy] = meshgrid(x,y); 
zz = r(5)*xx./xx;
obj.vertex = [xx(:) yy(:) zz(:)]'; % xx(:) devuelve todos los elementos del 
% meshgrid en un vector columna
obj.topol = delaunay(xx(:), yy(:))';
end