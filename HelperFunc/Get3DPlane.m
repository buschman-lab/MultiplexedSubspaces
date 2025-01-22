function [X,Y,Z] = Get3DPlane(x,scaleFactor)
%from Panichello and Buschman 2021 
%x is a n x 3D projection into low dimensional space
%to plot, use PATCH(X,Y,Z,C)
%get the eigenvalues
if nargin <2; scaleFactor = 25; end
eigvecsUp = pca(x);
mn = mean(x(:,1:3),1);
x = eigvecsUp(1,1:2).*scaleFactor;
y = eigvecsUp(2,1:2).*scaleFactor;
z = eigvecsUp(3,1:2).*scaleFactor;

X = [x -x]+mn(1);
Y = [y -y]+mn(2);
Z = [z -z]+mn(3);
end