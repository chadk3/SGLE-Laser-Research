
close all;
clear all;
clc;

t = [0:0.1:100];

theta10 = 0;
theta20 = 0;

omega1 = sqrt(.17);
omega2 = sqrt(.19);

theta1 = omega1*pi*t + theta10;
theta2 = omega2*pi*t + theta20;

% Create Torus

u = (0:0.01:1)*2*pi;
v = (0:0.01:1)*2*pi;

[U,V] = meshgrid(u,v);

% a is outer radius of torus
% b is inner radius of torus

a = 3;
b = 1;

% Parametric equations for torus

X = (a+b.*cos(V)).*cos(U);
Y = (a+b.*cos(V)).*sin(U);
Z = b.*sin(V);

figure
surf(X,Y,Z,'EdgeColor','none')
colormap(bone)
alpha(0.5)
axis equal
axis off
hold on

% sample points

X = (a+b.*cos(theta2)).*cos(theta1);
Y = (a+b.*cos(theta2)).*sin(theta1);
Z = b.*sin(theta2);

figure(1)
plot3(X,Y,Z,'k.','MarkerSize',6)
%set(gcf,'Position',[100 200 200 150])

