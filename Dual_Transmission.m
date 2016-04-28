close all
clear all
clc

% Initial Conditions
B = 1/3;
gamma = .1;

% Parameters for Q1
alpha11 = 0.1 * pi;
alpha12 = 0.554 * pi;
alpha13 = 0.23 * pi;
alpha1p = 0.43 * pi;
K1 = 0.1;

% Parameters for Q2
alpha21 = 0 * pi;
alpha22 = 0.84 * pi;
alpha23 = 0.1 * pi;
alpha2p = 0.45 * pi;
K2 = 0.1;

dt = .01;
In = [0:dt:30];

w1 = B * In * sin(2*(alpha11 - alpha1p));
Q1 = (1/2)*(exp(-1i*K1)*(cos(2*alpha12 - 2*alpha13 - alpha1p) + 1i*cos(2*alpha13 - alpha1p)) * (1i*cos(2*alpha11 - alpha1p - w1) - cos(alpha1p - w1)) + exp(1i*K1)*(sin(2*alpha12 - 2*alpha13 - alpha1p) - 1i*sin(2*alpha13 - alpha1p)) * (sin(alpha1p - w1) - 1i*sin(2*alpha11 - alpha1p - w1)));

w2 = B * In * sin(2*(alpha21 - alpha2p));
%Q2 = (1/2)*(exp(-1i*K2)*(cos(2*alpha22 - 2*alpha23 - alpha2p) + 1i*cos(2*alpha23 - alpha2p)) * (1i*cos(2*alpha21 - alpha2p - w2) - cos(alpha2p - w2)) + exp(1i*K2)*(sin(2*alpha22 - 2*alpha23 - alpha2p) - 1i*sin(2*alpha23 - alpha2p)) * (sin(alpha2p - w2) - 1i*sin(2*alpha21 - alpha2p - w2)));
Q2 = ones(size(In));

Q = Q1.*Q2;

Qp = diff(Q)/dt;
Qpp = diff(Qp)/dt;

delta = gamma - log(abs(Q(1,1)));
betta = real(Qp(1,1)/Q(1,1));
mu = real((Q(1,1)*Qpp(1,1) - Qp(1,1)^2)/Q(1,1)^2)/2;

Tc = -1*delta + betta * In;
Tcq = -1*delta + betta * In + mu * In.^2;
Ts = -1 *gamma + real(log(Q));


plot(In, Tc, 'b.', In, Tcq, 'r--', In, Ts, 'g-')
axis([0 30 -10 10])