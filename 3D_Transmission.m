
close all
clear all
clc

% Initial Conditions

B = 1/3;
gamma = 0.1;
K = 0.1;

alpha1 = 0.1 * pi;
alpha2 = 0.554 * pi;
alpha3 = 0.23 * pi;
alphap = 0.43 * pi;

dt = .01;
In = [0:dt:30];
dalpha = .01;
alpha1vec = [0:dalpha:2*pi];
alpha2vec = alpha2 * ones(size(alpha1vec));
alpha3vec = alpha3 * ones(size(alpha1vec));
alphapvec = alphap * ones(size(alpha1vec));

w = B * In.* sin(2*(alpha1 - alphap));
Q = (1/2)*(exp(-1i*K)*(cos(2*alpha2 - 2*alpha3 - alphap) + 1i*cos(2*alpha3 - alphap)) * (1i*cos(2*alpha1 - alphap - w) - cos(alphap - w)) + exp(1i*K)*(sin(2*alpha2 - 2*alpha3 - alphap) - 1i*sin(2*alpha3 - alphap)) * (sin(alphap - w) - 1i*sin(2*alpha1 - alphap - w)));

Ts = -1 *gamma + real(log(Q));

surf(In, alpha1vec, Ts)
