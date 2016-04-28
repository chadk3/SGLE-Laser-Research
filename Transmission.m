
close all
clear all
clc

% Initial Conditions

B = 1/3;
gamma = 0.1;
K = 0.1;

alpha1 = 0.1 * pi;
alpha2 = 0.82 * pi;
alpha3 = 0.1 * pi;
alphap = 0.44 * pi;


dt = .01;
In = [0:dt:30];
dalpha = .01;

% Changes one alpha from 0 to 2*pi in dalpha increments while holding the
% other parameters constant

for n = 0.75:dalpha:0.85
                
    alpha1 = n * pi;
    w = B * In * sin(2*(alpha1 - alphap));
    Q = (1/2)*(exp(-1i*K)*(cos(2*alpha2 - 2*alpha3 - alphap) + 1i*cos(2*alpha3 - alphap)) * (1i*cos(2*alpha1 - alphap - w) - cos(alphap - w)) + exp(1i*K)*(sin(2*alpha2 - 2*alpha3 - alphap) - 1i*sin(2*alpha3 - alphap)) * (sin(alphap - w) - 1i*sin(2*alpha1 - alphap - w)));

    Ts = -1 *gamma + real(log(Q));

    plot(In, Ts, 'g-')
    axis([0 30 -3 3])
    title(['Plot for alpha1 = ' num2str(n) '\pi'])
    pause
end

