close all
clear all
clc

% Initial Conditions
B = 1/3;

% Parameters for Q1
alpha11 = 0 * pi;
alpha12 = 0 * pi;
alpha13 = 0 * pi;
alpha1p = 0 * pi;
K1 = 0.1;
gamma = .1;

dt = .01;
In = [0:dt:30];
dalpha = .1;

% Constructs all possible transmission curves for every permutation of
% the four alphas from 0 to 2*pi in dalpha increments

for j = 0:dalpha:2
    
    alpha1p = j * pi;

    for k = 0:dalpha:2
    
        alpha13 = k * pi;
        
        for m = 0:dalpha:2
        
            alpha12 = m * pi;
    
            for n = 0:dalpha:2
                
                alpha11 = n * pi;
    
                w1 = B * In * sin(2*(alpha11 - alpha1p));
                Q = (1/2)*(exp(-1i*K1)*(cos(2*alpha12 - 2*alpha13 - alpha1p) + 1i*cos(2*alpha13 - alpha1p)) * (1i*cos(2*alpha11 - alpha1p - w1) - cos(alpha1p - w1)) + exp(1i*K1)*(sin(2*alpha12 - 2*alpha13 - alpha1p) - 1i*sin(2*alpha13 - alpha1p)) * (sin(alpha1p - w1) - 1i*sin(2*alpha11 - alpha1p - w1)));
                Qp = diff(Q)/dt;
                Qpp = diff(Qp)/dt;

                Ts = -1 *gamma + real(log(Q));

                plot(In, Ts, 'g-')
                axis([0 30 -3 3])
                title(['alpha1 = ' num2str(alpha11 / pi) '\pi, alpha2 = ' num2str(alpha12 / pi) '\pi, alpha3 = ' num2str(alpha13 / pi) '\pi, alphap = ' num2str(alpha1p / pi) '\pi' ])
                pause
            end
        end
    end
end
