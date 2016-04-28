close all
clear all
clc

% Initial Conditions

B = 1/3;
gamma = 0.1;
K = 0.1;

alpha11 = 0.1 * pi;
alpha21 = 0.554 * pi;
alpha31 = 0.23 * pi;
alphap1 = 0.43 * pi;

alpha12 = 0.1 * pi;
alpha22 = 0.44 * pi;
alpha32 = 0.13 * pi;
alphap2 = 0.63 * pi;

alpha13 = 0.1 * pi;
alpha23 = 0.4 * pi;
alpha33 = 0.7 * pi;
alphap3 = 0.3 * pi;

alpha14 = 0.1 * pi;
alpha24 = 0.8 * pi;
alpha34 = 0.8 * pi;
alphap4 = 0.8 * pi;

dt = .01;
In = [0:dt:30];
dalpha = .01;
min = 0;
max = 2;
count = 1;

for j = min:dalpha:max
                
    alpha1 = j * pi;
    alpha2 = alpha21;
    alpha3 = alpha31;
    alphap = alphap1;
    w = B * In * sin(2*(alpha1 - alphap));
    Q = (1/2)*(exp(-1i*K)*(cos(2*alpha2 - 2*alpha3 - alphap) + 1i*cos(2*alpha3 - alphap)) * (1i*cos(2*alpha1 - alphap - w) - cos(alphap - w)) + exp(1i*K)*(sin(2*alpha2 - 2*alpha3 - alphap) - 1i*sin(2*alpha3 - alphap)) * (sin(alphap - w) - 1i*sin(2*alpha1 - alphap - w)));

    Ts = -1 *gamma + real(log(Q));
    
    alpha11vec(count,:) = Ts;
    count = count + 1;
end

count = 1;

for k = min:dalpha:max
    
    alpha1 = k * pi;
    alpha2 = alpha22;
    alpha3 = alpha32;
    alphap = alphap2;
    w = B * In * sin(2*(alpha1 - alphap));
    Q = (1/2)*(exp(-1i*K)*(cos(2*alpha2 - 2*alpha3 - alphap) + 1i*cos(2*alpha3 - alphap)) * (1i*cos(2*alpha1 - alphap - w) - cos(alphap - w)) + exp(1i*K)*(sin(2*alpha2 - 2*alpha3 - alphap) - 1i*sin(2*alpha3 - alphap)) * (sin(alphap - w) - 1i*sin(2*alpha1 - alphap - w)));

    Ts = -1 *gamma + real(log(Q));
    
    alpha12vec(count,:) = Ts;
    count = count + 1;
end

count = 1;

for m = min:dalpha:max
    
    alpha1 = m * pi;
    alpha2 = alpha23;
    alpha3 = alpha33;
    alphap = alphap3;
    w = B * In * sin(2*(alpha1 - alphap));
    Q = (1/2)*(exp(-1i*K)*(cos(2*alpha2 - 2*alpha3 - alphap) + 1i*cos(2*alpha3 - alphap)) * (1i*cos(2*alpha1 - alphap - w) - cos(alphap - w)) + exp(1i*K)*(sin(2*alpha2 - 2*alpha3 - alphap) - 1i*sin(2*alpha3 - alphap)) * (sin(alphap - w) - 1i*sin(2*alpha1 - alphap - w)));

    Ts = -1 *gamma + real(log(Q));
    
    alpha13vec(count,:) = Ts;
    count = count + 1;
end

count = 1;

for n = min:dalpha:max
    
    alpha1 = n * pi;
    alpha2 = alpha24;
    alpha3 = alpha34;
    alphap = alphap4;
    w = B * In * sin(2*(alpha1 - alphap));
    Q = (1/2)*(exp(-1i*K)*(cos(2*alpha2 - 2*alpha3 - alphap) + 1i*cos(2*alpha3 - alphap)) * (1i*cos(2*alpha1 - alphap - w) - cos(alphap - w)) + exp(1i*K)*(sin(2*alpha2 - 2*alpha3 - alphap) - 1i*sin(2*alpha3 - alphap)) * (sin(alphap - w) - 1i*sin(2*alpha1 - alphap - w)));

    Ts = -1 *gamma + real(log(Q));
    
    alpha14vec(count,:) = Ts;
    count = count + 1;
end

for p = 1:(max/dalpha);
    plot(In, alpha11vec(p,:), 'g-', In, alpha12vec(p,:), 'b-', In, alpha13vec(p,:), 'r-', In, alpha14vec(p,:), 'k-')
    axis([0 30 -3 3])
    title(['Plot for alpha = ' num2str(p*dalpha) '\pi'])
    pause
end

