close all
clear all
clc

% Initial Conditions

B = 1/3;
gamma = 0.1;
K = 0.1;

alpha10 = 0.1 * pi;
alpha20 = 0.554 * pi;
alpha30 = 0.23 * pi;
alphap0 = 0.43 * pi;

dt = .01;
In = [0:dt:30];
dalpha = .01;
min = 0;
max = 2;
count = 1;

for j = min:dalpha:max
                
    alpha1 = j * pi;
    alpha2 = alpha20;
    alpha3 = alpha30;
    alphap = alphap0;
    w = B * In * sin(2*(alpha1 - alphap));
    Q = (1/2)*(exp(-1i*K)*(cos(2*alpha2 - 2*alpha3 - alphap) + 1i*cos(2*alpha3 - alphap)) * (1i*cos(2*alpha1 - alphap - w) - cos(alphap - w)) + exp(1i*K)*(sin(2*alpha2 - 2*alpha3 - alphap) - 1i*sin(2*alpha3 - alphap)) * (sin(alphap - w) - 1i*sin(2*alpha1 - alphap - w)));

    Ts = -1 *gamma + real(log(Q));
    
    alpha1vec(count,:) = Ts;
    count = count + 1;
end

count = 1;

for k = min:dalpha:max
    alpha1 = alpha10;
    alpha2 = k * pi;
    alpha3 = alpha30;
    alphap = alphap0;
    w = B * In * sin(2*(alpha1 - alphap));
    Q = (1/2)*(exp(-1i*K)*(cos(2*alpha2 - 2*alpha3 - alphap) + 1i*cos(2*alpha3 - alphap)) * (1i*cos(2*alpha1 - alphap - w) - cos(alphap - w)) + exp(1i*K)*(sin(2*alpha2 - 2*alpha3 - alphap) - 1i*sin(2*alpha3 - alphap)) * (sin(alphap - w) - 1i*sin(2*alpha1 - alphap - w)));

    Ts = -1 *gamma + real(log(Q));
    
    alpha2vec(count,:) = Ts;
    count = count + 1;
end

count = 1;

for m = min:dalpha:max
    alpha1 = alpha10;
    alpha2 = alpha20;
    alpha3 = m * pi;
    alphap = alphap0;
    w = B * In * sin(2*(alpha1 - alphap));
    Q = (1/2)*(exp(-1i*K)*(cos(2*alpha2 - 2*alpha3 - alphap) + 1i*cos(2*alpha3 - alphap)) * (1i*cos(2*alpha1 - alphap - w) - cos(alphap - w)) + exp(1i*K)*(sin(2*alpha2 - 2*alpha3 - alphap) - 1i*sin(2*alpha3 - alphap)) * (sin(alphap - w) - 1i*sin(2*alpha1 - alphap - w)));

    Ts = -1 *gamma + real(log(Q));
    
    alpha3vec(count,:) = Ts;
    count = count + 1;
end

count = 1;

for n = min:dalpha:max
    alpha1 = alpha10;
    alpha2 = alpha20;
    alpha3 = alpha30;
    alphap = n * pi;
    w = B * In * sin(2*(alpha1 - alphap));
    Q = (1/2)*(exp(-1i*K)*(cos(2*alpha2 - 2*alpha3 - alphap) + 1i*cos(2*alpha3 - alphap)) * (1i*cos(2*alpha1 - alphap - w) - cos(alphap - w)) + exp(1i*K)*(sin(2*alpha2 - 2*alpha3 - alphap) - 1i*sin(2*alpha3 - alphap)) * (sin(alphap - w) - 1i*sin(2*alpha1 - alphap - w)));

    Ts = -1 *gamma + real(log(Q));
    
    alphapvec(count,:) = Ts;
    count = count + 1;
end

for p = 1:(max/dalpha);
    plot(In, alpha1vec(p,:), 'g-', In, alpha2vec(p,:), 'b-', In, alpha3vec(p,:), 'r-', In, alphapvec(p,:), 'k-')
    axis([0 30 -3 3])
    title(['Plot for alpha = ' num2str(p*dalpha) '\pi'])
    legend('alpha1', 'alpha2', 'alpha3', 'alphap')
    pause
end

