
close all
clear all
clc

alpha1 = 0 * pi;
alpha2 = 0.82 * pi;
alpha3 = 0.1 * pi;
alphap = 0.44 * pi;

g0 = 1.0;
D = 0.4;
K = 0.1;
gamma = 0.1;
e0 = 1;
tau = 0.1;
A = 2/3;
B = 1/3;

Lt = 40;

%nt = 128, 256, 512, 1024, 2048, ...
nt = 128;
t2 = linspace(-Lt/2, Lt/2, nt+1);
t = t2(1:nt); 

z=0:20:500;

% Initial Conditions
%psi = sech(0.5*t);
%psit = fft(psi);
noise = 1;
psit = noise*(randn(1,nt) + 1i*rand(1,nt));
psi = ifft(psit);

kt = (2*pi/Lt)*[0:(nt/2-1) (-nt/2): -1].';

dalpha = .05;
count = 0;

for m = 0:dalpha:1
    
    alphap = m * pi;
    
    [tj, psitsol] = ode45('sgle_rhs', z, psit, [], kt, t, gamma, D, g0, e0, alpha1, alpha2, alpha3, alphap, B, K, tau);

    for j = 1:length(z)
        psisol(j,:) = ifft(psitsol(j,:));
    end
    
    count = count + 1;
    
    figure(count), waterfall(t', z', abs(psisol))
    title(['Plot for alphap = ' num2str(m) '\pi'])
end

