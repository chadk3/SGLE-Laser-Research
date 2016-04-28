
% Simulates the sinusoidal Ginzburg-Landau equation.

close all
clear all
clc

% g0 is the nondimensional pumping strength.
g0 = .75;

% D is the averaged group velocity dispersion of the cavity, and is 
% positive for anomalous dispersion and negative for normal dispersion.
D = 0.4;

% alpha1 is the angle of the left QWP, alphap the polarizer, alpha2 
% the right QWP, and alpha3 the HWP.

% alpha1 = 21.056419 * pi;
% alpha2 = 24.078937 * pi;
% alpha3 = 25.45597 * pi;
% alphap = 28.007656 * pi;

alpha1 = 0;
alpha2 = 0.82 * pi;
alpha3 = 0.1 * pi;
alphap = 0.45 * pi;

% K is the birefringence.
K = 0.1;

% gamma measures the distributed losses caused by the output coupling
% and the fiber attenuation.
gamma = 0.1;

% e0 is the saturating energy of the gain medium.
e0 = 1;

% tau characterizes the bandwidth of the pump.
tau = 0.1;

% A is the cross-phase modulations
A = 2/3;

% B is the four-wave mixing
B = 1/3;

Lt = 40;
%nt = 2048;
%nt = 1024;
%nt = 512;
%nt = 256;
%nt = 128;
nt = 64;
t2 = linspace(-Lt/2, Lt/2, nt+1);
t = t2(1:nt); 

z=0:100:2000;

% Initial Conditions
%psi = sech(0.5*t);
%psit = fft(psi);
noise = 1;
psit = noise*(randn(1,nt) + 1i*rand(1,nt));
psi = ifft(psit);


% Spectral k values
kt = (2*pi/Lt)*[0:(nt/2-1) (-nt/2): -1].';

[tj, psitsol] = ode45('sgle_rhs', z, psit, [], kt, t, gamma, D, g0, e0, alpha1, alpha2, alpha3, alphap, B, K, tau);

for j = 1:length(z)
    psisol(j,:) = ifft(psitsol(j,:));
end

maxpulse = max(abs(psisol(end,:)));
E = trapz(t,abs(psisol(end,:)).^2);

waterfall(t', z', abs(psisol))
title(['g0 = ' num2str(g0) ' alpha1 = ' num2str(alpha1/pi) '\pi alpha2 = ' num2str(alpha2/pi) '\pi alpha3 = ' num2str(alpha3/pi) '\pi alphap = ' num2str(alphap/pi) '\pi max pulse = ' num2str(maxpulse) ' E = ' num2str(E) ' score = ' num2str(maxpulse + E)])

figure(2)
plot(abs(psisol(end,:)),'k-o')

% Objective Function
obj = E/kurtosis(abs(psisol(end,:)));
obj2 = E/kurtosis(abs(fftshift(fft(abs(psisol(end,:))))));
obj3 = E/trapz(t,abs(psisol(end,:)).^4);




