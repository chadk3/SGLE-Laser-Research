
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

dalpha = .2;
count1 = 0;
count2 = 0;
count3 = 0;

for j = 0:dalpha:.4
    
    alpha2 = j * pi;
    
    for k = 0:dalpha:1
    
        alpha1 = k * pi;
    
        [tj, psitsol] = ode45('sgle_rhs', z, psit, [], kt, t, gamma, D, g0, e0, alpha1, alpha2, alpha3, alphap, B, K, tau);

        for m = 1:length(z)
            psisol(m,:) = ifft(psitsol(m,:));
            count1 = count1 + 1;
            
            if (2 <= count1 <= length(z))
                if(psisol(count1,:) == psisol(count1-1,:))
                    count2 = count2 + 1;
                end
            end
    
            if (count2 >= 5)
            
                count3 = count3 + 1;
                max = max(psisol(end,:));
                energy = trapz(t,abs(psi).^2);
            
                alphas(count3, :) = [alpha1, alpha2, alpha3, alphap, max, energy];
            end
        end
    end
end