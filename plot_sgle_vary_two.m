
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

z=0:40:1000;

% Initial Conditions
%psi = sech(0.5*t);
%psit = fft(psi);
noise = 1;
psit = noise*(randn(1,nt) + 1i*rand(1,nt));
psi = ifft(psit);

kt = (2*pi/Lt)*[0:(nt/2-1) (-nt/2): -1].';

dalpha = .01;

count = 0;

for j = 0:dalpha:1
    
    alpha2 = j * pi;
    
    for k = 0.1:dalpha:0.9
        
        alpha1 = k * pi;
    
        [tj, psitsol] = ode45('sgle_rhs', z, psit, [], kt, t, gamma, D, g0, e0, alpha1, alpha2, alpha3, alphap, B, K, tau);

        for m = 1:length(z)
            psisol(m,:) = ifft(psitsol(m,:));
        end
        
        maxpulse = max(abs(psisol(end,:)));
        E = trapz(t,abs(psisol(end,:)).^2);
        
        count = count + 1;
        
        filename = [ 'Figure' num2str(count) '.jpeg' ];
        
        figure(count), waterfall(t', z', abs(psisol))
        title(['alpha1 = ' num2str(alpha1/pi) '\pi alpha2 = ' num2str(alpha2/pi) '\pi alpha3 = ' num2str(alpha3/pi) '\pi alphap = ' num2str(alphap/pi) '\pi max pulse = ' num2str(maxpulse) ' E = ' num2str(E)])

        saveas(gcf, ['Fig' num2str(count), '.png'])
        
        close all
    end
end