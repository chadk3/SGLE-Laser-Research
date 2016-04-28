
close all
clear all
clc

%alpha1 = 0.4 * pi;
%alpha2 = 0.4 * pi;
%alpha3 = 0.1 * pi;
%alphap = 0.44 * pi;

alpha1 = 2.6775;
alpha2 = 1.5885;
alpha3 = 1.9970;
alphap = 2.9873;

%g0 = 1.0;
g0 = 2.4291;

D = 0.4;
K = 0.1;
gamma = 0.1;
e0 = 1;
tau = 0.1;
A = 2/3;
B = 1/3;

Lt = 40;
nt = 128;
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

waterfall(t', z', abs(psisol))

temp1 = abs(psisol(length(z),:))*10;
temp2 = abs(psisol(length(z)-1,:))*10;
temp3 = abs(psisol(length(z)-2,:))*10;
temp4 = abs(psisol(length(z)-3,:))*10;

temp1r = round(temp1)/10;
temp2r = round(temp2)/10;
temp3r = round(temp3)/10;
temp4r = round(temp4)/10;

if (isequal(temp1r,temp2r,temp3r,temp4r))

    if (temp1r(1,1) > 0)
        
        if(abs(temp1r(1,1) - temp1r(1,end)) < .3)
            index1 = 1;
            index2 = length(t);
            while(temp1r(1,index1) > 0)
                index1 = index1 + 1;
            end
        
            while(temp1r(1,index2) > 0)
                index2 = index2 - 1;
            end
        
            zerotest = 0; 
        
            for k = index1:index2
                if (temp1r(1,k) > 0)
                    zerotest = 1;
                end
            end
        
            if (zerotest == 0)
                
                index3 = 1;
                for m = index1:length(t)
                    temp1rnew(1,index3) = temp1r(1,m);
                    index3 = index3 + 1;
                end
                
                for p = 1:(index1-1)
                    temp1rnew(1,index3) = temp1r(1,p);
                    index3 = index3 + 1;
                end
                
                temp1r = temp1rnew;
            end
        end
    end
            
    wsum = sum(temp1r.*t);
    rsum = sum(temp1r);
    com = wsum/rsum;

    [maxp, index] = max(temp1r);

    if (abs(round(com) - round(t(index))) <= 1)
        success = 1
    end
end