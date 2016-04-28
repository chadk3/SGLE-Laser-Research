
close all;
clear all;
clc;

% Just need a loop to increment gain

g0 = 1;
D = 0.4;
K = 0.1;
gamma = 0.1;
e0 = 1;
tau = 0.1;
A = 2/3;
B = 1/3;

Lt = 60;
nt = 256; %or 512
t2 = linspace(-Lt/2, Lt/2, nt+1);
t = t2(1:nt);
kt = (2*pi/Lt)*[0:(nt/2-1) (-nt/2): -1].';

% Initial Conditions
noise = 1;
psit = noise*(randn(1,nt) + 1i*rand(1,nt));
psi = ifft(psit);

anscount = 1;
topscore = 0;

dt = 1;
objint_old = 0;
objint = 1;
dt_old = 0;
tol = 0.01; % tolerance to test objective functions

while abs(objint-objint_old) > tol
    
    ttor = [0:dt+dt_old:100];
    
    z=0:100:2000;
    
    alpha10 = 0;
    alpha20 = 0;
    alpha30 = 0;
    alphap0 = 0;
    
    omega1 = sqrt(.13);
    omega2 = sqrt(.17);
    omega3 = sqrt(.19);
    omegap = sqrt(.23);
    
    alpha1t = omega1*pi*ttor + alpha10;
    alpha2t = omega2*pi*ttor + alpha20;
    alpha3t = omega3*pi*ttor + alpha30;
    alphapt = omegap*pi*ttor + alphap0;
    
    for i = 1:length(ttor)
        
        alpha1 = mod(alpha1t(1,i),2);
        alpha2 = mod(alpha2t(1,i),2);
        alpha3 = mod(alpha3t(1,i),2);
        alphap = mod(alphapt(1,i),2);
        
        [tj, psitsol] = ode45('sgle_rhs', z, psit, [], kt, t, gamma, D, g0, e0, alpha1, alpha2, alpha3, alphap, B, K, tau);
        
        for j = 1:length(z)
            psisol(j,:) = ifft(psitsol(j,:));
        end
        
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
                    
                    while((index1 < length(t)) && (temp1r(1,index1) > 0))
                        index1 = index1 + 1;
                    end
                    
                    while((index2 > 0) && (temp1r(1,index2) > 0))
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
                
                %mode locks
                maxpulse = max(abs(psisol(end,:)));
                E = trapz(t,abs(psisol(end,:)).^2);
                score = maxpulse + E;
                
                answers(anscount,:) = [alpha1/pi alpha2/pi alpha3/pi alphap/pi maxpulse E score];
                
                if (score > topscore)
                    
                    topscore = score;
                    topanswer = [alpha1/pi alpha2/pi alpha3/pi alphap/pi maxpulse E score];
                    
                end
                
                %obj(anscount) = E/kurtosis(abs(fftshift(fft(abs(psisol(end,:))))));
                obj(anscount) = E/kurtosis(abs(psisol(end,:)));
                
                anscount = anscount + 1;
                
            end
            
        end
        
        waterfall(t', z', abs(psisol))
        
    end
    
    %objint = trapz(t,obj);
    objint = trapz(anscount-1,obj);
    
    if abs(objint - objint_old) > tol;
        dt = dt/2;
        objint_old = objint;
    end
    
end

% Xing's Objective Function
% obj=eng/kurtosis(abs(fftshift(fft(phi(:,end)))));
% Where ‘eng' is the pulse energy and
%'phi(:,end)' is the steady state pulse shape (which
% is the last snapshot)






