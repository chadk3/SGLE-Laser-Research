
close all
clear all
clc

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
z=0:40:1000;

noise = 1;
psit = noise*(randn(1,nt) + 1i*rand(1,nt));
psi = ifft(psit);

kt = (2*pi/Lt)*[0:(nt/2-1) (-nt/2): -1].';

Top = zeros(10,5);

% Top 10 [gain, alpha1, alpha2, alpha3, alphap]
Topvar(1,:) = [1 0.3 0 0.1 .44];
Topvar(2,:) = [2.7007 0.2068 0.5548 0.8793 0.5579];
Topvar(3,:) = [1 0.8 0.82 0.1 0.44];
Topvar(4,:) = [1 0.3511 0.054297 0.7087 0.99293];
Topvar(5,:) = [1 0.40238 0.08981 0.68377 0.54689];
Topvar(6,:) = [1 0.38185 0.63529 0.046564 0.23207];
Topvar(7,:) = [1 0.23506 0.27546 0.95161 0.34667];
Topvar(8,:) = [2.7969 0.6152 0.0269 0.3225 0.4638];
Topvar(9,:) = [1.8679 0.412816 0.01477 0.702765 0.506749];
Topvar(10,:) = [1.9600 0.6360 0.5256 0.2596 0.0512];

Topscore = zeros(10,3);

%[Max Pulse, Energy, Max Pulse + Energy]
Topscore(1,:) = [2.5849 3.5134 6.0983];
Topscore(2,:) = [2.7412 3.9676 6.7088];
Topscore(3,:) = [2.6032 3.5544 6.1576];
Topscore(4,:) = [2.5474 3.5107 6.0581];
Topscore(5,:) = [2.5368 3.4559 5.9927];
Topscore(6,:) = [2.5733 3.5498 6.1231];
Topscore(7,:) = [2.5983 3.4226 6.0209];
Topscore(8,:) = [2.7158 3.9508 6.6666];
Topscore(9,:) = [2.5505 3.4767 6.0273];
Topscore(10,:) = [2.6561 3.4096 6.0657];

count = 1;
count2 = 19;

for q = 1:50
    
    Toppert = zeros(10,5);

    for j = 1:10
        a = .8;
        b = 3.58;
        r = a + (b-a).*rand(1);
        pert = rand(1,5);
        Toppert(j,:) = pert;
        Toppert(j,1) = r;
    end

    for k = 1:10
        %g0 = 1;
        g0 = Toppert(k,1);
        alpha1 = Toppert(k,2) * pi;
        alpha2 = Toppert(k,3) * pi;
        alpha3 = Toppert(k,4) * pi;
        alphap = Toppert(k,5) * pi;
    
        [tj, psitsol] = ode45('sgle_rhs', z, psit, [], kt, t, gamma, D, g0, e0, alpha1, alpha2, alpha3, alphap, B, K, tau);
    
        for m = 1:length(z)
            psisol(m,:) = ifft(psitsol(m,:));
        end
        
        figure(1), waterfall(t', z', abs(psisol))
        
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
                
                [mintopscore,minindex] = min(Topscore(:,3)');
                
                if(score > mintopscore)
                    
                    Topscore(minindex,1) = maxpulse;
                    Topscore(minindex,2) = E;
                    Topscore(minindex,3) = score;
                    bestparam = [g0, alpha1/pi, alpha2/pi, alpha3/pi, alphap/pi];
                    Topvar(minindex,:) = bestparam;
                    
                    count = count + 1;
                    count2 = count2 + 1;
        
                    figure(count2), waterfall(t', z', abs(psisol))
                    title(['g0 = ' num2str(g0) ' alpha1 = ' num2str(alpha1/pi) '\pi alpha2 = ' num2str(alpha2/pi) '\pi alpha3 = ' num2str(alpha3/pi) '\pi alphap = ' num2str(alphap/pi) '\pi max pulse = ' num2str(maxpulse) ' E = ' num2str(E) ' Score = ' num2str(score)])

                    saveas(gcf, ['Found' num2str(count), '.png'])
   
                end
            end
        end
    end
    
    bestscore = max(Topscore(:,3))

end

Topscore
Topvar