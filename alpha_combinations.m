
close all;
clear all;
clc;

count = 1;

dalpha = 0.1;

for i = 0:dalpha:2
    for j = 0:dalpha:2
        for k = 0:dalpha:2
            for m = 0:dalpha:2
                
                alpha(count,:) = [i j k m];
                count = count + 1;
                
            end
        end
    end
end

alpha



