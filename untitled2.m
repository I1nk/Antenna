clear
clc

l =10;
offset = -l*.5;
N = 5;
for k = 1:N    
%find the mid point for the segment
    tempy = 2 * k;
    tempy = tempy - 1;
    tempy = tempy * l;
    tempy = tempy * 0.5;
    tempy = tempy / N;
    tempy = tempy + offset
end