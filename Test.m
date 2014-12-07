clear
clc
close all

gap = 1e-3;
lamda = 0.5;
lens = 0.85;
l = lamda*lens;
N = (l/gap) + (mod((l/gap),2) == 0); % number of sections
E = zeros(1,N); % magnetic frill generator
%E(ceil(N/2)) = 100/2;
a = 0.0001;
b = 2.3*a;
vs = 100;
offset = l * -0.5;
delz = l / N;
zvalue = zeros(1,N);

for k = 1:N
    
    %find the mid point for the segment and offset it to -l/2 to l/2
    tempy = 2 * k;
    tempy = tempy - 1;
    tempy = tempy * l;
    tempy = tempy * 0.5;
    tempy = tempy / N;
    tempy = tempy + offset;
    zvalue(k) = tempy;
    
    %find r1 and r2
    r1 = tempy * tempy;
    r2 = r1;
    r1 = r1 + (a*a);
    r2 = r2 + (b*b);
    r1 = sqrt(r1);
    r2 = sqrt(r2);
    
    %      i
    %Find E
    %      k
    E(k) = 2 * log(b/a);
    E(k) = 1 / E(k);
    E(k) = E(k) * ( (exp(-1j*2*pi*r1/lamda) / r1) ...
        - (exp(-1j*2*pi*r2/lamda) / r2) );
    
end

% get current distribution (normal conditions)
[I,z,cnd] = pfield(l/lamda,.0001/lamda,E,'e','d');
for i = 1:N
    I(i) = abs(I(i));
    z(i) = z(i) * lamda;
end

figure;
plot(zvalue,I)
grid on
figure;
plot(zvalue,E)
grid on