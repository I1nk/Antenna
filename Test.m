clear
clc

vs = 100;
gap = 1e-3;
lamda = 0.5;
lens = 0.85;
l = lamda*lens;
N = ceil(l/gap) + (mod(ceil(l/gap),2) == 0);
zvalue = zeros(1,N);
E = zeros(1,N);
imp = zeros(1,N);
a = .0001;
b = 2.3*a;
offset = l * -0.5;
delz = l / N;
%This is for the delta gap generator
%E(ceil(N/2)) = vs/gap;

%%
%Set up the magnetic frill generator
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
    E(k) = vs / E(k);
    E(k) = -E(k) * ( (exp(-1j*2*pi*r1/lamda) / r1) ...
        - (exp(-1j*2*pi*r2/lamda) / r2) );
    
end

%%
% get current distribution (normal conditions)
[I,z,cnd] = pfield(l/lamda,.0001/lamda,E,'e','d');

%Find the magnitude of the current
Ia = abs(I);
%Set up the impedance 
z = z * lamda;

%Find the impedance of the antenna
for k = 1:N
   if(I(k) == 0)
      imp(k) = 0;
   else
      imp(k) = E(k) * zvalue(k);
      imp(k) = imp(k) / I(k);
   end
end

%for matlab only
%imp = E .* z;
%imp = imp ./ I;

%%
%plot the information

figure;
plot(z,Ia)
xlim([-l/2 l/2])
grid on


figure;
plot(z,abs(imp))
xlim([-l/2 l/2])
grid on

