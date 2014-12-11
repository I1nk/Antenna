clear
clc


%This is going to be turned into a function to make running this file easier
%for the final report.

vs = 100;
gap = 1e-3;
lamda = 0.5;
lens = 0.85;
l = lamda*lens;
N = ceil(l/gap) + (mod(ceil(l/gap),2) == 0);
E = zeros(1,N);
imp = zeros(1,N);
a = .0001;
b = 2.3*a;
offset = l * -0.5;
delz = l / N;
%This is for the delta gap generator
E(ceil(N/2)) = vs/gap;


%%
% get current distribution (normal conditions)
[I,zvalue,cnd] = hfield(lens,a/lamda,E,'e','d');

%Find the magnitude of the current
Ia = abs(I);
%Set up the impedance 
zvalue = zvalue * lamda;

%Find the impedance of the antenna
for k = 1:N
   if(I(k) == 0)
      imp(k) = 0;
   else
      imp(k) = E(k) * abs(zvalue(k));
      imp(k) = imp(k) / I(k);
   end
end

%%
%plot the information

figure;
plot(zvalue,Ia)
xlim([-l/2 l/2])
grid on

%For testing only and not going to use this plot to turn in with the report.


