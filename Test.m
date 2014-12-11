clear
clc

vs = 100;
gap = 1e-3;
lamda = 0.5;
lens = .85;
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
    E(k) = 2 * log(2.3); %b/a = 2.3
    E(k) = vs / E(k);
    E(k) = -E(k) * ( (exp(-1j*2*pi*r1/lamda) / r1) ...
        - (exp(-1j*2*pi*r2/lamda) / r2) );
    
end

%%
% get current distribution (normal conditions)
[I,z,cnd] = pfield(lens,a/lamda,E,'e','d');

%Find the magnitude of the current
Ia = abs(I);
%Set up the impedance 
z = z * lamda;

%Find the impedance of the antenna
for k = 1:N
   if(I(k) == 0)
      imp(k) = 0;
   else
      imp(k) = E(k) * abs(zvalue(k));
      imp(k) = imp(k) / I(k);
   end
end

%for matlab only
%imp = E .* z;
%imp = imp ./ I;

%find the sum of the impedance of the antenna

imp_sum = sum(imp);

fprintf('Z = %f %fi ohms\n',real(imp_sum),imag(imp_sum));

%******************************************************
%     COMPUTATION OF THE INPUT IMPEDANCE
%******************************************************
zin=100/I(ceil(N/2))

beta = 2 * pi / lamda;
eta = 377;
etm=zeros(361);
etmm=etm(1,1:361);
hl = l/2;
dz = delz;
nm = N;
%******************************************************************
%     COMPUTATION OF AMPLITUDE RADIATION PATTERN OF THE ANTENNA
%******************************************************************
for i=1:181
    theta=(i-1.0)*pi/180;
    cth=cos(theta);
    sth=sin(theta);
      if abs(cth)<0.001
         ft=1;
      else
         ft=sin(beta*dz*cth*0.5)/(beta*dz*cth*0.5);
      end
    crt=0;
         for m=1:nm
             zm=hl-(m-0.5)*dz;
             crt=crt+exp(j*beta*zm*cth)*ft*I(m)*dz;
         end
    ptt=abs(crt)*sth*sth*eta*0.5;
    etmm(i)=ptt;
 end
amax=etmm(1);
for i=2:181
    if etmm(i)>=amax
       amax=etmm(i);
    end
end
for i=1:181
    ptt=etmm(i)/amax;
      if ptt<=0.00001
         ptt=0.00001;
         etmm(i)=20*log10(ptt);
      end
      etmm(i)=20*log10(ptt);
end


figure;
i=[1:181];
xi=i-1;
% Polar Plot
etmm1=[etmm(1:181),fliplr(etmm(1:180))];
q=polar_dB([0:360],etmm1,-40,0,4,'-');
set(q,'linewidth',1.5);
title('RADIATION PATTERN   vs   OBSERVATION ANGLE')


%%
%plot the information

figure;
plot(z,Ia)
xlim([-l/2 l/2])
grid on

%For testing only and not going to use this plot to turn in with the report.
figure;
plot(z,abs(imp))
xlim([-l/2 l/2])
grid on




