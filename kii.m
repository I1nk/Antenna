function [out] = kii(x,y,z,Nx,Ny,wy,lx)

   ep = 8.854e-12;
   k = 2 * pi * ep;
   k = 1 / k;

   a = lx/Nx;
   b = wy/Ny;
   
   temp = a * a;
   temp =  temp + (b * b);

   temp = sqrt(temp);

   %Doing this so this value won't need to be found two times
   temp1 = temp;
   
   temp = temp + b;
   temp = temp / a;
   temp = log(temp);
   temp = temp * a;
   
  
   temp1 = temp1 + a;
   temp1 = temp1 / b;
   temp1 = log(temp1);
   temp1 = temp1 * b;
   
   out = temp1 + temp;
   out = out * k;

end 
