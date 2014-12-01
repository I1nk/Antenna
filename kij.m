function [out] = kij(x,y,z,xp,yp,zp,s)

   %surface area
   
   ep = 8.854e-12;
   
   %dem constant
   k = 4 * pi * ep;

   %Find the differance sq for x y and z
   dx = x - xp;
   dx = dx * dx;
   
   dy = y - yp;
   dy = dy * dy;
   
   dz = z - zp;
   dz = dz * dz;

   
   sumsq = dx + dy + dz;

   sumes = sqrt(sumsq);

   out = k * sumes;
   out = s / out;

end
