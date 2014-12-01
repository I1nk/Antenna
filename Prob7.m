% Part A
%Clear the screen and memory
%clear 
clc

%Demensions of the metal plate
wy =  50e-3;
lx = 20e-3;
h = 2e-3;

%How to split up the metal plate
Nx = 15;
Ny = 25;
Nsq = Nx * Ny;
NsqBy2 = Nsq * 2;
Sarea = (wy/Ny)*(lx/Nx);

%Set up the array to hold the position of each 
%element in the metal plate
ele(NsqBy2).x = 0;
ele(NsqBy2).y = 0;
ele(NsqBy2).z = 0;

%Set up the equations for the loop for finding the postions
%of each element in the metal plate
index = 1;
inter_counter = 0;
outer_counter = 0;

%Work on the bottom plate
for outer = 1:Ny
   
   %Find the Y value
   tempy = wy;
   tempy = tempy + (2 * wy * outer_counter);
   tempy = tempy  / (2 * Ny);

   for inter = 1:Nx

      %Find the x value
      tempx = lx;
      tempx = tempx + (2 * lx * inter_counter);
      tempx = tempx / (2 * Nx);

      %Store it into the struct array
      ele(index).x = tempx;
      ele(index).y = tempy;
      ele(index).z = 0;
   
      %Increment the x value counter
      inter_counter = inter_counter + 1;
      index = index + 1;
   end
   
      %Reset the inter_counter 
      inter_counter = 0;
      %Increment the y value counter
      outer_counter = outer_counter + 1;
end


%Reset the counters for the top plate
outindex = 1;


%Work on the top plate
for outer = 1:Ny
   for inter = 1:Nx
      ele(index).z = h;
      ele(index).x = ele(outindex).x;
      ele(index).y = ele(outindex).y;
      outindex = outindex + 1;
      index = index + 1;
   end
end


%Get rid of the varaibles that is no longer needed for this MOM.
clear index inter inter_counter outer outer_counter outindex tempx tempy ans


%Set up the equations for the matrix
eq = zeros(NsqBy2,NsqBy2);
count = 1;
constantval = kii(ele(1).x, ele(1).y, ele(1).z,...
             Nx,Ny,wy,lx);

%Loop though the array making the system of equations
for outer = 1:NsqBy2
   for inter = 1:NsqBy2
      if(outer ~= inter)
         eq(outer,inter) = kij(ele(outer).x, ele(outer).y, ele(outer).z,...
                               ele(inter).x, ele(inter).y, ele(inter).z,...
                               Sarea);
      else
         eq(outer,inter) = constantval;
      end
   end
end

%Set up the voltage for solving the system of equations
Vtot = 2;
V1 = 1;
V = ones(NsqBy2,1);
V(Nsq+1:NsqBy2) = -V1;

%Solve the system of equations
Q = linsolve(eq, V);

%Find the C in F
q1 = sum((Q(1:Nsq))) * Sarea
C = q1 * 0.5

%Below is the equation for a parallel cap
actCap = ((50e-3)*(20e-3)*(8.853e-12)) / (2e-3)

