function [] = Project()
    % Pockington with magnetic frill generator
    start = 1;
%% Set up the varaibles for the Main project
%while loop here to srink the script down when it is not needed to be seen
while(start == 1)
    %Constants for the antenna
    %Thickness of the antenna (radius wise) in meters
    a = 1e-4;
    %Wavelength of the antenna
    lamda = 0.5;
    %Voltage in the gap in V
    vo = 100;
    %gap size in m
    gap = 1; 
    %Length of the antenna
    l = 0.85 * lamda; 
    %Length of the antenna with gap
    ltot = l + gap; 
    lp = l * 0.5;
    beta = 2 * pi;

    %Demensions of the antenna
    radius = a;
    radius_sq = radius * radius;

    %How to split up the antenna
    %Base the number of elements in the entrie antenna from the number of 
    %elements in the gap
    gapsize = 3;
    %The number of elements based on the gap number of elements
    %The antenna is split into 3 sections with the same number of elements
    NumberOfElements = 3 * gapsize;

    %Set up the array to hold the position of each 
    %element in the antenna
    ele(NumberOfElements).x = 0;
    ele(NumberOfElements).y = 0;
    ele(NumberOfElements).z = 0;
    start = 0;
end

%reset the while loop varaible
start = 1;

%%  Begin making the antenna array
%while loop here to srink the script down when it is not needed to be seen
while(start == 1)
    %Place the base of the antenna at the orign
    %Set up the equations for the loop for finding the postions
    %of each element in the antenna
    index = 1;
    
    %Build the bottom the antenna
    [ele, index] = build(ele, index, gapsize, lp, radius, 0);
    %build the gap of the antenna 
    [ele, index] = build(ele, index, gapsize, lp, radius, lp);
    %build the top of the antenna
    [ele, index] = build(ele, index, gapsize, lp, radius, (lp + gap));
    start = 0;
end

%reset the while loop varaible
start = 1;

%%  Build the equation matrix    


    
end

%Default builder to put together the antenna 
function [out, indexo] = build(ele, index, gapsize, l, radius, offset)

    indexo = index;
    out = ele;
    
    %Work on the bottom plate
    for for_index = 1:gapsize

       %Find the Y value
       tempy = 2 * for_index;
       tempy = tempy - 1;
       tempy = tempy * l;
       tempy = tempy * 0.5;
       tempy = tempy / gapsize;
       tempy = tempy + offset;

       %Find the x value
       tempx = radius;

       %Store it into the struct array
       out(indexo).x = tempx;
       out(indexo).y = tempy;
       out(indexo).z = 0;
       indexo = indexo + 1;

    end
    
end

