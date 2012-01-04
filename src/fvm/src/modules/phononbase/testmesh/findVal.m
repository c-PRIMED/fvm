function [ y ] = findVal( x1,x2,y1,y2,x)
%findVal Finds linearly interpolated value

slope=(y2-y1)/(x2-x1);
y=slope*(x-x1)+y1;
end

