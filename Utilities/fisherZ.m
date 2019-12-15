function z = fisherZ(r)

%Need to round to 3 sig digits because fisherZ(0.9999) > 3.8002 (the max
%value allowed by the transform). 
r = round(r,3);
r(r==1)  =  .999;
r(r==-1) = -.999;
z        = .5 .* log((1+r) ./ (1-r));