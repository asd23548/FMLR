function loc = index2loc(index,imagesize)
if length(imagesize) == 3
    d1 = imagesize(1);
    d2 = imagesize(2);
    d3 = imagesize(3);
    
    i = ceil(index/(d2*d3));
    R = index - (i-1) * d1*d2;
    j = ceil(R/d3);
    k = R-(j-1)*d3;
    loc = [i,j,k]';
end

if length(imagesize) == 2
    d1 = imagesize(1);
    d2 = imagesize(2);
    
    i = ceil(index/(d1));
    j = index - (i-1) * d2;
    loc = [i,j]';
end
