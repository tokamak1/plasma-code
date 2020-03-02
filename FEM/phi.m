function y = phi(n,x)
%phi Summary of this function goes here
%   FEM basic function

    switch n
        case -1
            y = [x.^1; x.^2/2; x.^3/3; x.^4/4];
        case 0
            y = [x.^0; x; x.^2; x.^3];
        case 1
            y = [0*x; x.^0; 2*x; 3*x.^2];
        case 2
            y = [0*x; 0*x; 2*x.^0; 6*x];
    end
    
end

