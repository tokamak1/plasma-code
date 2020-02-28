function y = phi(n,x)
%phi Summary of this function goes here
%   FEM basic function

    switch n
        case -1
            y = [0, 0, 0, 0];
        case 0
            y = [1, x, x^2, x^3];
        case 1
            y = [0, 1, 2*x, 3*x^2];
        case 2
            y = [0, 0, 2, 6*x];
    end
    
end

