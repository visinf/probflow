function [mu] = updateMu(mu0, kd, ks, Ix2, Ixy, Iy2, Itx, Ity, F1, F2, gsmD, gsmS, lamD, lamS, Ac, bc)
% Perform update step of flow estimates mu, i.e. solve equation system 
% A*x = b as given in equation (31) of the supplemental material

    [A,b] = setUpEquationSystem(mu0, kd, ks, Ix2, Ixy, Iy2, Itx, Ity, F1, F2, gsmD, gsmS, lamD, lamS, Ac, bc);
    
    mu = full(A\b); 

    mu = reshape(mu, size(mu0));
    
end

