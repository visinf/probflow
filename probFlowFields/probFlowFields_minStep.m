function [mu, sig1, sig2, kd1, kd2, ks] = probFlowFields_minStep(mu0, sig1, sig2, It, Ix, Iy, Ixx, Ixy, Iyx, Iyy, Ixt, Iyt, gsmD1, gsmD2, gsmS, nIter, sorSolver, energyOpt)
% Perform one minimization step for (mu,sig,k)
    
    % Initialize variables
    [n,m,~] = size(mu0);
    mu = mu0;
    
    % Normalize image derivatives
    [Ix2,IxIy,Iy2,IxIt,IyIt,It2,Ixx2,IxxIxy,Ixy2,IxxIxt,IxyIxt,Ixt2,Iyx2,IyxIyy,Iyy2,IyxIyt,IyyIyt,Iyt2] ...
        = normalizeDerivatives(It, Ix, Iy, Ixx, Ixy, Iyx, Iyy, Ixt, Iyt, energyOpt);
     
    % Initial update for k and sig
    [kd1,kd2,ks] = updateK(mu, mu0, sig1, sig2, Ix2, IxIy, Iy2, It2, Ixx2, IxxIxy, Ixy2, IxxIxt, IxyIxt, Ixt2, IxIt, IyIt, Iyx2, IyxIyy, Iyy2, IyxIyt, IyyIyt, Iyt2, gsmD1, gsmD2, gsmS, energyOpt);
    [sig1, sig2] = updateSig(kd1, kd2, ks, Ix2, Iy2, Ixx2, Ixy2, Iyx2, Iyy2, gsmD1, gsmD2, gsmS, energyOpt, n, m);

    % Alternating estimation of mu, sig and k
    for i = 1:nIter
        mu = updateMu(mu0, kd1, kd2, ks, Ix2, IxIy, Iy2, IxIt, IyIt, Ixx2, IxxIxy, Ixy2, IxxIxt, IxyIxt, Iyx2, IyxIyy, Iyy2, IyxIyt, IyyIyt, gsmD1, gsmD2, gsmS, sorSolver, energyOpt);
        [kd1,kd2,ks] = updateK(mu, mu0, sig1, sig2, Ix2, IxIy, Iy2, It2, Ixx2, IxxIxy, Ixy2, IxxIxt, IxyIxt, Ixt2, IxIt, IyIt, Iyx2, IyxIyy, Iyy2, IyxIyt, IyyIyt, Iyt2, gsmD1, gsmD2, gsmS, energyOpt);
        [sig1, sig2] = updateSig(kd1, kd2, ks, Ix2, Iy2, Ixx2, Ixy2, Iyx2, Iyy2, gsmD1, gsmD2, gsmS, energyOpt, n, m);
    end
    
end

