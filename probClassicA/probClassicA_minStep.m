function [mu, sig1, sig2, kd, ks] = probClassicA_minStep(mu0, muHat, sig1, sig2, It, Ix, Iy, F1, F2, lamD, lamS, lamC, gsmD, gsmS, nIter)
% Perform one minimization step for (mu,sig,k)

    % Initialize variables
    [n,m,~] = size(mu0);
    mu = mu0;
    
    % Generate matrix corresponding to coupling term E_c
    Ac = spdiags(2 * lamC * ones(2*n*m,1), 0, 2*n*m, 2*n*m);
    bc = 2 * lamC * muHat(:);
    
    % Calculate factors of image derivatives
    Ix2 = Ix.^2;
    Ixy = Ix .* Iy;
    Iy2 = Iy.^2;
    Itx = It .* Ix;
    Ity = It .* Iy;
    It2 = It.^2;
    
    % Initial update for k and sig
    [kd,ks] = updateK(mu, mu0, sig1, sig2, Ix2, Ixy, Iy2, It2, Itx, Ity, F1, F2 , gsmD, gsmS, lamD, lamS);
    [sig1, sig2] = updateSig(kd, ks, Ix2, Iy2, F1, F2, gsmD, gsmS, lamD, lamS, lamC);
    
    % Alternating estimation of mu, sig and k
    for i = 1:nIter
        mu = updateMu(mu0, kd, ks, Ix2, Ixy, Iy2, Itx, Ity, F1, F2, gsmD, gsmS, lamD, lamS, Ac, bc);
        [kd,ks] = updateK(mu, mu0, sig1, sig2, Ix2, Ixy, Iy2, It2, Itx, Ity, F1, F2, gsmD, gsmS, lamD, lamS);
        [sig1, sig2] = updateSig(kd, ks, Ix2, Iy2, F1, F2, gsmD, gsmS, lamD, lamS, lamC);
    end
    
    % Restrict flow increment to an absolute value of 1
    dmu = mu - mu0;
    dmu(dmu > 1)  = 1;
    dmu(dmu < -1) = -1;
    mu = mu0 + dmu;
    
    % Re-estimate k and sig according to bounded mu
    [kd, ks] = updateK(mu, mu0, sig1, sig2, Ix2, Ixy, Iy2, It2, Itx, Ity, F1, F2, gsmD, gsmS, lamD, lamS);
    [sig1, sig2] = updateSig(kd, ks, Ix2, Iy2, F1, F2, gsmD, gsmS, lamD, lamS, lamC);
end

