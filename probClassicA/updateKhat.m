function [kHat] = updateKhat(muHat, sigHat, neighbors, neighborsSig, gsmN, lamN)
% Perform update step of latent variables kHat, cf. equation (36) of the
% supplemental material
    
    % Substract current estimates from neighborhood values 
    neighbors = bsxfun(@minus, muHat(:)', neighbors);
    neighborsSig = bsxfun(@plus, sigHat(:)', neighborsSig);
    
    % Compute nonlocal term
    logSig_logPi = log(sqrt(gsmN.sigma)) - log(gsmN.pi);
    neighbors = neighbors.^2 + neighborsSig;
    kHat = -lamN * bsxfun(@plus, bsxfun(@rdivide, permute(neighbors, [3,1,2]), 2*gsmN.sigma), logSig_logPi);

    % Numerically stable computation of exp(expKhat)/sum(exp(expKhat))
    maxData = max(kHat,[],1);
    kHat = bsxfun(@minus,kHat,maxData);
    kHat = exp(kHat);
    zData = sum(kHat,1);
    kHat = bsxfun(@rdivide,kHat,zData);
end

