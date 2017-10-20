function [kd,ks] = updateK(mu, mu0, sig1, sig2, Ix2, Ixy, Iy2, It2, Itx, Ity, F1, F2, gsmD, gsmS, lamD, lamS)
% Perform update step of latent variables k, cf. equation (36) of the
% supplemental material
    
    % Initialize help variables
    [n,m,~] = size(mu);
    Ls = size(gsmS.sigma,1);
    ks = zeros(Ls,n*m,4);
    mu1 = mu(:,:,1);
    mu2 = mu(:,:,2);
    dmu1 = mu(:,:,1) - mu0(:,:,1);
    dmu2 = mu(:,:,2) - mu0(:,:,2);
    
    % Compute data term
    termD = It2 + 2 * Itx .* dmu1 + 2 * Ity .* dmu2 + Ix2 .* sig1 + Iy2 .* sig2 ...
        + Ix2 .* dmu1.^2 + 2 * Ixy .* dmu1 .* dmu2 + Iy2 .* dmu2.^2;
    
    % Determine exponential of k update
    expD1 = bsxfun(@rdivide,termD(:)',2*gsmD.sigma);
    expD2 = log(sqrt(gsmD.sigma)) - log(gsmD.pi);
    expD = - lamD * bsxfun(@plus,expD1,expD2);
    
    % Numerically stable computation of exp(expD)/sum(exp(expD))
    maxD = max(expD,[],1);
    expD = bsxfun(@minus,expD,maxD);
    expD = exp(expD);
    zD = sum(expD,1);
    kd = bsxfun(@rdivide,expD,zD);

    
    % Compute smoothness term 1
    termS = (F1 * mu1(:)).^2 + abs(F1) * sig1(:);
    
    % Determine exponential of k update
    expS1 = bsxfun(@rdivide,termS(:)',2*gsmS.sigma);
    expS2 = log(sqrt(gsmS.sigma)) - log(gsmS.pi);
    expS = - lamS * bsxfun(@plus,expS1,expS2);
    
    % Numerically stable computation of exp(expS)/sum(exp(expS))
    maxS = max(expS,[],1);
    expS = bsxfun(@minus,expS,maxS);
    expS = exp(expS);
    zS = sum(expS,1);
    ks(:,:,1) = bsxfun(@rdivide,expS,zS);
    
    
    % Compute smoothness term 2
    termS = (F2 * mu1(:)).^2 + abs(F2) * sig1(:);
  
    % Determine exponential of k update
    expS1 = bsxfun(@rdivide,termS(:)',2*gsmS.sigma);
    expS = - lamS * bsxfun(@plus,expS1,expS2);
    
    % Numerically stable computation of exp(expS)/sum(exp(expS))
    maxS = max(expS,[],1);
    expS = bsxfun(@minus,expS,maxS);
    expS = exp(expS);
    zS = sum(expS,1);
    ks(:,:,2) = bsxfun(@rdivide,expS,zS);
    

    % Compute smoothness term 3
    termS = (F1 * mu2(:)).^2 + abs(F1) * sig2(:);

    % Determine exponential of k update
    expS1 = bsxfun(@rdivide,termS(:)',2*gsmS.sigma);
    expS = - lamS * bsxfun(@plus,expS1,expS2);
    
    % Numerically stable computation of exp(expS)/sum(exp(expS))
    maxS = max(expS,[],1);
    expS = bsxfun(@minus,expS,maxS);
    expS = exp(expS);
    zS = sum(expS,1);
    ks(:,:,3) = bsxfun(@rdivide,expS,zS);
    

    % Compute smoothness term 4
    termS = (F2 * mu2(:)).^2 + abs(F2) * sig2(:);
   
    % Determine exponential of k update
    expS1 = bsxfun(@rdivide,termS(:)',2*gsmS.sigma);
    expS = - lamS * bsxfun(@plus,expS1,expS2);
    
    % Numerically stable computation of exp(expS)/sum(exp(expS))
    maxS = max(expS,[],1);
    expS = bsxfun(@minus,expS,maxS);
    expS = exp(expS);
    zS = sum(expS,1);
    ks(:,:,4) = bsxfun(@rdivide,expS,zS);
 
end

