function [kd1,kd2,ks] = updateK(mu, mu0, sig1, sig2, Ix2, IxIy, Iy2, It2, Ixx2, IxxIxy, Ixy2, IxxIxt, IxyIxt, Ixt2, IxIt, IyIt, ...
    Iyx2, IyxIyy, Iyy2, IyxIyt, IyyIyt, Iyt2, gsmD1, gsmD2, gsmS, energyOpt)
% Perform update step of latent variables k, cf. equation (36) of the
% supplemental material

    % Initialize help variables
    [n,m,~] = size(mu);
    lamD1 = energyOpt.lamD1;
    lamD2 = energyOpt.lamD2;
    lamS = energyOpt.lamS;
    filter_smooth = energyOpt.filter_smooth;
    dmu1 = mu(:,:,1) - mu0(:,:,1);
    dmu2 = mu(:,:,2) - mu0(:,:,2);
    Ls = size(gsmS.sigma,1);
    ks = zeros(Ls,n*m,2);

    % Compute data term 1
    termD = It2 + 2 * IxIt .* dmu1 + 2 * IyIt .* dmu2 + Ix2 .* sig1 + Iy2 .* sig2 ...
        + Ix2 .* dmu1.^2 + 2 * IxIy .* dmu1 .* dmu2 + Iy2 .* dmu2.^2;
    
    % Determine exponential of k update
    expD1 = bsxfun(@rdivide,termD(:)',2*gsmD1.sigma);
    expD2 = log(sqrt(gsmD1.sigma)) - log(gsmD1.pi);
    expD = - lamD1 * bsxfun(@plus,expD1,expD2);
    
    % Numerically stable computation of exp(expD)/sum(exp(expD))
    maxD = max(expD,[],1);
    expD = bsxfun(@minus,expD,maxD);
    expD = exp(expD);
    zD = sum(expD,1);
    kd1 = bsxfun(@rdivide,expD,zD);
    
    
    % Compute data term 2
    termD = Ixt2 + 2 * IxxIxt .* dmu1 + 2 * IxyIxt .* dmu2 + Ixx2 .* sig1 + Ixy2 .* sig2 ...
        + Ixx2 .* dmu1.^2 + 2 * IxxIxy .* dmu1 .* dmu2 + Ixy2 .* dmu2.^2;
    
    termD = termD + Iyt2 + 2 * IyxIyt .* dmu1 + 2 * IyyIyt .* dmu2 + Iyx2 .* sig1 + Iyy2 .* sig2 ...
        + Iyx2 .* dmu1.^2 + 2 * IyxIyy .* dmu1 .* dmu2 + Iyy2 .* dmu2.^2;
    
    expD1 = bsxfun(@rdivide,termD(:)',2*gsmD2.sigma);
    expD2 = log(sqrt(gsmD2.sigma)) - log(gsmD2.pi);
    expD = - lamD2 * bsxfun(@plus,expD1,expD2);
    maxD = max(expD,[],1);
    expD = bsxfun(@minus,expD,maxD);
    expD = exp(expD);
    zD = sum(expD,1);
    kd2 = bsxfun(@rdivide,expD,zD);
    
    
    % Compute smoothness term 1
    filt = filter_smooth{1};
    lamS_tmp = 0.5 * imfilter(lamS,abs(filt),'conv','replicate');
    smoothTerm1 = imfilter(mu(:,:,1),filt,'conv','replicate');
    smoothTerm2 = imfilter(mu(:,:,2),filt,'conv','replicate');
    sig1_tmp = imfilter(sig1,abs(filt),'conv','replicate');
    sig2_tmp = imfilter(sig2,abs(filt),'conv','replicate');
    termS = smoothTerm1.^2 + smoothTerm2.^2 + sig1_tmp + sig2_tmp;
    
    % Determine exponential of k update
    expS1 = bsxfun(@rdivide,termS(:)',2*gsmS.sigma);
    expS2 = log(sqrt(gsmS.sigma)) - log(gsmS.pi);
    expS = bsxfun(@plus,expS1,expS2);
    expS = bsxfun(@times, - lamS_tmp(:)', expS);
    
    % Numerically stable computation of exp(expS)/sum(exp(expS))
    maxS = max(expS,[],1);
    expS = bsxfun(@minus,expS,maxS);
    expS = exp(expS);
    zS = sum(expS,1);
    ks(:,:,1) = bsxfun(@rdivide,expS,zS);

    
    % Compute smoothness term 2
    filt = filter_smooth{2};
    lamS_tmp = 0.5 * imfilter(lamS,abs(filt),'conv','replicate');
    smoothTerm1 = imfilter(mu(:,:,1),filt,'conv','replicate');
    smoothTerm2 = imfilter(mu(:,:,2),filt,'conv','replicate');
    sig1_tmp = imfilter(sig1,abs(filt),'conv','replicate');
    sig2_tmp = imfilter(sig2,abs(filt),'conv','replicate');
    termS = smoothTerm1.^2 + smoothTerm2.^2 + sig1_tmp + sig2_tmp;
    
    % Determine exponential of k update
    expS1 = bsxfun(@rdivide,termS(:)',2*gsmS.sigma);
    expS = bsxfun(@plus,expS1,expS2);
    expS = bsxfun(@times, - lamS_tmp(:)', expS);
    
    % Numerically stable computation of exp(expS)/sum(exp(expS))
    maxS = max(expS,[],1);
    expS = bsxfun(@minus,expS,maxS);
    expS = exp(expS);
    zS = sum(expS,1);
    ks(:,:,2) = bsxfun(@rdivide,expS,zS);
 
end

