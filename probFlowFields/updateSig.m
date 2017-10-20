function [sig1, sig2] = updateSig(kd1, kd2, ks, Ix2, Iy2, Ixx2, Ixy2, Iyx2, Iyy2, gsmD1, gsmD2, gsmS, energyOpt, n, m)
% Perform update step of variances sig1 and sig2, cf. equation (34) of
% the supplemental material

    % Initialize help variables
    lamD1 = energyOpt.lamD1;
    lamD2 = energyOpt.lamD2;
    lamS = energyOpt.lamS;
    filter_smooth = energyOpt.filter_smooth;
    lamS1 = 0.5 * imfilter(lamS,abs(filter_smooth{1}),'conv','replicate');
    lamS2 = 0.5 * imfilter(lamS,abs(filter_smooth{2}),'conv','replicate');
    
    % Compute contribution of data term 1
    termD1 = bsxfun(@rdivide,kd1,2*gsmD1.sigma);
    termD1 = sum(termD1,1);
    termD1 = reshape(termD1, [n,m]);
    
    % Compute contribution of data term 2
    termD2 = bsxfun(@rdivide,kd2,2*gsmD2.sigma);
    termD2 = sum(termD2,1);
    termD2 = reshape(termD2, [n,m]);
    
    % Compute contribution of smoothness terms 1 & 2
    tmp = lamS1 .* reshape(sum(bsxfun(@rdivide, ks(:,:,1), 2*gsmS.sigma), 1),[n,m,1]);
    termS1 = zeros(n,m);
    termS1(1:end-1,:) = tmp(1:end-1,:);
    termS1(2:end,:) = termS1(2:end,:) + tmp(1:end-1,:);
    tmp = lamS2 .* reshape(sum(bsxfun(@rdivide, ks(:,:,2), 2*gsmS.sigma), 1),[n,m,1]);
    termS2 = zeros(n,m);
    termS2(:,1:end-1) = tmp(:,1:end-1);
    termS2(:,2:end) = termS2(:,2:end) + tmp(:,1:end-1);
    
    % Update sig1
    denominator = 2 * (termS1 + termS2) + 2 * lamD1 * termD1 .* Ix2 + 2 * lamD2 * termD2 .* (Ixx2 + Iyx2);
    sig1 = bsxfun(@rdivide, 1, denominator);
    
    % Compute contribution of smoothness terms 1 & 2
    tmp = lamS1 .* reshape(sum(bsxfun(@rdivide, ks(:,:,1), 2*gsmS.sigma), 1),[n,m,1]);
    termS1 = zeros(n,m);
    termS1(1:end-1,:) = tmp(1:end-1,:);
    termS1(2:end,:) = termS1(2:end,:) + tmp(1:end-1,:);
    tmp = lamS2 .* reshape(sum(bsxfun(@rdivide, ks(:,:,2), 2*gsmS.sigma), 1),[n,m,1]);
    termS2 = zeros(n,m);
    termS2(:,1:end-1) = tmp(:,1:end-1);
    termS2(:,2:end) = termS2(:,2:end) + tmp(:,1:end-1);
    
    % Update sig2
    denominator = 2 * (termS1 + termS2) + 2 * lamD1 * termD1 .* Iy2 + 2 * lamD2 * termD2 .* (Ixy2 + Iyy2);
    sig2 = bsxfun(@rdivide, 1, denominator);

end

