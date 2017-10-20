function [sig1, sig2] = updateSig(kd, ks, Ix2, Iy2, F1, F2, gsmD, gsmS, lamD, lamS, lamC)
% Perform update step of variances sig1 and sig2, cf. equation (34) of
% the supplemental material

    % Initialize help variables
    [n,m] = size(Ix2);
    npixels = n*m;
    F1t = F1';
    F2t = F2';
    
    % Compute contribution of data term
    termD = bsxfun(@rdivide,kd,2*gsmD.sigma);
    termD = sum(termD,1);
    termD = reshape(termD, [n,m]);
    
    % Compute contribution of smoothness terms 1 & 2
    termS1 = sum(bsxfun(@rdivide, ks(:,:,1), 2*gsmS.sigma), 1);
    termS1 = sum(abs(F1t) * spdiags(termS1(:), 0, npixels, npixels),2);
    termS1 = reshape(termS1,[n,m]);
    termS2 = sum(bsxfun(@rdivide, ks(:,:,2), 2*gsmS.sigma), 1);
    termS2 = sum(abs(F2t) * spdiags(termS2(:), 0, npixels, npixels),2);
    termS2 = reshape(termS2,[n,m]);
    
    % Update sig1
    denominator = 2 * lamS * (termS1 + termS2) + 2 * lamD * termD .* Ix2 + 2 * lamC;
    sig1 = bsxfun(@rdivide, 1, denominator);
    
    
    % Compute contribution of smoothness terms 3 & 4
    termS3 = sum(bsxfun(@rdivide, ks(:,:,3), 2*gsmS.sigma), 1);
    termS3 = sum(abs(F1t) * spdiags(termS3(:), 0, npixels, npixels),2);
    termS3 = reshape(termS3,[n,m]);
    termS4 = sum(bsxfun(@rdivide, ks(:,:,4), 2*gsmS.sigma), 1);
    termS4 = sum(abs(F2t) * spdiags(termS4(:), 0, npixels, npixels),2);
    termS4 = reshape(termS4,[n,m]);
    
    % Update sig2
    denominator = 2 * lamS * (termS3 + termS4) + 2 * lamD * termD .* Iy2 + 2 * lamC;
    sig2 = bsxfun(@rdivide, 1, denominator);

end

