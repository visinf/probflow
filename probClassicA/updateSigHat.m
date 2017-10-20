function [sigHat] = updateSigHat(kHat, gsmN, lamC, lamN, n, m)
% Perform update step of variance estimate sigHat, cf. equation (35) of the 
% supplemental material

    coeffK = 2 * lamN * sum(sum(bsxfun(@rdivide, kHat, 2*gsmN.sigma), 1), 2);
    coeffK = reshape(coeffK, [n m]);

    coeffSig = coeffK + 2 * lamC;
    sigHat = bsxfun(@rdivide, 1, coeffSig);

end

