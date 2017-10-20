function [muHat] = updateMuHat(mu, kHat, neighbors, gsmN, lamC, lamN)
% Perform update step of flow estimates muHat, cf. equation (33) of the 
% supplemental material

    [n,m,~] = size(mu);
    
    kHat = permute(sum(bsxfun(@rdivide, kHat, 2*gsmN.sigma),1),[2,3,1]);
    
    numerator = lamC * mu + lamN * reshape(sum(kHat .* neighbors,1), [n m]);
    denominator = lamC + lamN * reshape(sum(kHat,1), [n m]);
    
    muHat = numerator ./ denominator;

end

