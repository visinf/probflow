function [mu_new, sig1_new, sig2_new] = upsampleVariables(mu, sig1, sig2, newDim)
% Upsample the given variables to the new dimension newDim using bilinear
% interpolation

    mu_new = zeros([newDim 2]);
    imDim = size(mu(:,:,1));

    factor1 = newDim(2) / imDim(2);
    mu_new(:,:,1)  = imresize(mu(:,:,1), newDim, 'bilinear') * factor1;
    factor2 = newDim(1) / imDim(1);
    mu_new(:,:,2) = imresize(mu(:,:,2), newDim, 'bilinear') * factor2;

    sig1_new = imresize(sig1, newDim, 'bilinear') * factor1;

    sig2_new = imresize(sig2, newDim, 'bilinear') * factor2;
    

end
