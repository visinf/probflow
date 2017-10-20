function [lambda] = calculateLamS(im, kappa, derivFilter)
% Compute spatially varying trade-off parameter lambda based on image
% derivatives, cf. Revaud et al. [36]

    % Get image luminance
    luminance = (0.299 * im(:,:,1) + 0.587 * im(:,:,2) + 0.114 * im(:,:,3))/255;
    
    % Compute image derivatives
    Ix = imfilter(luminance, derivFilter,  'conv', 'replicate', 'same');
    Iy = imfilter(luminance, derivFilter',  'conv', 'replicate', 'same');

    % Calculate trade-off parameter lambda
    expLam = - kappa * sqrt(Ix .* Ix + Iy .* Iy);
    lambda = exp(expLam);
    
end

