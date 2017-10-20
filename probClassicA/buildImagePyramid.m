function [pyramid, plevels] = buildImagePyramid(images, plevels, pscale, pmin, pmax)
% Construct image pyramid following the approach of Sun et al. [43]
    
    if isempty(plevels)
        N1 = 1 + floor(log(max(size(images,1),size(images,2))/pmax)/log(1/pscale));
        N2 = 1 + floor(log(min(size(images,1),size(images,2))/pmin)/log(1/pscale));
    	plevels = min(N1, N2);
        if plevels < 0
            plevels = 1;
        end
    end
    
    sigmaFilter = sqrt(1/pscale) / sqrt(2);
    filter = fspecial('gaussian', 2*round(1.5*sigmaFilter) + 1, sigmaFilter);
    pyramid = compute_image_pyramid(images, filter, plevels, pscale);

end

