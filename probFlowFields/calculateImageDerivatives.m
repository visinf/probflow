function [It,Ix,Iy,Ixx,Ixy,Iyx,Iyy,Ixt,Iyt] = calculateImageDerivatives(im1,im2,flow,derivFilter)
% Calculate first- and second-order image derivatives, cf. Revaud et al. [36]
    
    % Warp second image, get out-of-boundary mask
    [im2_w, outOfBoundary] = warpIm_original(im2, flow);
    
    % Calculate temporal difference
    It = im2_w - im1; 
    
    % Average images to obtain averaged gradients
    tempAvImg = 0.5 * im1 + 0.5 * im2_w;

    % Compute first-order spatial derivatives
    Ix = imfilter(tempAvImg,derivFilter, 'replicate');
    Iy = imfilter(tempAvImg,derivFilter', 'replicate');
    
    % Compute second-order derivatives
    Ixx = imfilter(Ix,derivFilter, 'replicate');
    Ixy = imfilter(Ix,derivFilter', 'replicate');
    Iyx = Ixy;
    Iyy = imfilter(Iy,derivFilter', 'replicate');
    Ixt = imfilter(It,derivFilter, 'replicate');
    Iyt = imfilter(It,derivFilter', 'replicate');
        
    % Set derivatives of out-of-boundary pixels to zero
    It(outOfBoundary) = 0;
    Ix(outOfBoundary) = 0;
    Iy(outOfBoundary) = 0;
    Ixx(outOfBoundary) = 0;
    Ixy(outOfBoundary) = 0;
    Iyx(outOfBoundary) = 0;
    Iyy(outOfBoundary) = 0;
    Ixt(outOfBoundary) = 0;
    Iyt(outOfBoundary) = 0;    
    
end

function a = rectify(a,b)
    a(a<1) = 1;
    a(a>b) = b;
end

function [im_w, outOfBound] = warpIm_original(im, flow)
% Backward-warping of image im with flow
    
    % Initialize help variables
    [n, m, c] = size(im);
    im_w = zeros(n,m,c);

    % Compute image grid (x,y) and shifted grid (xx,yy)
    [x, y] = meshgrid(1:m, 1:n);
    xx = x + flow(:,:,1);        
    yy = y + flow(:,:,2);
    
    % Detect out-of-boundary pixels
    outOfBound = (xx<1) | (yy<1) | (xx>m) | (yy>n);
    outOfBound = repmat(outOfBound, [1,1,3]);
    
    % Determine neighboring pixels for bilinear interpolation
    x1 = floor(xx);
    y1 = floor(yy);
    dx = xx-x1;
    dy = yy-y1;
    
    % Set out-of-boundary to nearest neighbor values
    x2 = rectify(x1+1, m);
    x1 = rectify(x1, m);
    y2 = rectify(y1+1, n);
    y1 = rectify(y1, n);
    
    % Compute linear indices for bilinear interpolation
    idx0 = sub2ind([n, m], y, x);
    idx1 = sub2ind([n, m], y1, x1);
    idx2 = sub2ind([n, m], y1, x2);
    idx3 = sub2ind([n, m], y2, x1);
    idx4 = sub2ind([n, m], y2, x2);
    
    % Perform bilinear interpolation per channel
    for i = 1:c
        tmp = im(:,:,i);
        tmp(idx0) = tmp(idx1).*(1-dx).*(1-dy) + tmp(idx2).*dx.*(1-dy) + tmp(idx3).*(1-dx).*dy + tmp(idx4).*dx.*dy;
        im_w(:,:,i) = tmp;
    end
end

