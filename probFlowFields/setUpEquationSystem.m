function [A, b] = setUpEquationSystem(mu0, kd1, kd2, ks, Ix2, IxIy, Iy2, IxIt, IyIt, Ixx2, IxxIxy, Ixy2, IxxIxt, IxyIxt, Iyx2, IyxIyy, Iyy2, IyxIyt, IyyIyt, gsmD1, gsmD2, gsmS, energyOpt)
% Determine equation system that needs to be solved in order to update flow
% estimate mu

    [n,m,~] = size(mu0);
    npixels = n * m;
    
    lamD1 = energyOpt.lamD1;
    lamD2 = energyOpt.lamD2;
    lamS = energyOpt.lamS;
    filter_smooth = energyOpt.filter_smooth;
    
     % Determine contribution of data term 1
    termD = bsxfun(@rdivide,kd1,2*gsmD1.sigma);
    termD = sum(termD,1);
    termD = 2 * lamD1 * termD(:);

    tmp = termD.*Ix2(:);
    dxx = spdiags(tmp, 0, npixels, npixels);
    tmp = termD.*Iy2(:);
    dyy = spdiags(tmp, 0, npixels, npixels);
    tmp = termD.*IxIy(:);
    dxy = spdiags(tmp, 0, npixels, npixels);
    
    Ad = [dxx dxy; dxy dyy];
    bd = - [termD .* IxIt(:); termD .* IyIt(:)];
    
    % Determine contribution of data term 2
    termD = bsxfun(@rdivide, kd2, 2*gsmD2.sigma);
    termD = sum(termD, 1);
    termD = 2 * lamD2 * termD(:);
    
    tmp = termD .* (Ixx2(:) + Iyx2(:));
    dxx = spdiags(tmp, 0, npixels, npixels);
    tmp = termD .* (Ixy2(:) + Iyy2(:));
    dyy = spdiags(tmp, 0, npixels, npixels);
    tmp = termD .* (IxxIxy(:) + IyxIyy(:));
    dxy = spdiags(tmp, 0, npixels, npixels);

    % Compute data matrix Ad and vector bd
    Ad = Ad + [dxx dxy; dxy dyy];
    bd = bd + Ad * mu0(:) - [termD .* (IxxIxt(:) + IyxIyt(:)); termD .* (IxyIxt(:) + IyyIyt(:))];
    
    
    % Determine contribution of smoothness terms
    termS = bsxfun(@rdivide, ks, 2*gsmS.sigma);
    termS = sum(termS, 1);
    
    filt = filter_smooth{1};
    lamS1 = 0.5 * imfilter(lamS, abs(filt), 'conv', 'replicate');
    F1 = make_convn_mat(filt, [n m], 'valid', 'same');
    s1 = termS(:,:,1);
    s1 = spdiags(2 * lamS1(:) .* s1(:), 0, npixels, npixels);
    
    filt = filter_smooth{2};
    lamS2 = 0.5 * imfilter(lamS, abs(filt), 'conv', 'replicate');
    F2 = make_convn_mat(filt, [n m], 'valid', 'same');
    s2 = termS(:,:,2);
    s2 = spdiags(2 * lamS2(:) .* s2(:), 0, npixels, npixels);

    % Compute smoothness matrix As and vector bs
    As = [F1' * s1 * F1 + F2' * s2 * F2, sparse(npixels, npixels); sparse(npixels, npixels), F1' * s1 * F1 + F2' * s2 * F2];
    bs = 0;
    
    
    % Add different matrix and vector contributions
    A = Ad + As;
    b = bd + bs;
   
end

