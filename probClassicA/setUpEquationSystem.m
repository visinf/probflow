function [A, b] = setUpEquationSystem(mu0, kd, ks, Ix2, Ixy, Iy2, Itx, Ity, F1, F2, gsmD, gsmS, lamD, lamS, Ac, bc)
% Determine equation system needed to update flow estimate mu, cf.
% equations (31) and (32) of the supplemental material
    
    % Initialize help variables
    [n,m,~] = size(mu0);
    npixels = n * m;

    % Determine contribution of data term
    termD = bsxfun(@rdivide, kd, 2*gsmD.sigma);
    termD = sum(termD, 1);
    termD = 2 * lamD * termD(:);
    
    % Compute data matrix Ad and vector bd
    tmp = termD .* Ix2(:);
    dxx = spdiags(tmp, 0, npixels, npixels);
    tmp = termD .* Iy2(:);
    dyy = spdiags(tmp, 0, npixels, npixels);
    tmp = termD .* Ixy(:);
    dxy = spdiags(tmp, 0, npixels, npixels);
    
    Ad = [dxx dxy; dxy dyy];
    bd = Ad * mu0(:) - [termD .* Itx(:); termD .* Ity(:)];

    % Determine contribution of smoothness terms
    termS = bsxfun(@rdivide,ks,2*gsmS.sigma);
    termS = sum(termS, 1);
    s1 = termS(:,:,1);
    s1 = spdiags(2 * s1(:), 0, npixels, npixels);
    s2 = termS(:,:,2);
    s2 = spdiags(2 * s2(:), 0, npixels, npixels);
    s3 = termS(:,:,3);
    s3 = spdiags(2 * s3(:), 0, npixels, npixels);
    s4 = termS(:,:,4);
    s4 = spdiags(2 * s4(:), 0, npixels, npixels);
    
    % Compute smoothness matrix As
    As = lamS * [F1' * s1 * F1 + F2' * s2 * F2, sparse(npixels, npixels); ...
        sparse(npixels, npixels), F1' * s3 * F1 + F2' * s4 * F2];
    
    % Add different matrix and vector contributions
    A = Ad + As + Ac;
    b = sparse(bd + bc);
   
end

