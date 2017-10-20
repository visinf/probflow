function [a11,a12,a22,b1,b2,smoothHor,smoothVert] = setUpSOR(mu0, kd1, kd2, ks, Ix2, IxIy, Iy2, IxIt, IyIt, Ixx2, IxxIxy, Ixy2, IxxIxt, IxyIxt, Iyx2, IyxIyy, Iyy2, IyxIyt, IyyIyt, gsmD1, gsmD2, gsmS, energyOpt)
% Determine SOR variables needed to update flow estimate mu

    % Initialize help variables
    [n, m] = size(Ix2);
    lamD1 = energyOpt.lamD1;
    lamD2 = energyOpt.lamD2;
    lamS = energyOpt.lamS;
    filter_smooth = energyOpt.filter_smooth;
    smoothVert = zeros(n, m);
    smoothHor = zeros(n, m);
    
    % Determine contribution of data term 1
    termD = bsxfun(@rdivide, kd1, 2*gsmD1.sigma);
    termD = sum(termD,1);
    termD = reshape(2 * lamD1 * termD, [n,m]);

    a11 = termD.*Ix2;
    a22 = termD.*Iy2;
    a12 = termD.*IxIy;
    
    b1 = - termD .* IxIt;
    b2 = - termD .* IyIt;

    
    % Determine contribution of data term 2
    termD = bsxfun(@rdivide, kd2, 2*gsmD2.sigma);
    termD = sum(termD, 1);
    termD = reshape(2 * lamD2 * termD, [n,m]);
    
    a11 = a11 + termD .*(Ixx2+Iyx2);
    a22 = a22 + termD .*(Ixy2+Iyy2);
    a12 = a12 + termD .*(IxxIxy+IyxIyy);

    b1 = b1 - termD .* (IxxIxt + IyxIyt);  
    b2 = b2 - termD .* (IxyIxt + IyyIyt);
    

    % Determine contribution of smoothness term 1
    termS = bsxfun(@rdivide,ks,2*gsmS.sigma);
    termS = sum(termS, 1);
    
    s1 = reshape(2 * termS(:,:,1), [n,m]);
    filt = filter_smooth{1};
    lamS_tmp = 0.5 * imfilter(lamS, abs(filt), 'conv', 'replicate');
    smoothVert(1:end-1,:) = s1(1:end-1,:) .* lamS_tmp(1:end-1,:);
    
    tmp = smoothVert .* imfilter(mu0(:,:,1), filt, 'conv', 'replicate');
    b1(1:end-1,:) = b1(1:end-1,:) + tmp(1:end-1,:);
    b1(2:end,:) = b1(2:end,:) - tmp(1:end-1,:);
    tmp = smoothVert .* imfilter(mu0(:,:,2),filt,'conv','replicate');
    b2(1:end-1,:) = b2(1:end-1,:) + tmp(1:end-1,:);
    b2(2:end,:) = b2(2:end,:) - tmp(1:end-1,:); 
    
    
    % Determine contribution of smoothness term 2
    s2 = reshape(2 * termS(:,:,2), [n,m]);
    filt = filter_smooth{2};
    lamS_tmp = 0.5 * imfilter(lamS,abs(filt),'conv','replicate'); 
    smoothHor(:,1:end-1) = s2(:,1:end-1) .* lamS_tmp(:,1:end-1);
    
    tmp = smoothHor .* imfilter(mu0(:,:,1),filt,'conv','replicate');
    b1(:,1:end-1) = b1(:,1:end-1) + tmp(:,1:end-1);
    b1(:,2:end) = b1(:,2:end) - tmp(:,1:end-1);
    tmp = smoothHor .* imfilter(mu0(:,:,2),filt,'conv','replicate');
    b2(:,1:end-1) = b2(:,1:end-1) + tmp(:,1:end-1);
    b2(:,2:end) = b2(:,2:end) - tmp(:,1:end-1); 
   
end

