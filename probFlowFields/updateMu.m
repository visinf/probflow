function [mu] = updateMu(mu0, kd1, kd2, ks, Ix2, IxIy, Iy2, IxIt, IyIt, Ixx2, IxxIxy, Ixy2, IxxIxt, IxyIxt, Iyx2, IyxIyy, Iyy2, IyxIyt, IyyIyt, gsmD1, gsmD2, gsmS, sorSolver, energyOpt)
% Perform update step of flow estimates mu, i.e. solve equation system 
% A*x = b as given in equation (31) of the supplemental material

    % Initialize help variable
    dmu0 = zeros(size(mu0));

    % Use SOR solver if required else utilize Matlab backslash
    % operation to solve equation system
    if sorSolver
        [a11,a12,a22,b1,b2,smoothHor,smoothVert] = setUpSOR(mu0, kd1, kd2, ks, Ix2, IxIy, Iy2, IxIt, IyIt, Ixx2, IxxIxy, Ixy2, IxxIxt, IxyIxt, Iyx2, IyxIyy, Iyy2, IyxIyt, IyyIyt, gsmD1, gsmD2, gsmS, energyOpt);
        [du, dv] = sor_probflow(single(dmu0(:,:,1)'), single(dmu0(:,:,2)'), single(a11'), single(a12'), single(a22'), single(b1'), single(b2'), single(smoothHor'), single(smoothVert'), single(30), single(1.9));
        dmu = cat(3,du',dv');
        mu = mu0 + dmu;
    else
        [A,b] = setUpEquationSystem(mu0, kd1, kd2, ks, Ix2, IxIy, Iy2, IxIt, IyIt, Ixx2, IxxIxy, Ixy2, IxxIxt, IxyIxt, Iyx2, IyxIyy, Iyy2, IyxIyt, IyyIyt, gsmD1, gsmD2, gsmS, energyOpt);
        mu = full(A\b);
        mu = reshape(mu, size(mu0));
    end
    
    
      
end

