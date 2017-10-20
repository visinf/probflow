function [Ix2,IxIy,Iy2,IxIt,IyIt,It2,Ixx2,IxxIxy,Ixy2,IxxIxt,IxyIxt,Ixt2,Iyx2,IyxIyy,Iyy2,IyxIyt,IyyIyt,Iyt2] ...
    = normalizeDerivatives(It, Ix, Iy, Ixx, Ixy, Iyx, Iyy, Ixt, Iyt, energyOpt)
% Normalization of image derivatives following Zimmer et al. [49]

    % Initialize parameters
    epsD1 = energyOpt.epsD1;
    epsD2 = energyOpt.epsD2;

    % Normalize first-order image derivatives 
    Ix2 = Ix.^2;
    IxIy = Ix .* Iy;
    Iy2 = Iy.^2;
    IxIt = It .* Ix;
    IyIt = It .* Iy;
    It2 = It.^2;
    
    normalization = Ix2 + Iy2 + epsD1;
    Ix2 = sum(Ix2./normalization,3);
    IxIy = sum(IxIy./normalization,3);
    Iy2 = sum(Iy2./normalization,3);
    IxIt = sum(IxIt./normalization,3);
    IyIt = sum(IyIt./normalization,3);
    It2 = sum(It2./normalization,3);
    
    
    % Normalize second-order image derivatives 
    Ixx2 = Ixx.^2;
    IxxIxy = Ixx .* Ixy;
    Ixy2 = Ixy.^2;
    IxxIxt = Ixx .* Ixt;
    IxyIxt = Ixy .* Ixt;
    Ixt2 = Ixt .* Ixt;
    
    normalization = Ixx2 + Ixy2 + epsD2;
    Ixx2 = sum(Ixx2./normalization,3);
    IxxIxy = sum(IxxIxy./normalization,3);
    Ixy2 = sum(Ixy2./normalization,3);
    IxxIxt = sum(IxxIxt./normalization,3);
    IxyIxt = sum(IxyIxt./normalization,3);
    Ixt2 = sum(Ixt2./normalization,3);
    
    Iyx2 = Iyx.^2;
    IyxIyy = Iyx .* Iyy;
    Iyy2 = Iyy.^2;
    IyxIyt = Iyx .* Iyt;
    IyyIyt = Iyy .* Iyt;
    Iyt2 = Iyt .* Iyt;
    
    normalization = Iyx2 + Iyy2 + epsD2;
    Iyx2 = sum(Iyx2./normalization,3);
    IyxIyy = sum(IyxIyy./normalization,3);
    Iyy2 = sum(Iyy2./normalization,3);
    IyxIyt = sum(IyxIyt./normalization,3);
    IyyIyt = sum(IyyIyt./normalization,3);
    Iyt2 = sum(Iyt2./normalization,3);

end

