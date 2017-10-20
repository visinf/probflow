function [mu,sig,kd1,kd2,ks] = probFlowFields(frame1,frame2,muInit,params)
% Jointly estimate optical flow and variances between images frame1 and
% frame2 using ProbFlowFields
    
%% Parameters
    if ~exist('params', 'var')
        params = struct;
    end
    
    % Params data linearization
    nWarp = getParam(params, 'nWarp', 5);
    
    % Params energy function
    energyOpt = struct;
    energyOpt.lamD1 = getParam(params, 'lamD1', 0);
    energyOpt.lamD2 = getParam(params, 'lamD2', 0.01);
    energyOpt.lamS = getParam(params, 'lamS', 0.42);
    energyOpt.epsD1 = getParam(params, 'epsD1', 0.01);
    energyOpt.epsD2 = getParam(params, 'epsD2', 0.01);
    energyOpt.filter_smooth = getParam(params, 'filter_smooth', {[1;-1],[1,-1]});
    kappa = getParam(params, 'kappa', 0.01);
    
    % Params optimization
    nIter = getParam(params, 'nIter', 5);
    initSig = getParam(params, 'initSig', 1e-7);
    sorSolver = getParam(params, 'sorSolver', true);
    
    % Params image pre-processing
    sigmaFilt = getParam(params, 'sigmaFilt', 1.1);
    
    % Params GSM penalties
    pathGSM = getParam(params, 'pathGSM', './gsm_params');
    nameGSMd1 = getParam(params, 'nameGSMd', 'L10_probFlowFields_data');
    nameGSMd2 = getParam(params, 'nameGSMn', 'L10_probFlowFields_data');
    nameGSMs = getParam(params, 'nameGSMs', 'L10_probFlowFields_smooth');
    
    % Params image derivatives
    h = getParam(params, 'h', [1 -8 0 8 -1]/12);
    

%% Initialization

    % Initialize gsm penalties
    gsmD1 = load(sprintf('%s/%s',pathGSM,nameGSMd1));
    gsmD2 = load(sprintf('%s/%s',pathGSM,nameGSMd2));
    gsmS = load(sprintf('%s/%s',pathGSM,nameGSMs));

    % Pre-process images 
    im1 = imgaussfilt(double(frame1),sigmaFilt,'FilterSize',2*(floor(3*sigmaFilt)+1)+1); 
    im2 = imgaussfilt(double(frame2),sigmaFilt,'FilterSize',2*(floor(3*sigmaFilt)+1)+1);
    
    % Pre-compute spatially varying trade-off parameter lamS
    energyOpt.lamS = energyOpt.lamS * calculateLamS(im1,kappa,h);
    
    % Initialize estimates
    mu = muInit;
    [n,m,~] = size(im1);
    sig1 = initSig * ones(n,m);
    sig2 = initSig * ones(n,m);


%% Flow calculation

    % Perform nWarp warping steps
    for l = 1:nWarp
        
        % Compute image derivatives
        mu0 = mu;
        [It,Ix,Iy,Ixx,Ixy,Iyx,Iyy,Ixt,Iyt] = calculateImageDerivatives(im1,im2,mu0,h);
        
        % Update step (mu,sig,k)
        [mu, sig1, sig2, kd1, kd2, ks] = probFlowFields_minStep(mu0, sig1, sig2, It, Ix, Iy, Ixx, Ixy, Iyx, Iyy, Ixt, Iyt, gsmD1, gsmD2, gsmS, nIter, sorSolver, energyOpt);
        
    end


%% Concatenate results
    
    sig = cat(3,sig1,sig2);


end

