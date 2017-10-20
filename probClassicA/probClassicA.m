function [mu,sig,kd,ks,muHat,sigHat,kHat] = probClassicA(frame1,frame2,muInit,params)
% Jointly estimate optical flow and variances between images frame1 and
% frame2 using ProbClassicA

%% Parameters

    if ~exist('params', 'var')
        params = struct;
    end
    
    % Params data linearization
    nWarp = getParam(params, 'nWarp', 10);
    
    % Params energy function
    lamD = getParam(params, 'lamD', 0.0959);
    lamS = getParam(params, 'lamS', 0.1000);
    lamC = getParam(params, 'lamC', logspace(log10(1e-4), log10(1e2), nWarp));
    lamN = getParam(params, 'lamN', 0.0053);

    % Params optimization (mu,sig,k)
    nIter = getParam(params, 'nIter', 5);
    initSig = getParam(params, 'initSig', 1e-7);
    
    % Params optimization (muHat,sigHat,kHat)
    nbhSize = getParam(params, 'nbhSize', 5);
    nIterHat = getParam(params, 'nIterHat', 3);

    % Params image pre-processing
    strucTextDec = getParam(params, 'strucTextDec', true);
    
    % Params GSM penalties
    pathGSM = getParam(params, 'pathGSM', './gsm_params');
    nameGSMd = getParam(params, 'nameGSMd', 'L5_probClassicA_data');
    nameGSMs = getParam(params, 'nameGSMs', 'L5_probClassicA_smooth');
    nameGSMn = getParam(params, 'nameGSMn', 'L5_probClassicA_nonlocal');

    % Params image derivatives
    h = getParam(params, 'h', [1 -8 0 8 -1]/12);
    tempAv = getParam(params, 'tempAv', 0.5);
    interMeth = getParam(params, 'interMeth', 'bi-cubic');

    % Params image pyramids
    minPixelPyramid = getParam(params, 'minPixelPyramid', 6);
    maxPixelPyramid = getParam(params, 'maxPixelPyramid', 16);
    pyramidScale1 = getParam(params, 'pyramidScale1', 0.5);
    pyramidLevels1 = getParam(params, 'pyramidLevels1', []);
    pyramidScale2 = getParam(params, 'pyramidScale2', 0.8);
    pyramidLevels2 = getParam(params, 'pyramidLevels2', 2);

%% Initialization

    % Initialize gsm penalties
    gsmD3 = load(sprintf('%s/%s',pathGSM,nameGSMd));
    gsmS3 = load(sprintf('%s/%s',pathGSM,nameGSMs));
    gsmN = load(sprintf('%s/%s',pathGSM,nameGSMn));

    % Convert images to grayscale and apply structure texture
    % decomposition if required
    im1 = double(rgb2gray(frame1));
    im2 = double(rgb2gray(frame2));
    images = cat(3,im1,im2);
    
    if strucTextDec 
        images = structure_texture_decomposition_rof(images,1/8, 100, 0.95);
    end

    % Build image pyramid
    [pyramid_images, pyramidLevels1] = buildImagePyramid(images, pyramidLevels1, pyramidScale1, minPixelPyramid, maxPixelPyramid);
    currentImages = pyramid_images{pyramidLevels1};
    im1 = currentImages(:,:,1);
    [n,m] = size(im1);

    % Initialize estimates
    mu = zeros(n,m,2);
    factorU = m / size(muInit,2);
    mu(:,:,1) = imresize(muInit(:,:,1), [n,m])*factorU;
    factorV = n / size(muInit,1);
    mu(:,:,2) = imresize(muInit(:,:,2), [n,m])*factorV;
    
    muHat = zeros(n,m,2);
    if ~isempty(nbhSize)
        muHat(:,:,1) = denoise_LO(mu(:,:,1), nbhSize, lamC(1)/lamN, nIterHat);
        muHat(:,:,2) = denoise_LO(mu(:,:,2), nbhSize, lamC(1)/lamN, nIterHat);
    end

    sig1 = initSig * ones(n,m);
    sig2 = initSig * ones(n,m);
    sigHat1 = initSig * ones(n,m);
    sigHat2 = initSig * ones(n,m);

%% Flow calculation

    % Iterate over graduated non-convexity levels
    for gnc = 1:3
        if gnc == 1
            % Use quadratic penalties for gnc = 1
            gsmD = struct;
            gsmD.sigma = 1;
            gsmD.pi = 1;
            
            gsmS = struct;
            gsmS.sigma = 1;
            gsmS.pi = 1;
            
            pyramidLevels = pyramidLevels1;

        elseif gnc == 2
            % Use mixture of quadratic and non-convex penalties for gnc = 2
            gsmD = struct;
            gsmD.sigma = [gsmD3.sigma;1];
            gsmD.pi = [gsmD3.pi*0.5;0.5];
            
            gsmS = struct;
            gsmS.sigma = [gsmS3.sigma;1];
            gsmS.pi = [gsmS3.pi*0.5;0.5];
            
            % Build image pyramid with different image scaling and pyramid
            % levels for remaining gnc steps
            pyramid_images = buildImagePyramid(images, pyramidLevels2, pyramidScale2, minPixelPyramid, maxPixelPyramid);
            pyramidLevels = pyramidLevels2;

        elseif gnc == 3
            % Use non-convex penalties for gnc = 3
            gsmD = gsmD3;
            gsmS = gsmS3;

        end
        
        % Iterate over pyramid levels
        for plevel = pyramidLevels:-1:1

            currentImages = pyramid_images{plevel};
            im1 = currentImages(:,:,1);
            [n,m] = size(im1);

            % Upsample variables to current size [n, m]
            [mu,sig1,sig2] = upsampleVariables(mu,sig1,sig2,[n m]);
            [muHat,sigHat1,sigHat2] = upsampleVariables(muHat,sigHat1,sigHat2,[n m]);

            % Initialization filter matrices
            F1 = make_convn_mat([-1;1], [n m], 'valid', 'same');
            F2 = make_convn_mat([-1,1], [n m], 'valid', 'same');

            % Perform nWarp warping steps per pyramid level
            for l = 1:nWarp
                currentLamC = lamC(l);   

                % Calculate image derivatives
                mu0 = mu;
                [It, Ix, Iy] = partial_deriv(currentImages, mu0, interMeth, h, tempAv);

                % Update step (mu,sig,k)
                [mu, sig1, sig2, kd, ks] = probClassicA_minStep(mu0, muHat, sig1, sig2, It, Ix, Iy, F1, F2, lamD, lamS, currentLamC, gsmD, gsmS, nIter);

                % Update step (muHat,sigHat,kHat)
                if ~isempty(nbhSize)
                    muHatInit = mu(:,:,1);
                    sigHatInit = sigHat1;
                    [muHat1,sigHat1,kHat1] = probClassicA_minStepHat(mu(:,:,1), muHatInit, sigHatInit, currentLamC, lamN, gsmN, nbhSize, nIterHat);

                    muHatInit = mu(:,:,2);
                    sigHatInit = sigHat2;
                    [muHat2,sigHat2,kHat2] = probClassicA_minStepHat(mu(:,:,2), muHatInit, sigHatInit, currentLamC, lamN, gsmN, nbhSize, nIterHat);
                    
                    muHat = cat(3,muHat1,muHat2);
                end
            end

        end

    end

%% Concatenate results

    sig = cat(3,sig1,sig2);
    sigHat = cat(3,sigHat1,sigHat2);
    kHat = cat(4,kHat1,kHat2);

end

