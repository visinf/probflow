function [muHat, sigHat, kHat] = probClassicA_minStepHat(mu, muHatInit, sigHatInit, lamC, lamN, gsmN, medianSize, nIter)
% Perform one minimization step for (muHat,sigHat,kHat)

    % Define help constants
    [n,m,~] = size(mu);
    hMedianSize = (medianSize - 1)/2;    
    centerPixel = (medianSize^2  - 1)/2 + 1;

    % Initialize variables
    muHat = muHatInit;
    sigHat = sigHatInit;
    
    % Precompute neighborhood values for muHat
    muHat_old = padarray(muHat, hMedianSize*[1 1], 'symmetric', 'both'); % symmetric padding to obtain a neighborhood for every pixel
    neighbors = im2col(muHat_old, medianSize*[1 1],'sliding'); % obtain neighborhoods of size [medianSize medianSize]
    neighbors = [neighbors(1:centerPixel-1,:); neighbors(centerPixel+1:end,:)]; % remove center pixel as it is not included in the neighborhood

    % Precompute neighborhood values for sigHat
    sigHat_old = padarray(sigHat, hMedianSize*[1 1], 'symmetric', 'both');
    neighborsSig = im2col(sigHat_old, medianSize*[1 1],'sliding');
    neighborsSig = [neighborsSig(1:centerPixel-1,:); neighborsSig(centerPixel+1:end,:)];

    % Initial update for kHat and sigHat
    kHat = updateKhat(muHat, sigHat, neighbors, neighborsSig, gsmN, lamN);
    sigHat = updateSigHat(kHat, gsmN, lamC, lamN, n, m);

    % Alternating estimation of muHat, sigHat and kHat
    for j = 1:nIter
        
        % Precompute neighborhood values for muHat
        muHat_old = padarray(muHat, hMedianSize*[1 1], 'symmetric', 'both');
        neighbors = im2col(muHat_old, medianSize*[1 1],'sliding');
        neighbors = [neighbors(1:centerPixel-1,:); neighbors(centerPixel+1:end,:)];
        
        % Precompute neighborhood values for sigHat
        sigHat_old = padarray(sigHat, hMedianSize*[1 1], 'symmetric', 'both');
        neighborsSig = im2col(sigHat_old, medianSize*[1 1],'sliding');
        neighborsSig = [neighborsSig(1:centerPixel-1,:); neighborsSig(centerPixel+1:end,:)];

        muHat = updateMuHat(mu, kHat, neighbors, gsmN, lamC, lamN);
        kHat = updateKhat(muHat, sigHat, neighbors, neighborsSig, gsmN, lamN);
        sigHat = updateSigHat(kHat, gsmN, lamC, lamN, n, m);   
    end

end

