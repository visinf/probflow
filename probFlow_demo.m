%% Demo file showing the application of ProbClassicA and ProbFlowFields
% Images are taken from Middlebury [2] and Sintel [11] optical flow
% benchmarks.
% Initial flow estimate for ProbFlowFields is computed with sparse
% FlowFields matches [1] and EpicFlow interpolation [36].


% Please modify the below path to add the directory with code of Sun et al. [43]!
addpath(genpath('/path/to/flow_code/utils'));

%% Examplary use of ProbFlowFields

addpath(genpath('./ProbFlowFields'))

% Load necessary data (Sintel training, temple_3, frames 49/50) 
frame1 = imread('./data/probFlowFields_frame1.png');
frame2 = imread('./data/probFlowFields_frame2.png');
gt = readFlowFile('./data/probFlowFields_gt.flo');
muInit = readFlowFile('./data/probFlowFields_muInit.flo'); % Loading of pre-computed FlowFields matches

% Set parameters
params = struct;
params.sorSolver = true; % requires previous compilation of sor.c

% Compute estimates and variances
disp('ProbFlowFields in process...')
[mu,sig,~,~,~] = probFlowFields(frame1,frame2,muInit,params);

% Show initial estimate, ProbFlowFields estimate and ground truth flow
img_muInit = flowToColor(muInit); 
figure('name','ProbFlowFields initial estimate'); imshow(img_muInit,'Border','tight'); pos = get(gca, 'position');
img_mu = flowToColor(mu); 
figure('name','ProbFlowFields estimate'); imshow(img_mu,'Border','tight');
img_gt = flowToColor(gt); 
figure('name','ProbFlowFields ground truth'); imshow(img_gt,'Border','tight');

% Depict uncertainty map
uncertainty = sum(log(sig),3); 
figure('name','ProbFlowFields uncertainty'); imshow(img_muInit,'Border','tight'); imagesc(uncertainty); colormap(jet); axis off; set(gca,'position',pos);
disp('ProbFlowFields done, processing paused')
pause

%% Examplary use of ProbClassicA
close all
addpath(genpath('./ProbClassicA'))

% Load necessary data (Middlebury training, RubberWhale) 
frame1 = imread('./data/probClassicA_frame1.png');
frame2 = imread('./data/probClassicA_frame2.png');
gt = readFlowFile('./data/probClassicA_gt.flo');
muInit = zeros(size(gt)); % Initialization with zero

% Compute estimates and variances
disp('ProbClassicA in process... (this may take some time)')
[mu,sig,~,~] = probClassicA(frame1,frame2,muInit);

% Show initial estimate, ProbFlowFields estimate and ground truth flow
img_mu = flowToColor(mu); 
figure('name','ProbClassicA estimate'); imshow(img_mu,'Border','tight'); pos = get(gca, 'position');
img_gt = flowToColor(gt); 
figure('name','ProbClassicA ground truth'); imshow(img_gt,'Border','tight');

% Depict uncertainty map
uncertainty = sum(log(sig),3); 
figure('name','ProbClassicA uncertainty'); imshow(img_gt,'Border','tight'); imagesc(uncertainty); colormap(jet); axis off; set(gca,'position',pos);