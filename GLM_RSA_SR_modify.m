clear; close all; clc;
set(groot, 'defaultfigurewindowstyle', 'docked')
set(groot, 'defaultaxesfontsize', 20)
set(groot, 'defaultlinelinewidth', 2)

boolean = @logical;

TR = 2;
frames_run = 168;
frames_trial = 14;
trial_length = TR*frames_trial;

nCate = 24;

label_file = './label/Label_prep.nii.gz';
D = load('bold_sub-0061.mat');
bold = D.bold;
behav = array2table(D.behav, 'variablenames', D.behav_label);
confound = array2table(D.confounds, 'variablenames', D.confounds_label);

nVoxel = size(bold,2);

Plotting = true;
%% Bold z transform
boldz = nan(size(bold)); % (voxel x run) wise z transform
nRun = size(bold, 1)/frames_run;

for iRun = 1:nRun
    idx = ((iRun-1)*frames_run+1):(iRun)*frames_run;
    boldz(idx,:) = zscore(bold(idx,:), [], 1);
end


%% pre-setting for GLM

Spc = 0.5; % should be TR/integer
n_down = TR/Spc;

hrf = spm_hrf(Spc);

tmp = 0:Spc:trial_length;
x = tmp(1:end-1)';

tmp = zeros(size(x));
box_stim = tmp; box_stim(x<1.5) = 1;
box_choice_e = tmp; box_choice_e((6<=x)&(x<7.5)) = 1;
box_choice_l = tmp; box_choice_l((12<=x)&(x<13.5)) = 1;
box_est = tmp; box_est((18<=x)&(x<=22)) = 1;

conv_hrf = @(y,hrf) conv([zeros(length(hrf)-1,1); y(:)], hrf(:), 'valid');

if Plotting
    xb = (-1:Spc:-Spc)'; zb = zeros(size(xb));
    xa = (trial_length:Spc:trial_length+10)'; za = zeros(size(xa));
    x2 = [xb; x; xa];
    
    figure()
    plot(x2, [zb; box_stim; za], 'color', 'red'); xlim([-1, 40]); ylim([-0.1, 1.1]); hold on;
    plot(x2, [zb; box_choice_e; za], 'color', [0.8,0.8,0]); ylim([-0.1, 1.1]);
    plot(x2, [zb; box_est; za], 'color', 'blue'); ylim([-0.1, 1.1]);
    
    plot(x2, conv_hrf([zb; box_stim; za], hrf), 'color', 'red');  ylim([-0.1, 1.1]);
    plot(x2, conv_hrf([zb; box_choice_e; za], hrf), 'color', [0.8,0.8,0]);  ylim([-0.1, 1.1]);
    plot(x2, conv_hrf([zb; box_est; za], hrf), 'color', 'blue'); ylim([-0.1, 1.1]);
    
    legend('Stimulus', 'Choice', 'Estimation')
    
    
end

%% Nuisance regressor : discrete cosine transform basis (for highpass filter)

cutoff = 0.005;

nT = (frames_run*TR)/(1/cutoff);

basis_rep = 0:0.5:nT;
DCT_basis = nan(frames_run, length(basis_rep));

for i = 1:length(basis_rep)
    y = cos(linspace(0, 2*pi*basis_rep(i), frames_run+1));
    DCT_basis(:,i) = y(1:end-1);
end


if Plotting
    figure();
    plot(DCT_basis);
    title('Discrete cosine transform basis')
end



%% Nuisance regressor : WM+CSF+Motion in fMRIprep result

confounds = table2array(confound);
nNui = size(DCT_basis,2) + size(confounds,2);
nuisance_regressor = zeros(size(bold,1),nNui*nRun);

for iRun = 1:nRun
    idx1 = ((iRun-1)*frames_run+1):(iRun)*frames_run;
    idx2 = ((iRun-1)*nNui+1):(iRun)*nNui;
    
    conf_run = zscore(confounds(idx1,:), [], 1);
    % plot(conf_run); legend(D.confounds_label)
    x = [DCT_basis, conf_run];
    nuisance_regressor(idx1,idx2) = x;
    
end

if Plotting
    subplot(1,2,1)
    imagesc(x); colormap(jet); set(gca, 'clim', [-4,4])
    subplot(1,2,2)
    imagesc(nuisance_regressor)
    
end
%% make design matrix

% stimulus regressor matrix

stim = behav.stim;
cate = discretize(stim, nCate);
stim_c = zeros(length(stim), nCate);

for iC = 1:nCate
    stim_c(:,iC) = cate == iC;
end

% imagesc(stim_c')
stim_design = kron(stim_c, box_stim); % stimulus box car

timing = D.behav(:,strcmp(D.behav_label,'Timing'));
choice_design = kron(timing == 1, box_choice_e) +...
    kron(timing == 2, box_choice_l); %choice box car

estim_design = kron(ones(size(D.behav,1),1), box_est); %estimation box car

design_box = [stim_design, choice_design, estim_design];
% design_box = [stim_design];

design_conv = zeros(size(design_box));


for iCol = 1:size(design_box,2)
    design_conv(:,iCol) = conv_hrf(design_box(:,iCol), hrf);
end

design_fin = downsample(design_conv, n_down);

if Plotting
    figure()
    subplot(1,2,1);
    imagesc(design_box(1:1000,:))
    subplot(1,2,2);
    imagesc(design_conv(1:1000,:))
    
    figure()
    subplot(2,1,1)
    plot(design_box(:,1)); xlim([1900, 2100])
    subplot(2,1,2)
    plot(design_conv(:,1)); xlim([1900, 2100])
    
    
    X2 = [10*design_fin, nuisance_regressor];
    figure(); imagesc(X2); colormap(jet); set(gca, 'clim',[-5,5])
    
end

%% estimating coefficient (activity pattern)

Y = boldz;
X = [design_fin, nuisance_regressor];

coef = X\Y;
coef_stim = coef(1:size(stim_design,2),:);
 
% mdl = fitlm(X, Y(:,1));
% disp(mdl.Coefficients(2,:));

% estimating T value

pred = X*coef;
Err = Y-pred;
dof = size(X,1)-size(X,2);
MSE = sum(Err.^2,1)/dof;
V = diag(inv(X'*X))*MSE; % variance of estimated coefficients of multiple linear regression
se_coef = sqrt(V);
se_coef_stim = se_coef(1:size(stim_design,2),:);
T_stim = coef_stim./se_coef_stim;

%% Extract label index

label_nifti = niftiread(label_file);
coords = D.coords;

label = nan(1, nVoxel);

for i = 1:length(label)
    label(i) = label_nifti(coords(i,1), coords(i,2),coords(i,3));
end

%% Mean activity of ROI (trial)

close all;

ROI_name = {'V1','V2','MT','Motor'};
ROI_idx = {ismember(label, [1,2]), ismember(label, [3,4]),...
    ismember(label, [5,6]), ismember(label, [7,8])};

for i = 1:length(ROI_name)
    vox_idx = ROI_idx{i};
    y = boldz(:,vox_idx);
    y_trial =reshape(y, frames_trial,[],size(y,2));
    
    y_mean1 = mean(y_trial(:,timing == 1,:), [2,3]);
    y_mean2 = mean(y_trial(:,timing == 2,:), [2,3]);
    
    if Plotting
        figure()
        plot(y_mean1); hold on;
        plot(y_mean2);
        legend('Early', 'Late', 'Location','Northeastoutside')
        title(ROI_name{i});
    end
    % saveas(gcf,['HW1_' ROI_name{i} '.png']);
end

%% Conceptual RDM

Conc1 = nan(nCate,nCate);
for i = 1:nCate
    for j = 1:nCate
        Conc1(i,j) = min(abs(i-j), nCate-abs(i-j));
        %          Hypo1(i,j) = 1-cos((i-j)/nCate*pi*2);
    end
end

Conc2 = ones(nCate,nCate);
Conc2(1:round(nCate/2),1:round(nCate/2))= 0;
Conc2(round(nCate/2)+1:end,round(nCate/2)+1:end)= 0;

if Plotting
    subplot(1,2,1)
    imagesc(Conc1); axis square
    title('Abs ori difference')
    subplot(1,2,2)
    imagesc(Conc2); axis square
    title('Categorical')    
end
% saveas(gcf,'HW3_Conceptual RDM.png');

%% Differnt types of RDM
close all;
for idx = 1:4
vox_idx = ROI_idx{idx};

activ_ROI = coef_stim(:,vox_idx);
Conc = Conc1;

Methods = {'Euclidean', 'Correlation'};

for i = 1:length(Methods)
    
    RDM = squareform(pdist(activ_ROI, Methods{i})); % pdist is pairwise distance function in matlab
    RDM(boolean(eye(size(RDM)))) = nan;
    figure()
    subplot(1,2,1)
    imagesc(RDM);  colorbar; axis square;
    title(char([ROI_name{idx} '-' Methods{i}]))
    
    subplot(1,2,2)
    RC = corr([Conc(:), RDM(:)], 'type', 'kendall', 'rows', 'complete');
    plot(Conc(:), RDM(:), 'o'); axis square
    title(sprintf('Rank correlation : %.04f', RC(2,1)))
    
    % saveas(gcf,['HW2_' ROI_name{idx} '_' Methods{i} '.png']);
end
end

%% Mean pattern subtraction
for idx = 1:4
vox_idx = ROI_idx{idx};

activ_ROI = coef_stim(:,vox_idx);
activ_ROI_demean = activ_ROI-mean(activ_ROI,1);

RDM = 1-corr(activ_ROI');
RDM(boolean(eye(size(RDM)))) = nan;
figure()
subplot(1,2,1);
imagesc(RDM);  colorbar;axis square
title([ROI_name{idx} '-Correlation (Raw)'])
subplot(1,2,2);
plot(Conc(:), RDM(:), 'o'); axis square
RC = corr([Conc(:), RDM(:)], 'type', 'kendall', 'rows', 'complete');
title(sprintf('Rank correlation : %.04f', RC(2,1)))


RDM = 1-corr(activ_ROI_demean');
RDM(boolean(eye(size(RDM)))) = nan;
figure()
subplot(1,2,1);
imagesc(RDM);  colorbar;axis square
title(sprintf([ROI_name{idx} '-Correlation \n (Mean pattern subtraction)']))
subplot(1,2,2);
plot(Conc(:), RDM(:), 'o'); axis square
RC = corr([Conc(:), RDM(:)], 'type', 'kendall', 'rows', 'complete');
title(sprintf('Rank correlation : %.04f', RC(2,1)))

% saveas(gcf,['HW2_' ROI_name{idx} '_subtracted.png']);
end


%% T-statistics based
activ_ROI1 = coef_stim(:,vox_idx);
activ_ROI2 = T_stim(:,vox_idx);

activ_ROI_demean1 = activ_ROI1-mean(activ_ROI1,1);
activ_ROI_demean2 = activ_ROI2-mean(activ_ROI2,1);

RDM = 1-corr(activ_ROI_demean1');
RDM(boolean(eye(size(RDM)))) = nan;
figure()
subplot(1,2,1);
imagesc(RDM);  colorbar;axis square
title(sprintf('Correlation (Betas)'))
subplot(1,2,2);
plot(Conc(:), RDM(:), 'o'); axis square
RC = corr([Conc(:), RDM(:)], 'type', 'kendall', 'rows', 'complete');
title(sprintf('Rank correlation : %.04f', RC(2,1)))


RDM = 1-corr(activ_ROI_demean2');
RDM(boolean(eye(size(RDM)))) = nan;
figure()
subplot(1,2,1);
imagesc(RDM);  colorbar;axis square
title(sprintf('Correlation (T)'))
subplot(1,2,2);
plot(Conc(:), RDM(:), 'o'); axis square
RC = corr([Conc(:), RDM(:)], 'type', 'kendall', 'rows', 'complete');
title(sprintf('Rank correlation : %.04f', RC(2,1)))


%% RDM of each ROI
close all;

for cnum = 1:2
    if cnum==1
    Conc = Conc1;
    elseif cnum==2
    Conc = Conc2;
    end

for i = 1:length(ROI_idx)
    vox_idx = ROI_idx{i};

    activ_ROI2 = T_stim(:,vox_idx);
    activ_ROI_demean2 = activ_ROI2-mean(activ_ROI2,1);
    
    RDM = 1-corr(activ_ROI_demean2');
    RDM(boolean(eye(size(RDM)))) = nan;
    figure()
    subplot(1,2,1);
    imagesc(RDM);  colorbar;axis square
    title(strcat(ROI_name{i}, '-Concept: H', string(cnum)))
    subplot(1,2,2);
    plot(Conc(:), RDM(:), 'o'); axis square
    RC = corr([Conc(:), RDM(:)], 'type', 'kendall', 'rows', 'complete');
    title(sprintf('Rank correlation : %.04f', RC(2,1)))

% saveas(gcf,strcat('HW3_', ROI_name{i}, '-Concept-H', string(cnum), '.png'));
end
end



%% plot (Y= dissimilarity / X= Conceptual dissimilarity)































