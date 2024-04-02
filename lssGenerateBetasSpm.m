function images = lssGenerateBetasSpm(subject, spmDir, outDir, includeConditions, settings)
% FORMAT images = lssGenerateBetasSpm(subject, spmDir, outDir, includeConditions, settings)
% This function takes an existing first-level SPM.mat file uses it to
% create one of two possible models: multi-regressor and multi-model.
% The multi-regressor approach estimates a single model with all trials
% represented by individual regressors. The multi-model approach estimates
% a model for each individual trial, setting the first regressor to the
% individual trial and all other regressors to be the same as the original
% model. Beta images are then moved and renamed in a single betas
% directory. The multi-regressor approach is similar to that described in
% Rissman et al. 2004 NI, and the multi-model approach is similar to the
% LS-S approach described in Turner et al. 2012 NI.
% This function is integrated with newLSS_correlation for beta-series
% functional connectivity with the multi-model approach through
% batch_newLSS.
%
%
% Inputs:
% subject:               Subject ID. String.
% spmDir:                Path to folder containing SPM.mat file. String.
% outDir:                Path to output directory, where generated files
%                        will be saved. String.
% ignoreConditions:      Conditions to be ignored. Set to NONE if you do
%                        not want to ignore any conditions. Cell array of
%                        strings.
% settings:              Additional settings. Structure.
% settings.model:        Which model type you wish to run: Rissman beta
%                        series (1) or LS-S multiple models (2). Double.
% settings.overwrite:    Overwrite any pre-existing files (1) or not (0).
%                        Double.
% settings.deleteFiles:  Delete intermediate files (1) or not (0). Double.
%
% Outputs:
% images:                Cell array of 4D images generated  by function in
%                        format images{conds}{sessImages}
%
% Requirements: SPM8, cellstrfind (Matlab function written by Taylor Salo)
%               
%
% Created by Maureen Ritchey 121010
% Modified by Taylor Salo 140806-141216 according to adjustments suggested 
% by Jeanette Mumford. LS-S now matches Turner et al. 2012 NI, where the
% design matrix basically matches the original design matrix (for a given
% block), except for the single trial being evaluated, which gets its own
% regressor. Also now you can ignore multiple conditions. I also added some
% overwrite stuff so that it won't re-run existing subjects if you don't
% want it to.

%% MAIN CODE
% Load pre-existing SPM file containing model information
fprintf('\nLoading previous model for %s:\n%s\n', subject, [spmDir, '\SPM.mat']);
if exist([spmDir, '\SPM.mat'],'file')
    OrigSpm = load([spmDir, '\SPM.mat']);
else
    error('Cannot find SPM.mat file.');
end

if ~exist(outDir, 'dir')
    fprintf('\nCreating directory:\n%s\n', outDir);
    mkdir(outDir)
end

temporaryRootFolder = ['/run/shm/' getenv('USER') '/'];
tempDir = [temporaryRootFolder subject '/'];

% Get model information from SPM file
fprintf('\nGetting model information...\n');
files = OrigSpm.SPM.xY.P;
fprintf('Modeling %i timepoints across %i sessions.\n', size(files, 1), length(OrigSpm.SPM.Sess));
% Make trial directory
betaDir = [outDir '\betas\'];
if ~exist(betaDir, 'dir')
    mkdir(betaDir)
end

% MULTI-MODEL APPROACH
if settings.model == 2
    spm_jobman('initcfg')
    spm('defaults', 'FMRI');
    parpool
%     matlabpool open
    
    if settings.useTempFS
        [nRows, nCols] = size(files);
        index = 1:nCols;
        for iRow = 1:nRows
            index = intersect(index, find(OrigSpm.SPM.xY.P(1, :) == OrigSpm.SPM.xY.P(iRow, :)));
        end

        indexDiff = diff(index)==1;
        consecutiveIndex = find([false, indexDiff]~=[indexDiff, false]);
        lastConsecutive = consecutiveIndex(2:2:end);
        newFiles = cellstr(files);
        
        for iFile = 1:nRows
            oldPath = fileparts(files(iFile, 1:lastConsecutive));
            newFile = strrep(files(iFile, :), oldPath, tempDir);
            newFiles{iFile} = newFile;
            newDataDir = fileparts(newFile);
            if ~exist(newDataDir, 'dir')
                mkdir(newDataDir);
            end
            system(['cp ' files(iFile, 1:end-2) ' ' newFile(1:end-2)]);
        end
    end
    
    % Loop across sessions
    for iSess = 1:length(OrigSpm.SPM.Sess)
        rows = OrigSpm.SPM.Sess(iSess).row;
        if settings.useTempFS
            sessFiles = newFiles(rows);
        else
            sessFiles = files(rows', :);
            sessFiles = cellstr(sessFiles);
        end
        
        % original: 6 motion regressors
        covariates = OrigSpm.SPM.Sess(iSess).C.C;
        % modified: 6 motion regressors + spike regressors
        originalNames = cell(1, length(OrigSpm.SPM.Sess(iSess).U));
        originalOnsets = cell(1, length(OrigSpm.SPM.Sess(iSess).U));
        originalDurations = cell(1, length(OrigSpm.SPM.Sess(iSess).U));
        for jCond = 1:length(OrigSpm.SPM.Sess(iSess).U)
            originalNames{jCond} = OrigSpm.SPM.Sess(iSess).U(jCond).name{1};
            originalOnsets{jCond} = OrigSpm.SPM.Sess(iSess).U(jCond).ons;
            originalDurations{jCond} = OrigSpm.SPM.Sess(iSess).U(jCond).dur;
        end
        includeConditions = originalNames;  %-- added in order to ensure that this runs even for runs without a certain condition
        [lssNames, lssOnsets, lssDurations] = lssMakeVectors(originalNames, originalOnsets, originalDurations, includeConditions);
        
        for jCond = 1:length(includeConditions)
            % ||: OR operator for scalar condition
            if settings.overwrite || ~exist([betaDir '4D_' includeConditions{jCond} '_Sess' sprintf('%03d', iSess) '.nii'], 'file')
                % As long as the current condition is in includeConditions,
                % set up a model for each individual trial.
                parfor kTrial = 1:length(lssOnsets{jCond})
                    singleName = lssNames{jCond}{kTrial}{1};
                    names = lssNames{jCond}{kTrial};
                    onsets = lssOnsets{jCond}{kTrial};
                    durations = lssDurations{jCond}{kTrial};
                    
                    % Make trial directory
                    if settings.useTempFS
                        trialDir = [tempDir 'Sess' sprintf('%03d', iSess) '\' singleName '\'];
                    else
                        trialDir = [betaDir 'Sess' sprintf('%03d', iSess) '\' singleName '\'];
                    end
                    if ~exist(trialDir,'dir')
                        mkdir(trialDir)
                    end

                    % Save regressor onset files
                    regFile = [trialDir 'st_regs.mat'];
                    parsave(regFile, names, onsets, durations);

                    covFile = [trialDir 'st_covs.txt'];
                    dlmwrite(covFile, covariates, '\t');

                    % Create matlabbatch for creating new SPM.mat file
                    matlabbatch = createSpmBatch(trialDir, OrigSpm.SPM);
                    matlabbatch = addSessionToBatch(matlabbatch, 1, sessFiles, regFile, covFile, OrigSpm.SPM);
                    
                    % Run matlabbatch to create new SPM.mat file using SPM batch tools
                    if settings.overwrite || ~exist([trialDir 'beta_0001.img'], 'file')
                        fprintf('\nCreating SPM.mat file:\n%s\n\n', [trialDir 'SPM.mat']);
                        spm_jobman('serial', matlabbatch);
                        runBatches = 1;
                    else
                        runBatches = 0;
                    end

                    if runBatches
                        fprintf('\nEstimating model from SPM.mat file.\n');
                        spmFile = [trialDir 'SPM.mat'];
                        matlabbatch = estimateSpmFile(spmFile);
                        spm_jobman('serial', matlabbatch);

                        % Copy first beta image to beta directory
                        NewSpm = load(spmFile);
                        stopThat = false
                        for mBeta = 1:length(NewSpm.SPM.Vbeta)
                            if ~isempty(strfind(NewSpm.SPM.Vbeta(mBeta).descrip, singleName)) && ~stopThat
                                betaFile = [trialDir NewSpm.SPM.Vbeta(mBeta).fname];
                                stopThat = true
                            end
                        end
                                
%                         if strfind(betaFile, '.img')
%                             [~, betaFileName, ~] = fileparts(betaFile);
%                             system(['cp ' betaFile ' ' betaDir 'Sess' sprintf('%03d', iSess) '_' singleName '.img']);
%                             system(['cp ' trialDir betaFileName '.hdr ' betaDir 'Sess' sprintf('%03d', iSess) '_' singleName '.hdr']);
%                         elseif strfind(betaFile, '.nii')
%                             system(['cp ' betaFile ' ' betaDir 'Sess' sprintf('%03d', iSess) '_' singleName '.nii']);
%                         end

                        % Discard extra files, if desired.
                        if settings.useTempFS
                            system(['rm -rf ' trialDir]);
                        end
                    end
                end
            end
        end
        
        % Make 4D image for each condition of interest in block.
%         for jCond = 1:length(includeConditions)
%             condVols = dir([betaDir 'Sess' sprintf('%03d', iSess) '_' includeConditions{jCond} '*.img']);
%             if isempty(condVols)
%                 condVols = dir([betaDir 'Sess' sprintf('%03d', iSess) '_' includeConditions{jCond} '*.nii']);
%             end
% 
%             cellVols = struct2cell(condVols);
%             cellVols = cellVols(1, :);
%             for kVol = 1:length(cellVols)
%                 cellVols{kVol} = [betaDir cellVols{kVol} ',1'];
%             end
%             images{jCond}{iSess} = [betaDir '4D_' includeConditions{jCond} '_Sess' sprintf('%03d', iSess) '.nii'];
%             matlabbatch{1}.spm.util.cat.name = [betaDir '4D_' includeConditions{jCond} '_Sess' sprintf('%03d', iSess) '.nii'];
%             matlabbatch{1}.spm.util.cat.vols = cellVols';
%             matlabbatch{1}.spm.util.cat.dtype = 0;
%             
%             if settings.overwrite || ~exist([betaDir '4D_' includeConditions{jCond} '_Sess' sprintf('%03d', iSess) '.nii'], 'file')
%                 save([betaDir '3Dto4D_jobfile.mat'], 'matlabbatch');
%                 spm_jobman('run', matlabbatch);
%             else
%                 fprintf('Exists: %s\n', [betaDir '4D_' includeConditions{jCond} '_Sess' sprintf('%03d', iSess) '.nii']);
%             end
%         end
    end
    
    % Delete data folder, if desired.
    if settings.useTempFS
        system(['rm -rf ' tempDir]);
    end
%     matlabpool close

% MULTI-REGRESSOR APPROACH
elseif settings.model == 1
    spmFile = fullfile(outDir, 'SPM.mat');
    counter = 1;

    % Loop across sessions
    for iSess = 1:length(OrigSpm.SPM.Sess)
        rows = OrigSpm.SPM.Sess(iSess).row;
        sessFiles = files(rows', :);
        sessFiles = cellstr(sessFiles);
        covariates = OrigSpm.SPM.Sess(iSess).C.C;

        onsets = {};
        durations = {};
        names = {};

        for jCond = 1:length(OrigSpm.SPM.Sess(iSess).U)
            % Check for special condition names to lump together
            if ~cellstrfind(OrigSpm.SPM.Sess(iSess).U(jCond).name{1}, includeConditions, '')
                onsets = [onsets OrigSpm.SPM.Sess(iSess).U(jCond).ons'];
                durations = [durations OrigSpm.SPM.Sess(iSess).U(jCond).dur'];
                singleName = [OrigSpm.SPM.Sess(iSess).U(jCond).name{1}];
                names = [names singleName];
                counter = counter + 1;
            % Otherwise set up a regressor for each individual trial
            else
                includeConditions{length(includeConditions) + 1} = OrigSpm.SPM.Sess(iSess).U(jCond).name{1};
                for kTrial = 1:length(OrigSpm.SPM.Sess(iSess).U(jCond).ons)
                    onsets = [onsets OrigSpm.SPM.Sess(iSess).U(jCond).ons(kTrial)];
                    durations = [durations OrigSpm.SPM.Sess(iSess).U(jCond).dur(kTrial)];
                    singleName = [OrigSpm.SPM.Sess(iSess).U(jCond).name{1} '_' num2str(kTrial)];
                    names = [names singleName];
                    counter = counter + 1;
                end
            end
        end

        % Save regressor onset files
        if settings.overwrite || ~exist(spmFile, 'file')
            fprintf('Saving regressor onset files for Session %i: %i trials included\n', iSess, length(names));
            regFile = [outDir 'st_regs_session_' num2str(iSess) '.mat'];
            save(regFile, 'names', 'onsets', 'durations');

            % Save covariates (e.g., motion parameters) that were specified
            % in the original model
            covFile = [outDir 'st_covs_session_' num2str(iSess) '.txt'];
            dlmwrite(covFile, covariates, '\t');

            % Create matlabbatch for creating new SPM.mat file
            if iSess == 1
                matlabbatch = createSpmBatch(outDir, OrigSpm.SPM);
            end
            matlabbatch = addSessionToBatch(matlabbatch, iSess, sessFiles, regFile, covFile, OrigSpm.SPM);
        end
    end
    
    % Run matlabbatch to create new SPM.mat file using SPM batch tools
    if settings.overwrite || ~exist(fullfile(outDir, 'SPM.mat'), 'file')
        fprintf('\nCreating SPM.mat file:\n%s\n', fullfile(outDir, 'SPM.mat'));
        spm_jobman('initcfg')
        spm('defaults', 'FMRI');
        spm_jobman('serial', matlabbatch);
        clear matlabbatch

        fprintf('\nEstimating model from SPM.mat file.\n');
        matlabbatch = estimateSpmFile(spmFile);
        spm_jobman('serial', matlabbatch);
    else
        fprintf('Exists: %s\n', [outDir 'SPM.mat']);
    end
    
    clear matlabbatch
    NewSpm = load(spmFile);
    includeConditions = unique(includeConditions);
    for iCond = 1:length(includeConditions)
        counter = 1;
        for jBeta = 1:length(NewSpm.SPM.Vbeta)
            if strfind(NewSpm.SPM.Vbeta(jBeta).descrip, includeConditions{iCond})
                cellVols{counter} = fullfile(NewSpm.SPM.swd, [NewSpm.SPM.Vbeta(jBeta).fname ',1']);
                counter = counter + 1;
            end
        end
        images{iCond}{1} = [betaDir '4D_' includeConditions{iCond} '.nii'];
        matlabbatch{iCond}.spm.util.cat.name = [betaDir '4D_' includeConditions{iCond} '.nii'];
        matlabbatch{iCond}.spm.util.cat.vols = cellVols';
        matlabbatch{iCond}.spm.util.cat.dtype = 0;
        clear cellVols
    end
    if settings.overwrite || ~exist([betaDir '4D_' includeConditions{end} '.nii'], 'file')
        save([betaDir '3Dto4D_jobfile.mat'], 'matlabbatch');
        spm_jobman('run', matlabbatch);
    else
        fprintf('Exists: %s\n', [betaDir '4D_' includeConditions{end} '.nii']);
    end
    clear spmVar matlabbatch newSPM

% If you didn't set settings.model properly.
else
    error('Specify model type as 1 or 2');
end

clear SPM
end

%% SUBFUNCTIONS
function matlabbatch = createSpmBatch(outDir, SPM)
% FORMAT matlabbatch = createSpmBatch(outDir, SPM)
% Subfunction for initializing the matlabbatch structure to create the SPM
% file.
%
%
% 140311 Created by Maureen Ritchey
matlabbatch{1}.spm.stats.fmri_spec.dir = {outDir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = SPM.xBF.UNITS;
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = SPM.xY.RT;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = SPM.xBF.T;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = SPM.xBF.T0;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = SPM.xBF.Volterra;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
if isempty(SPM.xM.VM)
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
else
    matlabbatch{1}.spm.stats.fmri_spec.mask = {SPM.xM.VM.fname};
end
matlabbatch{1}.spm.stats.fmri_spec.cvi = SPM.xVi.form;
end

function matlabbatch = addSessionToBatch(matlabbatch, iSess, sessFiles, regFile, covFile, SPM)
% FORMAT matlabbatch = addSessionToBatch(matlabbatch, iSess, sessFiles, regFile, covFile, SPM)
% Subfunction for adding sessions to the matlabbatch structure.
%
%
% 140311 Created by Maureen Ritchey
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).scans = sessFiles; %fix this
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).multi = {regFile};
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).multi_reg = {covFile};
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).hpf = SPM.xX.K(iSess).HParam;
end

function matlabbatch = estimateSpmFile(spmFile)
% FORMAT matlabbatch = estimateSpmFile(spmFile)
% Subfunction for creating a matlabbatch structure to estimate the SPM
% file.
%
%
% 140311 Created by Maureen Ritchey
matlabbatch{1}.spm.stats.fmri_est.spmmat = {spmFile};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
end

function [lssNames, lssOnsets, lssDurations] = lssMakeVectors(originalNames, originalOnsets, originalDurations, includeConditions)
% FORMAT [lssNames, lssOnsets, lssDurations] = lssMakeVectors(originalNames, originalOnsets, originalDurations, includeConditions)
% Uses SPM-format vectors (in variables corresponding to names, onsets, and
% durations) to create cell arrays of names, onsets, and durations for each
% LS-S model for conditions of interest.
%
%
% input
% originalNames:        Cell array of condition names for a single block.
% originalOnsets:       Cell array of trial onset vectors for a single block.
% originalDurations:    Cell array of trial duration vectors for a single block.
% includeConditions:    A cell array of events we wish to model with the
%                       beta LSS method.
% 
% output
% lssNames:             Cell array of LS-S condition names for a single
%                       block. Format lssNames{includedCondition}{conditionName}
% lssOnsets:            Cell array of LS-S condition onsets for a single
%                       block. Format lssOnsets{includedCondition}{trialVersion}(trial)
% lssDurations:         Cell array of LS-S condition durations for a single
%                       block. Format lssDurations{includedCondition}{trialVersion}(trial)
for iCond = 1:length(includeConditions)
    % Determine where conditions of interest are in vectors.
    % Setdiff reorders conditions, otherConditions must be reordered.
    otherConditionsIdx = ~strcmp(includeConditions{iCond}, originalNames);
    [otherConditions, originalOrder] = setdiff(originalNames, includeConditions{iCond});
    [~, sortedOrder] = sort(originalOrder);
    otherConditions = otherConditions(sortedOrder);
    includeConditionIdx = find(~otherConditionsIdx);
    
    % Check that condition of interest has more than one trial.
    % If condition A only has one trial, you don't need both ConditionA_001
    % and Other_ConditionA, because Other_ConditionA would be empty.
    if ~isempty(setdiff(originalOnsets{includeConditionIdx}, originalOnsets{includeConditionIdx}(1)))
        for jOnset = 1:length(originalOnsets{includeConditionIdx})
            % Create list of condition names
            % (e.g. ConditionA_001, Other_ConditionA, ConditionB, ConditionC, etc.)
            lssNames{iCond}{jOnset} = [{[originalNames{includeConditionIdx} '_' sprintf('%03d', jOnset)]...
                                        ['Other_' originalNames{includeConditionIdx}]}...
                                       otherConditions];
            
            % Single trial
            lssOnsets{iCond}{jOnset}{1} = originalOnsets{includeConditionIdx}(jOnset);
            lssDurations{iCond}{jOnset}{1} = originalDurations{includeConditionIdx}(jOnset);
            
            % Other trials of same condition
            lssOnsets{iCond}{jOnset}{2} = originalOnsets{includeConditionIdx};
            lssOnsets{iCond}{jOnset}{2}(jOnset) = [];
            lssDurations{iCond}{jOnset}{2} = originalDurations{includeConditionIdx};
            lssDurations{iCond}{jOnset}{2}(jOnset) = [];

            % Other conditions
            counter = 3; % A counter adjusts around the skipped condition.
            for kCond = find(otherConditionsIdx)
                lssOnsets{iCond}{jOnset}{counter} = originalOnsets{kCond};
                lssDurations{iCond}{jOnset}{counter} = originalDurations{kCond};
                counter = counter + 1;
            end
        end
    else
        % Single trial
        lssNames{iCond}{1} = [{[originalNames{includeConditionIdx} '_' sprintf('%03d', 1)]} otherConditions];
        lssOnsets{iCond}{1}{1} = originalOnsets{includeConditionIdx}(1);
        lssDurations{iCond}{1}{1} = originalDurations{includeConditionIdx}(1);
        
        % Other conditions
        conditionCounter = 2; % A counter adjusts around the skipped condition.
        for kCond = find(otherConditionsIdx)
            lssOnsets{iCond}{1}{conditionCounter} = originalOnsets{kCond};
            lssDurations{iCond}{1}{conditionCounter} = originalDurations{kCond};
            conditionCounter = conditionCounter + 1;
        end
    end
end
end

function parsave(outFile, names, onsets, durations)
save(outFile, 'names', 'onsets', 'durations');
end
