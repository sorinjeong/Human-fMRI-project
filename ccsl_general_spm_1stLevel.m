function [] = ccsl_general_spm_1stLevel(root_path, glm_name, subjID, session, events, contrasts, TR, deriv, mthresh, smoothing, skip_smoothing, remove_unzipped, movereg_names, task_name, num_dummy)
% formerly...(root_path, subjID, session, events, smoothing, TR, glm_name, num_runs, mthresh)
% Matlab code for running first-level analysis (BIDS data), multiple
% sessions, for any task/study...!
% :Inputs:
% 1. root_path: root path of a dataset.
% fmriprepped files will be in .../root_path/derivatives/fmriprep
% event files will be in
% .../root_path/sub-xxx/ses-yy/func/
% 2. glm_name: glm name (e.g., glm_name = 'glm2')
% 3. subjID: subject ID (e.g., subjID = 'sub-001')
% 4. session: session (visit) number (e.g., session = 1)
% 5. events: all event info should be saved into events
% 6. contrasts: all contrasts to be used in this study (only t-test
%    contrasts are supported yet)
% 7. TR: TR (e.g., TR = 1)
% 8. deriv: temporal & dispersion derivatives (e.g., deriv = [0 0])
% 9. mtresh: treshold for brain acvitation (e.g., mtresh = 0.8 #SPM default)
% 10. smoothing: spatial smoothing in mm (e.g., smoothing = [8 8 8])
% 11. skip_smoothing: skip smoothing if ssub-xxx file exists? (e.g.,
% skip_smoothing = 1)
% 12. skip_unzip: skip unzipping if *preproc_bold.nii file exists? (e.g.,
% skip_unzip = 1)
% 13. movereg_names: column names of 6 movement regressors in the tsv file

% set defaults
if ~exist('root_path','var'), error('root_path is missing!'); end
if ~exist('glm_name','var'), error('glm_name is missing!'); end
if ~exist('subjID','var'), error('subjID is missing!'); end
if ~exist('session','var'), error('session number is missing!'); end
if ~exist('events','var'), error('events is missing!'); end
if ~exist('contrasts','var'), error('contrasts is missing!'); end
if ~exist('TR','var'), error('TR is missing!'); end
if ~exist('deriv','var'), deriv=[0 0]; disp('Use default: deriv=[0 0]'); end
if ~exist('mthresh','var'), mthresh=0.2; disp('Use default: mthresh=0.2'); end
if ~exist('skip_smoothing','var'), skip_smoothing=1; disp('Use default: skip smoothing if already done'); end
if ~exist('remove_unzipped','var'), remove_unzipped=0; disp('Use default: unzipped nii files are not removed'); end
% if ~exist('skip_unzip','var'), skip_unzip=0; disp('Use default: do not skip unzipping even if already done'); end
if ~exist('movereg_names','var'), movereg_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'}; disp('Use default names for movement regressors: trans_* & rot_*'); end
if ~exist('task_name','var'), error('TR is missing!'); end
if ~exist('num_dummy','var'), num_dummy=0; disp('Use default: no dummy scans to exclude'); end


%% set ID, def_path
%defThres = 0.2;   % default threshold
%stimDuration = 0;   % stimulus duration (0 sec = delta function);
%smoothing = [8 8 8];  % smoothing
%currApproach = ['tom2007_d' num2str(stimDuration)];  % current approach..

%[data_path ID] = fileparts(pwd);  % e.g., ID = '576-D-1';
%data_path = '/mnt/nfs/proj/visuperc/wahn/data/';  % CHECK THIS PATH!!
%def_path = fullfile(data_path, ID);

disp( ['ID = ' subjID ] );
disp( ['pwd = ' pwd ])


%% check number of runs (num_runs) from 'events'
num_runs = length(events);

%% Path containing data
% path for fmriprep data 
if isnumeric(session)
    nii_path = fullfile(root_path, 'derivatives', 'fmriprep', subjID,  ['ses-' num2str(session, '%02.f')], 'func'); % assume nii.gz file is in '.../func'
else
    nii_path = fullfile(root_path, 'derivatives', 'fmriprep', subjID,  ['ses-' session], 'func'); % assume nii.gz file is in '.../func'
end

% path for confounding factors
conf_file = dir(fullfile(nii_path, ['*task-', task_name, '*confounds_timeseries.tsv']));   % find the file

%% create "R" variable from movement_regressor matrix and save
for rIdx = 1:num_runs
    [data1, header1, ] = tsvread(fullfile(nii_path, conf_file(rIdx).name));
    
    trans_x=strmatch(movereg_names{1},header1,'exact');
    trans_y=strmatch(movereg_names{2},header1,'exact');
    trans_z=strmatch(movereg_names{3},header1,'exact');
    rot_x=strmatch(movereg_names{4},header1, 'exact');
    rot_y=strmatch(movereg_names{5},header1,'exact');
    rot_z=strmatch(movereg_names{6},header1,'exact');
    
    % remove the first row, 26-31 columns --> movement regressors    R_mov1 = fillmissing(R_mov1, ?nearest?);    R_mov1(~isfinite(R_mov1))=0;
%     R = data1(2:end, [trans_x,trans_y,trans_z,rot_x,rot_y,rot_z]);
    R = data1(num_dummy+2:end, [trans_x,trans_y,trans_z,rot_x,rot_y,rot_z]); % remove dummy scan
    R = fillmissing(R, 'nearest');
    R(~isfinite(R)) = 0;
    
    %R = data1(2:end, (end-5):end);  % remove the first row, 26-31 columns --> movement regressors
    tmp_move_file_path = fullfile(nii_path, ['movement_regressors_run' num2str(rIdx) '.mat']);
    save(tmp_move_file_path, 'R')
    move_path_run{rIdx} = tmp_move_file_path;
end
disp('confound variables movement path defined')

%% Initialise SPM defaults
spm('defaults', 'FMRI');
spm_jobman('initcfg');

%%
%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------

% create a directory where data will be saved
if isnumeric(session) 
    save_path = fullfile(root_path, 'derivatives', glm_name, 'lev-1st', subjID,  ['ses-' num2str(session, '%02.f')]);
else
    save_path = fullfile(root_path, 'derivatives', glm_name, 'lev-1st', subjID,  ['ses-' session]);
end

mkdir(save_path);

% delete SPM.mat file if it exists already
if exist( fullfile( save_path, 'SPM.mat') )
    fprintf('\n SPM.mat exists in this directory. Overwriting SPM.mat file! \n\n')
    delete( fullfile( save_path, 'SPM.mat') )
end

%% gunzip all nii.gz files & Smoothing 
%% 1) check if smoothed *.nii files exist
nii_smoothed_files = dir(fullfile(nii_path, ['s' subjID '*task-' task_name '*preproc_bold.nii']));

if length(nii_smoothed_files) == num_runs && skip_smoothing == 1 %IF smoothe file already exist, 
    disp([num2str(num_runs) ' smoothed files already exist. Zipping and Smoothing will be skipped'])
    % Caution! If some nii.gz are unzipped but some are not, this may cause
    % an error! 
%% 2) If not, unzip nii.gz files
else 
    disp('Unzipping all functional image files for SPM analysis...')
    gunzip(fullfile(nii_path, [subjID '*preproc_bold.nii.gz']))
    disp('All functional image files are unzipped for SPM analysis')
    
%% 3) Run smoothing using the unzipped files 
    % These lines show how to convert a 4D file into multiple 3D files
    % And then smooth each 3D files...
    % But you can just smooth a single 4D file with SPM12..!
    for rIdx = 1:num_runs
        matlabbatch = [];  % clear matlabbatch..
        scanFiles = []; % clear scanFiles
        %epipath = fullfile(nii_path, 'func');  % location of the preprocessed files
        %tmpFiles = dir(fullfile(nii_path, 'sub-*run-01*preproc_bold.nii'));   % find the file
        tmpFiles = dir(fullfile(nii_path, ['sub-*task-' task_name '*run-' num2str(rIdx) '*preproc_bold.nii']));   % find the file
        % Here lines differ for 3D vs. 4D %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This is for 4D multiple files
        % get header information to read a 4D file
        tmpHdr = spm_vol( fullfile(nii_path, tmpFiles.name) );
        f_list_length = size(tmpHdr, 1);  % number of 3d volumes
        for jx = 1:f_list_length
            scanFiles{jx,1} = [nii_path '/' tmpFiles.name ',' num2str(jx) ] ; % add numbers in the end
            % End of difference for 3D vs. 4D %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        scanFiles = scanFiles(num_dummy+1:length(scanFiles),:); % remove dummy scans
        
        matlabbatch{1}.spm.spatial.smooth.data = scanFiles;
        matlabbatch{1}.spm.spatial.smooth.fwhm = smoothing;
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run', matlabbatch)
        disp(['run ' num2str(rIdx) ' smoothing is complete'])
                
        matlabbatch = [];  % clear matlabbatch..
    end
%% 4) Remove unzipped files 
    if remove_unzipped == 1 
        delete(fullfile(nii_path, [subjID '*task-' task_name '*preproc_bold.nii']))
        disp(['removed ' num2str(num_runs) ' nii files'])
    else 
        disp(['SKIP removing ' num2str(num_runs) ' nii files'])
    end
end

%% Specify design matrix
matlabbatch{1}.spm.stats.fmri_spec.dir = { save_path };
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
%%
for rIdx=1:num_runs
    scanFiles = [];  % reset it
    % rescan files
    % assume file name contains 'ssub*run-*1_*preproc_bold.nii'
    tmpFiles = dir(fullfile(nii_path, ['ssub-*task-' task_name '*run-' num2str(rIdx) '*preproc_bold.nii']));   % find the file
    % Here lines differ for 3D vs. 4D %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is for 4D multiple files
    % get herder information to read a 4D file
    tmpHdr = spm_vol( fullfile(nii_path, tmpFiles.name) );
    f_list_length = size(tmpHdr, 1);  % number of 3d volumes
    for jx = 1:f_list_length
        scanFiles{jx,1} = [nii_path '/' tmpFiles.name ',' num2str(jx) ] ; % add numbers in the end
        % End of difference for 3D vs. 4D %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    scanFiles = scanFiles(num_dummy+1:length(scanFiles),:); %remove dummy scans
    matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).scans = scanFiles;
    
    
    % use 'events' 
    num_conds = length(events{rIdx});  % number of conditions
    for cond_idx = 1:num_conds
        matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).cond(cond_idx).name = events{rIdx}{cond_idx}.name;
        matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).cond(cond_idx).onset = events{rIdx}{cond_idx}.onset;
        matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).cond(cond_idx).duration = events{rIdx}{cond_idx}.duration;
        matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).cond(cond_idx).tmod = 0; % use SPM default
        
        % check parametric modulators
        num_pmods = length(events{rIdx}{cond_idx}.pmod); % # of PMs at this event onset
        if num_pmods > 0
            for pm_idx = 1:num_pmods
                matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).cond(cond_idx).pmod(pm_idx).name = events{rIdx}{cond_idx}.pmod(pm_idx).name;
                matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).cond(cond_idx).pmod(pm_idx).param = events{rIdx}{cond_idx}.pmod(pm_idx).param;
                matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).cond(cond_idx).pmod(pm_idx).poly = 1;    
            end
        else
            % empty cells for PMs (pmod)
            matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).cond(cond_idx).pmod = struct('name', {}, 'param', {}, 'poly', {});
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).cond(cond_idx).orth = 0;   % don't orthogonalize PM regressors
    end   
    % Remaining details...
    matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).multi_reg = move_path_run(rIdx); %{move_path_run1};
    matlabbatch{1}.spm.stats.fmri_spec.sess(rIdx).hpf = 128;
end

%%
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = deriv; %[0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = mthresh;   % threshold
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

%% run categorical model specification
disp([subjID ' - specifying categorial model']);

spm_jobman('run', matlabbatch)
disp('categorical model is specified')

%% categorical model estimation
matlabbatch = [];
matlabbatch{1}.spm.stats.fmri_est.spmmat = { fullfile(save_path, 'SPM.mat') };
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
disp([subjID ' - estimating categorial model']);

spm_jobman('run', matlabbatch)
disp('categorical model is estimated')

%% create contrasts
matlabbatch = [];
matlabbatch{1}.spm.stats.con.spmmat = { fullfile(save_path, 'SPM.mat') };
disp([subjID ' - generating contrast']);
% check the number of contrasts
num_contrasts = length(contrasts);
% construct batch editor w/ 'contrasts' variable
for contIdx = 1:num_contrasts
    matlabbatch{1}.spm.stats.con.consess{contIdx}.tcon.name = contrasts{contIdx}.name;
    matlabbatch{1}.spm.stats.con.consess{contIdx}.tcon.convec = contrasts{contIdx}.convec;
    matlabbatch{1}.spm.stats.con.consess{contIdx}.tcon.sessrep = 'none';
end
matlabbatch{1}.spm.stats.con.delete = 1;  % delete existing contrasts? yes
spm_jobman('run', matlabbatch)
disp([glm_name ': contrasts are generated'])

end

