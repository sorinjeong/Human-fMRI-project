clear all; close all; clc
addpath(genpath('./toolbox'))

%% load data

load('./data4mm.mat')
bold        = data.bold.bold; % bold signal (N time points X voxels)
dim         = data.bold.dim(1:3); % x,y,z size of the 3d fMRI image 
coords      = data.bold.coords; % x,y,z coordinates of each voxel (x,y,z X voxels)
inds        = sub2ind(dim,coords(1,:),coords(2,:),coords(3,:)); % a value for each x,z,y coordinate
runs        = data.bold.runs; % (N time points X runs)
timepoint   = data.bold.timepoint; % (N time points X 1:6 time points)
ntotalvox   = length(coords); % total voxel number
stm         = data.behavior.stim; % stimulus -1:0:1 = S:M:L (trials X stimulus size)
chc         = data.behavior.choice; % choice -1:1 = 'small(left hand)':'large(right hand)' (trial X choice)

%% hyperparameters

radi        = 1; % a searchlight is a cubic with size (2*radius+1)^3
voxthre     = 10; % the minimum number of voxels for training SVM
timepoints  = [4]; % fMRI time points for analyze
zthre       = [30 32]; % the minimum z-coordinate of fMRI image for analyze, [1, 40]
nvox        = (2*radi+1)^3; % the maximum number of voxels within a searchlight
neffectvox  = sum(coords(3,:)>=zthre(1) & coords(3,:)<=zthre(2)); % the number of voxels to be analyzed

isavepath = ['./results/radi' num2str(radi)];
if isempty(dir(isavepath)) 
    mkdir(isavepath)
end

%% searchlight 

% specify time point 

for itp = timepoints
    
    Indtimepoint    = timepoint==itp;
    
    perfmap         = NaN(ntotalvox,1);
    perfmap_run     = NaN(ntotalvox,8);
    tic
    
    cv = 1;
    for iv = 1:ntotalvox
        
        % specify the searchlight center
        
        ix = coords(1,iv);
        iy = coords(2,iv);
        iz = coords(3,iv);
        
        if iz >= zthre(1) && iz <= zthre(2) && (iz-radi)>0 && (iz+radi)<dim(3)
            
            % specify the voxels coordinate of the searchlight
            
            xs = (ix-radi):(ix+radi);
            ys = (iy-radi):(iy+radi);
            zs = (iz-radi):(iz+radi);            
            isl = zeros(dim);
            for i = xs
                for j = ys
                    for k = zs
                        isl(i,j,k) = 1;
                    end
                end
            end            
            iind = sub2ind(dim,find(isl==1));
           % [x(:,1),x(:,2),x(:,3)] = [ind2sub(dim,find(isl==1));]; %%% cube에 들어간 27개 voxel들을 index말고 좌표로 볼 수 있음
            
            % extract BOLD signal of the searchlight
            
            isl2bold    = NaN(size(iind));
           % nvinsl      = length(nvox); % 왜있는지 모르겠다고 지우심
            for i = 1:nvox
                jind = find(iind(i) == inds);
                if ~isempty(jind)
                    isl2bold(i) = jind;
                end
            end
            isl2bold(isnan(isl2bold)) = [];
            
            %data마다 좀 다름
            ibold       = bold(:,isl2bold);
            mbold       = mean(ibold,2); %이상훈교수님 랩에서는 2 deimension마다 mean을 구함. bold signal마다 NaN이 섞여있는 경우가 있기 때문에 미리 mean해서 골라놓는 것이라고 함
            nvoxinsl    = length(ibold(1,:)); %현재 cube에 들어가있는 graymatter 수가 몇개인지
            
            if nvoxinsl > voxthre
                
                % specify BOLD signal of the time point and normalize BOLD
                % signal %내가 관심있는 timepoint(4)만 해당하는 index 뽑아냄
                
                itp_ind     = timepoint==itp; 
                itp_runs    = runs(itp_ind);
                itp_mbold   = mbold(itp_ind);
                itp_bold    = ibold(itp_ind,:);
                itp_zbold   = (itp_bold - repmat(mean(itp_bold,'omitnan'),208,1))./repmat(std(itp_bold,'omitnan'),208,1); %z-score NaN value섞여있어서 수식으로 직접 구함 꼭해줘야함
                
                % SVM train and test on the searchlight (one-run-leave-out)
                
                svms = NaN(200,1);
                for itestrun = 1:8
                    indtest     = itp_runs==itestrun & ~isnan(itp_mbold);
                    indtrain    = itp_runs~=itestrun & ~isnan(itp_mbold);
                    
                    trainbold   = itp_zbold(indtrain,:);
                    testbold    = itp_zbold(indtest,:);
                    trainlabel  = chc(indtrain);
                    testlabel   = chc(indtest);
                    
                    model   = fitcsvm(trainbold,trainlabel,'KernelFunction','linear');
                    isvm    = predict(model, testbold);
                    
                    coef = [model.Beta' model.Bias];
                                        
                    svms(indtest)               = isvm;
                    perfmap_run(iv,itestrun)    = mean(isvm==testlabel,'omitnan');
                end                
                iperf       = mean(svms==chc,'omitnan');
                perfmap(iv) = iperf;      
                
                % real time performance check
                
                timepervox = toc/cv;
                expectedelapse = (neffectvox-cv)*timepervox/60;                
                fprintf('itp=%d | vox=%d/%d | perf=%.2f | remain=%.fmin \n',itp,cv,neffectvox,iperf,expectedelapse)
                
                cv = cv + 1;
            end
            
        end
    end
    
    % save the decoding results at the time point
    
    save([isavepath '/result_' num2str(itp) '.mat'],'perfmap','perfmap_run')
    
    % test significance of the decoding for each searchlight
    
    [~,p,~,tstat] = ttest(perfmap_run',0.5,'tail','right');
    t = tstat.tstat;
    
    pmap = NaN(dim);
    tmap = NaN(dim);
    accmap = NaN(dim);
    
    for iv = 1:ntotalvox
        
        x = coords(1,iv);
        y = coords(2,iv);
        z = coords(3,iv);
        
        pmap(x,y,z) = p(iv)^-1;
        tmap(x,y,z) = t(iv);
        accmap(x,y,z) = perfmap(iv);
        
    end
    
    % save the statistics in NIFTI format
    
    nii0 = load_nii('./wGM.nii');
    nii0.img = pmap;
    save_nii(nii0,[isavepath '/map_p_' num2str(itp) '.nii'])
    nii0.img = tmap;
    save_nii(nii0,[isavepath '/map_t_' num2str(itp) '.nii'])
    nii0.img = accmap;
    save_nii(nii0,[isavepath '/map_acc_' num2str(itp) '.nii'])
    
end