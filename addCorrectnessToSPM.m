% FileList = {'CL121121_1','CL121122_1','CL121128_1','CL121227_1','CL130107_1','CL130109_1','CL130114_2','CL130116_2',...
%     'CL130121_2','CL130122_1','CL130130_1','CL130219_1','CL130220_1','CL130225_2','CL130226_1','CL130227_1'};

FileList = {'CL130227_1'};

%%
for fi = 1:numel(FileList)
    filenameo=FileList{fi};
    filefolder= ['Y:\EPhysRawData\fmri_oppa_analysis\' filenameo '\fmri\spm\spm7_concatenated'];
    addpath(filefolder)

    %% Import the file
    Root =[ filefolder '\'];
    fileToRead1='SPM.mat';
    SPM = importdata([Root fileToRead1]);
    spm = SPM;

    addpath('C:\Users\sorin\Documents\MATLAB\23.04.03_correctness');
    JaeminExcel = readmatrix (['23.04.03_correctness\' strcat(filenameo,'.xlsx')]);


%% variable
ctr_correct_phase1 = ones(length(spm.Sess.U(4).ons),1);
ocpr_correct_phase1 = ones(length(spm.Sess.U(1).ons),1);

ctr_correct_phase2=[];ocpr_correct_phase2=[];ctr_correct_phase3=[];ocpr_correct_phase3=[];
for i = 1:2:length(JaeminExcel)
    ctr_correct_phase2 = [ctr_correct_phase2; JaeminExcel(i,6)];
    ctr_correct_phase3 = [ctr_correct_phase3; JaeminExcel(i,7)];
end

    for j = 2:2: length(JaeminExcel)
        if ismissing(JaeminExcel(j,6)) == 1
            ocpr_correct_phase2 = ocpr_correct_phase2;
            ocpr_correct_phase3 = ocpr_correct_phase3;
        else
            ocpr_correct_phase2 = [ocpr_correct_phase2; JaeminExcel(j,6)];
            ocpr_correct_phase3 = [ocpr_correct_phase3; JaeminExcel(j,7)];
        
        end
    end



%% add to spm 
spm.Sess.U(1).ons(:,2) = ocpr_correct_phase1;
spm.Sess.U(2).ons(:,2) = ocpr_correct_phase2;
spm.Sess.U(3).ons(:,2) = ocpr_correct_phase3;
spm.Sess.U(4).ons(:,2) = ctr_correct_phase1;
spm.Sess.U(5).ons(:,2) = ctr_correct_phase2;
spm.Sess.U(6).ons(:,2) = ctr_correct_phase3;

%% save
if ~isfolder (['C:\Users\sorin\Documents\MATLAB\23.04.03_correctness\' filenameo])
mkdir(['C:\Users\sorin\Documents\MATLAB\23.04.03_correctness\' filenameo]); 
end
save(['C:\Users\sorin\Documents\MATLAB\23.04.03_correctness\' filenameo '\SPM'], "spm");


end





