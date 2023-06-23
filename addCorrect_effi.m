% FileList = {'CL121121_1','CL121122_1','CL121128_1','CL121227_1','CL130107_1','CL130109_1','CL130114_2','CL130116_2',...
%     'CL130121_2','CL130122_1','CL130130_1','CL130219_1','CL130220_1','CL130225_2','CL130226_1','CL130227_1'};

FileList = {'CL121122_1'};

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

    datafolder=['D:\internship\MATLAB\23.06.23_OPPA-corr_effi\rawdata'];

    addpath(datafolder);

    ctr_Table = 'TrialInfo_CTRL_B_ori.xlsx';
    ocpr_Table = 'TrialInfo_EXP_B_ori.xlsx';
    
%% control VS ocpr
for ctr_ocpr=1:2
    if ctr_ocpr == 1
       JE = ctr_Table;
       correct_phase1 = ones(length(spm.Sess.U(4).ons),1);

    else
        JE = ocpr_Table;
       correct_phase1 = ones(length(spm.Sess.U(1).ons),1);
    end

%% correctness

Table = readmatrix([datafolder '\' filenameo '\' JE],"Range",'A2:O41');
Table(ismissing(Table))=0;
% JE_correct_ans = Table(:,2);
JE_correct_ph2 = Table(:,5);
JE_correct_ph3 = Table(:,6);
JE_effidx = Table(:,14);
       

%% efficiency and correctness Grouping
Group=[];
for ro = 1:height(Table)
    if Table(ro,13)-Table(ro,12) < 2.8
        Group(ro,1) = 2;
    elseif (JE_effidx(ro) >= mean(JE_effidx)) && (JE_correct_ph3(ro)==1)
        Group(ro,1) = 1;
    else
        Group(ro,1) = 0;
    end
end


       %% add to spm 
       if ctr_ocpr == 1
           spm.Sess.U(4).ons(:,2) = ones(length(spm.Sess.U(4).ons),1);
           spm.Sess.U(5).ons(:,2) = JE_correct_ph2;
           spm.Sess.U(6).ons(:,2) = JE_effidx;
           spm.Sess.U(6).ons(:,3) = JE_correct_ph3;
           spm.Sess.U(6).ons(:,4) = Group;
       else
           spm.Sess.U(1).ons(:,2) = ones(length(spm.Sess.U(1).ons),1);
           spm.Sess.U(2).ons(:,2) = JE_correct_ph2;
           spm.Sess.U(3).ons(:,2) = JE_effidx;
           spm.Sess.U(3).ons(:,3) = JE_correct_ph3;        
           spm.Sess.U(3).ons(:,4) = Group;    
       end
end

%% save
if ~isfolder (['D:\internship\MATLAB\23.06.23_OPPA-corr_effi\' filenameo])
mkdir(['D:\internship\MATLAB\23.06.23_OPPA-corr_effi\' filenameo]); 
end
save(['D:\internship\MATLAB\23.06.23_OPPA-corr_effi\' filenameo '\SPM'], "spm");


end





