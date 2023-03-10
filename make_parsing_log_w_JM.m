FileList = {'CL121121_1','CL121122_1','CL121128_1','CL121227_1','CL130107_1','CL130109_1','CL130114_2','CL130116_2',...
    'CL130121_2','CL130122_1','CL130130_1','CL130219_1','CL130220_1','CL130225_2','CL130226_1','CL130227_1'};


for fi = 1:numel(FileList)
    filenameo=FileList{fi};
filefolder= ['Y:\EPhysRawData\fmri_oppa_analysis\' filenameo];
addpath(filefolder)

%% Import the file
Root =[ filefolder '\'];
fileToRead1=strcat(filenameo,'B.log');
DELIMITER = ' ';
HEADERLINES = 50000;
no_trials=80;

rawData1 = importdata([Root fileToRead1], DELIMITER, HEADERLINES);
log_data=rawData1;

%% sample LOG file 300 rows
LOG = log_data(1:300,:);
LogStr = string(LOG);

    %% test whether the log contains each event name
    EventName = ["Pause", "causeevent timeframe"];
    VarName = ["Pause", "CTF"];
    coordinate_pat = ["X=", "Y=", "Z="];

    S=struct;
    Others=[];E=[];N=[];
for i = 1:height(timeNevent)
    for j = 1: length(EventName)
    S.(VarName{j}) = []; 

    contain_TF = contains(timeNevent.event{i},optionalPattern(lettersPattern | digitsPattern) + EventName{j} + optionalPattern(lettersPattern | digitsPattern), "IgnoreCase",true);
    if contain_TF == 1
       E= timeNevent.time(i) ;

       digit = extract(timeNevent.event(i),digitsPattern);
       if ~isempty(digit)
       E= [E, digit];
       end

    end
    
    if ~isempty(E)
        E= [E; E];
    elseif contains(timeNevent.event(i),coordinate_pat) ==1
          continue;
    else
        Others = [Others; timeNevent.time(i), timeNevent.event(i)];
    end
    S.Others = Others;
    S.(VarName{j}) = E;
    end


end
