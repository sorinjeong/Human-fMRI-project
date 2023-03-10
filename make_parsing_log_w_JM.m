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

%% make time and event table using function
event_pat = "["+digitsPattern(4) + "." + digitsPattern(2) + "] " + optionalPattern(lettersPattern+": ")+optionalPattern(lettersPattern+":  ");

timeNevent = table; event_pre = []; time_pre = [];
for s =1:length(LogStr)
   [time,event] = get_time_event_from_log(LogStr(s));
   time_pre(s) = time;
   event_pre{s} = event;
end
processed_event = squeeze(split(string(event_pre), event_pat));
processed_event = [processed_event(:,2)];

timeNevent.time = time_pre';
timeNevent.event = processed_event;

    %% test whether the log contains each event name
    EventName = ["Pause", "causeevent timeframe", "Pause_num", "causeevent timeframe_num"];
    VarName = ["Pause", "CTF", "Pause_num", "CTF_num"];
    coordinate_pat = ["X=", "Y=", "Z="];
    % other_event_name = ["Others", "Others_note"]

    % V=[]; V1=[];V2=[];V_name1=[]; V_name2=[]; Others=[];
    % V=table(EventName); V=splitvars(V);
    S=struct;
    Others=[];E=[];N=[];

    for j = 1: length(EventName)
        S.(VarName{j}) = [];
%         message_num = strcat(EventName{j},'_num');
%         index_message_num = find(contains(EventName,message_num));

        for i = 1:height(timeNevent)
            

            %         contain_TF = strcmpi(timeNevent.event{i},char(event_name(j)));
contain_TF = contains(timeNevent.event{i},optionalPattern(lettersPattern | digitsPattern) + EventName{j} + optionalPattern(lettersPattern | digitsPattern), "IgnoreCase",true);
            if contain_TF == 1
                E= [E; timeNevent.time(i)] ;
                %            V(i,"EventName_"+j) = timeNevent.time(i);

                digit = extract(timeNevent.event(i),digitsPattern);
                if ~isempty(digit)
                    %              V(i,"EventName_"+index_message_num) = digit;
                 E= [E; timeNevent.time(i), digit];
                end

            elseif contains(timeNevent.event(i),coordinate_pat) ==1
                continue;

            else

                %         OthersTime=[OthersTime; timeNevent.time(i)];
                %         OthersNote=[OthersNote; timeNevent.event(i)];
                Others = [Others; timeNevent.time(i), timeNevent.event(i)];
                %             V.OthersTime(i) = timeNevent.time(i);
                %             V.OthersNote(i) = timeNevent.event(i);
            end



        end
        S.(VarName{j}) = E
%         S.(VarName{index_message_num}) = N
        S.Others = Others


    end
end

