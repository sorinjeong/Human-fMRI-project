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
event_name = ["Pause", "causeevent_timeframe"];
coordinate_pat = ["X=", "Y=", "Z="];

for j = 1: length(event_name)
message_num = char(strcat(event_name(j),'_num'));
index_message_num = find(contains(event_name,message_num));

    for i = 1:height(timeNevent)
        contain_TF = strcmp(timeNevent.event(i),event_name(j))
        if contain_TF == 1;
           var(j) = timeNevent.time(i) ;
           event_name{j} = [];
           event_name{j} = [event_name{j}; var(j)];

           if ~isempty(extract(timeNevent.event(i),digitsPattern))
               var(index_message_num) = extract(timeNevent.event(i),digitsPattern);
               event_name{index_message_num} = [];
               event_name{index_message_num} = [event_name{index_message_num}; var(index_message_num)];
           end
        elseif contains(timeNevent.event(i),coordinate_pat) ==1;
            continue;

        else
            var(Others) = timeNevent.time(i);
            var(Others_note) = timeNevent.event(i);
        end
    end
end

