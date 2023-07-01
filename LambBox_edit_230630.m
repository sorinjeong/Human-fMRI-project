%% subject numbering , folder root
for subname = 1:1
    Subjects{subname} = sprintf('Sub%.15g', subname);
end
Events=struct;EventName = ["Session","Lap","Trial","Context","Direction","Location","Association","Object","Choice","Correct","RT","isTimeout"];
for v=1:length(EventName); if ~isfield(Events,EventName{v}); Events.(EventName{v}) = [];end;end
for fi = 1:numel(Subjects)
    subjectName=Subjects{fi};
    % Root = ['Z:\E-Phys Analysis\fMRI_ocat\'];
   Root = ['D:\internship\MATLAB\23.06.30_LambBox'];
    % filefolder= [Root 'PilotData'];
    filefolder=Root;
    addpath(filefolder)
    cd([filefolder '\' subjectName])
    % savefolder= [Root 'PilotData_analyzed\'];
    savefolder= [Root '\'];
    addpath(savefolder)

%% Import the file
    fileToRead1= dir([subjectName '_LogTable.mat']);
    fileToRead2= dir([subjectName '_Trial1to32.mat']);
    load(fileToRead1.name);
    load(fileToRead2.name);
   

%% Make LambBox
for v=1:length(EventName)
if 4<v<9
    Events.(EventName{v}) = [Events.(EventName{v}); LogTable{:,v}];
elseif v==9
    
elseif v==10 || v== 11
    Events.(EventName{v}) = [Events.(EventName{v}); LogTable{:,v+2}];
elseif v==12
    Events.(EventName{v}) = zeros(length(LogTable),1);
    for l=1:length(LogTable); if Events(l,10) == 2;Events.isTimeout(l) = 1;end;end
else
     %Session
    subnum=repmat(fi,32,1);
    Events.Session=[Events.Session; subnum];
     %Lap 
    Events.Lap=[Events.Lap; LogTable.Lap];
     %Trial
    Events.Trial=[Events.Trial; Trial1to32]; 
     %Context
    Events.Choice=[Events.Choice; ContextNum];
end
end
%% 
% %save
% Events= orderfields(Events,EventName);
% LambBox=struct2table(Events);
% save([savefolder subjectName '\' subjectName '_LambBox'], "LambBox");
% writetable(LambBox,[savefolder subjectName '\' subjectName '_LambBox.xlsx']);
% 
% 





end