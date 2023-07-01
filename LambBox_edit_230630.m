%% subject numbering , folder root
for subname = 1:1
    Subjects{subname} = sprintf('Sub%.15g', subname);
end
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
Events=struct('Session',subnum,'Lap',LogTable.Lap,'Trial',Trial1to32,'Context',ContextNum,'Direction',LogTable.Direction,...
    'Location',LogTable.Location,'Association',LogTable.Association,'Object',LogTable.Obj_ID,'Choice',[],...
    'Correct',LogTable.Decision,'RT',LogTable.Duration,'isTimeout',zeros(height(LogTable),1));
    Events.Choice(find(LogTable.Choice=="A"))=1;
    Events.Choice(find(LogTable.Choice=="B"))=2;
    Events.Choice(find(LogTable.Choice==missing))=missing;
    Events.Choice=Events.Choice';
for l=1:height(LogTable); if Events.Correct(l) == 2;Events.isTimeout(l) = 1;end;end


%% save / for 1 subject
Events= orderfields(Events,EventName);
LambBox=struct2table(Events);
save([savefolder subjectName '\' subjectName '_LambBox'], "LambBox");
writetable(LambBox,[savefolder subjectName '\' subjectName '_LambBox.xlsx']);

%% save /all subjects
total_lambBox=[]; total_lambBox=[total_lambBox, LambBox];
save([savefolder 'AllsubLambBox'],"total_lambBox");
writetable(total_lambBox,[savefolder 'AllsubLambBox.xlsx']);

end