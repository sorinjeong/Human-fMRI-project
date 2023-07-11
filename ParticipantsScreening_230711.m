%% root 수기지정

cd('Z:\E-Phys Analysis\fMRI_ocat\PilotData_analyzed\ver_230711\');
load('Allsub_NumLogTable.mat');T=total_NumLogTable;clear("total_NumLogTable");
%Session info
SubInfoFile = readtable('OCAT subject info (pilot).xlsx','ReadRowNames',false);
SubInfoFile = renamevars(SubInfoFile(2:19,[1 3 4 5]),["Var1","Var3","Var4","Var5"],["Session","PASS","Sex","Age"]);

%% session별 나누기
Subs = unique(T.Session,"rows","sorted");

SpliT=[];
for i=1:length(Subs)
    varName = sprintf('SUB_%.15g',Subs(i));
    SpliT.(varName) = T(find(T.Session(:) == Subs(i)),:);
% end; clear("varName", "i");


%% Context-Object 조합

C = SpliT.(varName).Context(find(SpliT.(varName).Correct == 1));
O = SpliT.(varName).Object(find(SpliT.(varName).Correct == 1));
A=[C O]


A(find(O==4,1))