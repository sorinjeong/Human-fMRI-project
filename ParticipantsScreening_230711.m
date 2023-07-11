%% root 수기지정

cd('Z:\E-Phys Analysis\fMRI_ocat\PilotData_analyzed\ver_230711\');
load('Allsub_NumLogTable.mat');T=total_NumLogTable;clear("total_NumLogTable");
%Session info
SubInfoFile = readtable('OCAT subject info (pilot).xlsx','ReadRowNames',false);
SubInfoFile = renamevars(SubInfoFile(2:19,[1 3 4 5]),["Var1","Var3","Var4","Var5"],["Session","PASS","Sex","Age"]);

%% session별 나누기
Subs = unique(T.Session,"rows","sorted");

SpliT=[]; To=T; To((find(To.isTimeout==1)),:)=[];
for i=1:length(Subs)
    varName = sprintf('SUB_%.15g',Subs(i));
    SpliT.(varName) = To((To.Session(:) == Subs(i)),:);
% end; clear("varName", "i");


%% Context-Object 조합

C = SpliT.(varName).Context(find(SpliT.(varName).Association == 1));
O = SpliT.(varName).Object(find(SpliT.(varName).Association == 1));
A=[C O];
Asso= sort([A(find(O==4,1)),4; A(find(O==5,1)),5; A(find(O==6,1)),6; A(find(O==7,1)),7]);

N=num2str(Asso(1:4,2)');SubInfoFile.Combi_FFCC=zeros(height(SubInfoFile),1);
N(N ==' ')=[]; Combination = str2num(N);
SubInfoFile.Combi_FFCC(find(SubInfoFile.Session == Subs(i)))=Combination;

end; clear("varName", "i");

%% Data Group (Pass/Fail , Correct/Overall)

DataGroup = struct("PASS",[],"FAIL",[],"Correct",[],"Overall",[]);
DataGroup.Overall = SpliT;
P=SubInfoFile.Session(find(SubInfoFile.PASS == 1));
F=SubInfoFile.Session(find(SubInfoFile.PASS == 0));
for n=1:length(P); DataGroup.PASS.(sprintf('SUB_%.15g', P(n))) = SpliT.(sprintf('SUB_%.15g', P(n))); end
for n=1:length(F); DataGroup.FAIL.(sprintf('SUB_%.15g', F(n))) = SpliT.(sprintf('SUB_%.15g', F(n))); end


Tc=To;Tc((find(Tc.Correct==0)),:)=[];corr_Percent = [];

for i=1:length(Subs)
    varName = sprintf('SUB_%.15g',Subs(i));
    DataGroup.Correct.(varName) = Tc((Tc.Session(:) == Subs(i)),:);

%% Performance-Bias
Screening.(varName)= struct("Accuracy",[],"RT",[],"Bias_Lap",[],"Bias",[]);

ButtonA = length(find(DataGroup.PASS.(varName).Choice==1));
ButtonB = length(find(DataGroup.PASS.(varName).Choice==2));

Screening.(varName).Bias = (ButtonA-ButtonB) / height(DataGroup.PASS.(varName));

Screening.(varName).Bias_Lap = struct;
for j=1:8
lapnum = find(DataGroup.PASS.(varName).Lap==j);

lapnumchoice = DataGroup.PASS.(varName).Choice(lapnum);
ButtonA_Lap = length(find(lapnumchoice==1));
ButtonB_Lap = length(find(lapnumchoice==2));

Screening.(varName).Bias_Lap.(sprintf('Lap%.15g', j)) = (ButtonA_Lap-ButtonB_Lap) / height(lapnumchoice);


%% Performance-Accuracy
corr_Percent = [corr_Percent; ((height(DataGroup.Correct.(varName)))/32)*100];
bias = [bias; Screening.(varName).Bias];


end;end
SubInfoFile = addvars(SubInfoFile, corr_Percent, bias); 

















