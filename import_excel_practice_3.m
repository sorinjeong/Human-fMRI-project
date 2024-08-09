addpath(genpath('C:\Users\dutls\Documents\LeeLab'))
cd('C:\Users\dutls\Documents\LeeLab')
%%
sheet_names=sheetnames("event_TR.xlsx");
numsheets=numel(sheet_names);
for i= 1: numsheets
c_sbj=sprintf('sub-%02d',i);
curr_T= readtable("event_TR.xlsx","Sheet",c_sbj);
%%
first_scan=find(curr_T.Var2==0,1,'last');
curr_T.scan_num=[zeros(first_scan-1,1);(1:height(curr_T)-first_scan+1)'];
writetable(curr_T,"event_TR.xlsx","Sheet",c_sbj)
end