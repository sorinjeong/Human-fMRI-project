FileList = {'CL121121_1','CL121122_1','CL121128_1','CL121227_1','CL130107_1','CL130109_1','CL130114_2','CL130116_2',...
    'CL130121_2','CL130122_1','CL130130_1','CL130219_1','CL130220_1','CL130225_2','CL130226_1','CL130227_1'};

% FileList = {'CL130227_1'};

%%
excelname = 'C:\Users\sorin\Documents\MATLAB\23.04.03_correctness\TotalSubject.xls';
for fi = 1:numel(FileList)
    filenameo=FileList{fi};

    addpath('C:\Users\sorin\Documents\MATLAB\23.04.03_correctness');
    JaeminExcel = readtable (['23.04.03_correctness\' strcat(filenameo,'.xlsx')]);

writetable(JaeminExcel, excelname,'Sheet',filenameo);
end

