FileList = {'CL121121_1','CL121122_1','CL121128_1','CL121227_1','CL130107_1','CL130109_1','CL130114_2','CL130116_2',...
    'CL130121_2','CL130122_1','CL130130_1','CL130219_1','CL130220_1','CL130225_2','CL130226_1','CL130227_1'};
% subs = {'17';'53';'1006'; '1007'; '1016'; '2006'; '2007'; '2016'};
subs = [17 53 1006 1007 1016 2006 2007 2016];
for fi = 1:numel(FileList)
    filenameo=FileList{fi};

    filefolder= ['Y:\EPhysRawData\fmri_oppa_analysis\' filenameo];
    addpath(filefolder)
    X=importdata('MR_all.mat');
    Xnew=importdata('MR_seg.mat');
    
size(X)
size(Xnew)
%%
M_total=[];
for sb= 1:length(subs)
   
    linearidx = find(Xnew==subs(sb));
    siz = [size(X,1), size(X,2), size(X,3)];
    [I,J,K] = ind2sub(siz,linearidx);
    idx=[I,J,K];


    A = [];
    for t=1:size(X,4)
        for i=1:size(idx,1)
            x=idx(i,1); y=idx(i,2); z=idx(i,3);
            A(i)=X(x,y,z,t);
        end

        M_total(sb,t)=mean(A);


    end

end
save(['C:\Users\sorin\Documents\MATLAB\processed\' filenameo],"M_total")
end



