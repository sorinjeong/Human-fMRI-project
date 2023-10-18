addpath(genpath('C:\fmriLecture\ds000005-download'));

for i=1:16
if i<10
    sbj = strcat("sub-0",string(i));
    ses = strcat("ses-0",string(i));
else
    sbj = strcat("sub-",string(i));
    ses = strcat("ses-",string(i));
end
up_dir = fullfile('C:\fmriLecture\ds000005-download',sbj);


cd(up_dir)
mkdir(fullfile(up_dir,ses))
cd(fullfile(up_dir,ses))
movefile(fullfile(up_dir, 'anat'))
movefile(fullfile(up_dir, 'func'))

end
