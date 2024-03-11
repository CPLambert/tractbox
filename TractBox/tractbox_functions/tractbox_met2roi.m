function tractbox_met2roi(job)
%% function tractbox_met2roi(job)
% Map some measures to binary ROIs
%_______________________________________________________________________
% Version History:
% Version 1.0, February 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

%Data check
if numel(job.met2roiinput.roiimg) ~= numel(job.met2roiinput.metric)
    disp('Error - Uneven inputs - Aborting');
    return
end

for i=1:numel(job.met2roiinput.roiimg)
    disp(char(strcat('Processing Subject:',num2str(i),'/',num2str(numel(job.met2roiinput.roiimg)))));
    N=nifti(job.met2roiinput.roiimg{i});
    Ni=N.dat(:,:,:);
    Z=zeros(N.dat.dim);
    Z(Ni>0)=job.met2roiinput.metric(i);
    [aa,bb,cc] = fileparts(job.met2roiinput.roiimg{i});
    filename = fullfile(aa,[bb,'_met2roi-',job.met2roiinput.targetname,cc]);
    N.dat.fname = filename;
    N.dat(:,:,:) = Z;
    create(N);
end
end
