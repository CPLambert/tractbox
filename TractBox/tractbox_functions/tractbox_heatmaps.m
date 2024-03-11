function tractbox_heatmaps (job)
%% function tractbox_heatmaps
% Wrapper for SPM GUI to make similarity maps
%_______________________________________________________________________
% Version History:
% Version 1.0, January 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

options = tractbox_defaults;
disp(options.tractbox.version);
disp('Running Heatmap: Calculate')

%Data check
col(1)=numel(job.tractboxfile);
for k = 1:numel(job.regiondata)
    col(k+1)=numel(job.regiondata(k).targetimg);
end

if numel(unique(col))>1
    disp('Error - Unequal Inputs - Aborting');
    return
end

%% Update options files if required:
if numel(unique(col)) == 1
    for i=1:numel(job.tractboxfile)
        disp(char(strcat('Processing Subject:',num2str(i),'/',num2str(numel(job.tractboxfile)))));
        dotimport=job.tractboxfile{i};
        target = cell(numel(col)-1,1);
        label = target;
        for ii = 1:numel(job.regiondata)
            target{ii,1}=job.regiondata(ii).targetimg{i};
            label{ii,1}=job.regiondata(ii).targetname;
        end
        tractbox_similarity_map(dotimport,target,label,options);
    end
end
end