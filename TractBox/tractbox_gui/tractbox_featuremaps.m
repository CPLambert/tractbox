function tractbox_featuremaps(job)
%% function tractbox_featuremaps(job)
% Take SPM batch input and pass to tractbox_create_featureimage(dotimport,options)
% Note I have not added the option to be more selective over map output in the GUI
% at the moment. High-level users can simply modify defaults if this is a
% pressing issue.
%_______________________________________________________________________
% Version History:
% Version 1.0, January 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

options = tractbox_defaults;

%% Update options files if required:
for ii=1:numel(job.featureinput)
    if isfield(job.featureinput{ii},'tractboxfile')
        dotimport=job.featureinput{ii}.tractboxfile;
    end

    if isfield(job.featureinput{ii},'edstepsize')
        options.featuremap.edstep =job.featureinput{ii}.edstepsize;
    end

    if isfield(job.featureinput{ii},'savecoor')
        options.save.correlation =job.featureinput{ii}.savecoor;
    end

    if isfield(job.featureinput{ii},'saveed')
        options.save.ed =job.featureinput{ii}.saveed;
    end

    if isfield(job.featureinput{ii},'threshold')
        options.featuremap.thr =job.featureinput{ii}.threshold;
    end

end

%% Make maps
for i = 1:numel(dotimport)
    tractbox_create_featureimage(dotimport{i},options);
end
end