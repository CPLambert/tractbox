function tractbox_importfsldot(job)
%% function tractbox_importfsldot
% This will take fsl .dot output and import it into a data array for the
% tractbox tools (makes things a little quicker, and starts to fold in
% BIDS)
%_______________________________________________________________________
% Version History:
% Version 1.0, December 2021
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

for i=1:numel(job)

%% Prepare a few things
root=tractbox_import_struct;
rootjson=tractbox_import_json;
cfile=[];tfile=[];sfile=[];

for ii=1:numel(job(i).newdot)
    if isfield(job(i).newdot{ii},'subname')
        root.subject=job(i).newdot{ii}.subname;
    end

        if isfield(job(i).newdot{ii},'seedname')
        root.seed=job(i).newdot{ii}.seedname;
        if contains(root.seed,' ')
            root.seed(strfind(root.seed,32))='-';
        end
        end

        if isfield(job(i).newdot{ii},'fsldot')
        root.dot=job(i).newdot{ii}.fsldot;
        end

        if isfield(job(i).newdot{ii},'output')
        root.output=job(i).newdot{ii}.output;
        end

        if isfield(job(i).newdot{ii},'fsltract')
        tfile=job(i).newdot{ii}.fsltract;
        end

        if isfield(job(i).newdot{ii},'fslcoord')
        cfile=job(i).newdot{ii}.fsltract;
        end

        if isfield(job(i).newdot{ii},'fslseedspace')
        sfile=job(i).newdot{ii}.fslseedspace;
        end
end

if isempty(sfile) || isempty(root.dot) || isempty(root.seed) || isempty(root.subject)
    disp('Error - Not enough imputs supplied');return
end

%% Now fill in the blanks

if isempty(root.output)
[root.output,~,~]=fileparts(root.dot);
end

if isempty(cfile)
cfile=spm_select('FPList',root.output,'^coords_for_.*');
        if isempty(cfile),disp('Error - seed coords not found');return;end
end

if isempty(tfile)
tfile=spm_select('FPList',root.output,'^tract_space_coords_for_.*');
        if isempty(tfile),disp('Error - tract_space_coords not found');return;end
end

%% Now a few basics for making images later
root.seedspace = sfile;N=nifti(sfile);
root.dim = N.dat.dim;
root.mat = N.mat;
root.mat0 = N.mat0;
root.seedcoord = fsl2spmcoordinates(load(cfile),root.dim); %Note to self, in the past this required an offset. This now seems to be fixed
root.tractcoord = fsl2spmcoordinates(load(tfile),root.dim);

%% Now lets save:
filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_tractbox-import']);
save(filename,'root');
filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_tractbox-import.json']);
spm_jsonwrite(filename,rootjson,struct('indent','  '));
clear root rootjson;
end
end



