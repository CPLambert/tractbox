function tractbox_tseg2_thalamus(job)
%% function tractbox_tseg_thalamus(job)
% Thalamic segmentation as described in Lambert et al., 2017 Defining thalamic
% nuclei and topographic connectivity gradients in vivo. Needs the
% feature-maps to be pre-caclulated. Will produce small volume
% segmentations, rigidly aligned to thalamic population average space (i.e.
% for warping/VBM) but will also move the outputs back to individual
% subject native space. Will also create a maximum likelihood image that
% combines all of the nuclei.
%_______________________________________________________________________
% Version History:
% Version 1.0, December 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

%% Get defaults, find data
clc
options = tractbox_defaults;
mask = cellstr(options.tseg.thalamus.mask);
tpm = cellstr(options.tseg.thalamus.tpm);
disp(fullfile(['Starting',32,datestr(now)]));
disp(options.tractbox.version);

    for k=1:numel(job.jointdir)
        l = cellstr(spm_select('FPListRec',job.jointdir{k},'feature.*.nii'));
        %% Figure out some naming
        [aa,~,~]=fileparts(l{1});
        if contains(aa,'/'),dash=strfind(aa,'/');else,dash=strfind(aa,'\');end
        output = aa(1:dash(end-1));

        %% Create some directories for the data
        tmpdir = fullfile(output,'tmp');mkdir(tmpdir); %Temp
        newop = fullfile(output,options.tseg.thalamus.folder);mkdir(newop); %Output

        %% move FMs
        for ii=1:numel(l)
            [aa,bb,cc]=fileparts(l{ii});
            newfile = fullfile(tmpdir,[bb,cc]);
            copyfile(l{ii},newfile);
            TS.INPUT{ii}=newfile;
        end


        %% Lets copy the TPM and mask and coregister
        ttpm{1} = fullfile(output,'tmp','tpm.nii');
          ttpm{2} = fullfile(output,'tmp','wtpm.nii');
        copyfile(tpm{1},ttpm{1});

        tmsk{1} = fullfile(output,'tmp','mask.nii');
         tmsk{2} = fullfile(output,'tmp','wmask.nii');
        copyfile(mask{1},tmsk{1});
        warpthal2native(TS.INPUT{1},ttpm,tmsk,job.thwarp(k),tmpdir) 
        %% Note in the original paper we used warps. This should probably work as well but be quicker

        disp('> Segment');
        thalseg(TS,ttpm,options);

        %% Now lets clean everything up, and move back to native
        % First, lets store the average orientated segmentations:
        disp('> Cleaning up files in thalamus orientated space');
        thalspacecleanup(tmpdir,newop,tmsk{2},options);
    end

if ~options.tseg.thalamus.debug
    rmdir(tmpdir,'s');
end
disp(fullfile(['Finished',32,datestr(now)]));
end


% -------------------------------------------------------------------------
%% (2) Rigidly move native data to HCP average template (also will impose
% small volume bounding box)

function warpthal2native(spac,ttpm,tmsk,tdef,op) %Coreg step moving native data to thal-population (HCP population average thalamus)
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = tdef;
matlabbatch{1}.spm.util.defs.comp{1}.inv.space ={spac};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = [ttpm(1);tmsk(1)];
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {op};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 0;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
spm_jobman('run',matlabbatch);
end

% -------------------------------------------------------------------------
%% (3) Use the previously calculate thalamus TPMs in SPM segment as per the
% original paper:

function thalseg(TS,ttpm,options)
% Lets now set up the segment batch
for i = 1:numel(TS.INPUT)
    matlabbatch{1}.spm.spatial.preproc.channel(i).vols = TS.INPUT(i);
    matlabbatch{1}.spm.spatial.preproc.channel(i).biasreg = options.tseg.thalamus.seg.biasreg;
    matlabbatch{1}.spm.spatial.preproc.channel(i).biasfwhm = options.tseg.thalamus.seg.biasfwhm;
    matlabbatch{1}.spm.spatial.preproc.channel(i).write = options.tseg.thalamus.seg.write;
end

for i = 1:12
    matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm = {fullfile([ttpm{2},',',num2str(i)])};
    matlabbatch{1}.spm.spatial.preproc.tissue(i).ngaus = options.tseg.thalamus.seg.ngaus ;
    matlabbatch{1}.spm.spatial.preproc.tissue(i).native = options.tseg.thalamus.seg.native;
    matlabbatch{1}.spm.spatial.preproc.tissue(i).warped = options.tseg.thalamus.seg.warped;
end

matlabbatch{1}.spm.spatial.preproc.warp.mrf = options.tseg.thalamus.seg.mrf;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = options.tseg.thalamus.seg.cleanup;
matlabbatch{1}.spm.spatial.preproc.warp.reg = options.tseg.thalamus.seg.reg ;
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = options.tseg.thalamus.seg.fwhm;
matlabbatch{1}.spm.spatial.preproc.warp.samp =  options.tseg.thalamus.seg.samp;
matlabbatch{1}.spm.spatial.preproc.warp.write = options.tseg.thalamus.seg.write;
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
    NaN NaN NaN];
spm_jobman('run',matlabbatch);
end

% -------------------------------------------------------------------------
% (4) Clean-up and move the thalamic space data:
function [segs,pfile,rpfile] = thalspacecleanup(tmpdir,newop,mask,options)
segs=cellstr(spm_select('FPListRec',tmpdir,'^c.*.nii'));segs(1:3)=[];
tmaps = fullfile(newop,'tseg-thalamus-ao');mkdir(tmaps);

N=nifti(segs{1});A=zeros(N.dat.dim);B=A;
M=nifti(mask);Mi=M.dat(:,:,:);Mi = double(Mi>0.5);

for i=1:numel(segs)
    [~,bb,~]=fileparts(segs{i});
    dash=strfind(bb,'_feature');
    newfile=fullfile(tmaps,[bb(3:dash),'space-thalorient_tseg-c',num2str(i),'-',options.tseg.thalamus.labels{i},'.nii']);
  

    N=nifti(segs{i});
    Ni=N.dat(:,:,:).*Mi;
    Ni(isnan(Ni))=0;
    Ni(Ni<0.1)=0;
%     Z=zeros(size(Ni));
%     Z(Mi==1)=Ni(Mi==1);
    N.dat(:,:,:)=Ni;
    create(N);
    A(Ni>B)=i;
    B(Ni>B)=Ni(Ni>B);
      copyfile(segs{i},newfile);
end

pfile=fullfile(tmpdir,[bb(3:dash),'space-thalorient_tseg-combined-nuclei.nii']);
N.dat.fname = pfile;
N.dat.dtype='FLOAT32';
N.dat(:,:,:) = A;
create(N);

newfile=fullfile(tmaps,[bb(3:dash),'space-thalorient_tseg-combined-nuclei.nii']);
copyfile(pfile,newfile);

rpfile=fullfile(tmpdir,['r',bb(3:dash),'space-thalorient_tseg-combined-nuclei.nii']);
end

% -------------------------------------------------------------------------

