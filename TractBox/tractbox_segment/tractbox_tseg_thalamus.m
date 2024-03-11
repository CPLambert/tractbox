function tractbox_tseg_thalamus(job)
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
if ~isfield(job,'jointdir')
    for k=1:numel(job.leftdir)
        disp(fullfile(['Running Thalamic Tractographic Segmentation: Subject',32,num2str(k),'/',num2str(numel(job.leftdir))]));
        l = cellstr(spm_select('FPListRec',job.leftdir{k},'feature.*.nii'));
        r = cellstr(spm_select('FPListRec',job.rightdir{k},'feature.*.nii'));

        %% Figure out some naming
        [aa,~,~]=fileparts(l{1});
        if contains(aa,'/'),dash=strfind(aa,'/');else,dash=strfind(aa,'\');end
        output = aa(1:dash(end-1));

        %% Create some directories for the data
        tmpdir = fullfile(output,'tmp');mkdir(tmpdir); %Temp
        newop = fullfile(output,options.tseg.thalamus.folder);mkdir(newop); %Output

        %% First lets normalise left and right featuremaps, combine and copy

        % Create an array for later
        TS.FM=cell(size(l));
        TS.INPUT=cell(size(l));
        TS.TMP = cell(size(l));
        TS.MASK = {};

        disp('> Combining Feature Maps');
        TS = combinefmaps(l,r,newop,tmpdir,TS);

        %% For speed/accuracy, the segmentation will proceed as follows:
        % (1) Copy the input data to tmp
        % (2) Align with template space, this will naturally impose the bounding
        % box for faster SHOOT etc., and sort out any space/matric mischeif
        % (3) Delete the superfluous data
        % (4) Create some natives back to individual space.

        %% Lets copy the TPM and mask and coregister
        ttpm{1} = fullfile(output,'tmp','tpm.nii');
        tmsk{1} = fullfile(output,'tmp','mask.nii');
        copyfile(tpm{1},ttpm{1});copyfile(mask{1},tmsk{1});

        %% Note in the original paper we used warps. This should probably work as well but be quicker

        disp('> Segment');
        native2thal(tmsk,TS,options);
        thalseg(TS,ttpm,options);

        %% Now lets clean everything up, and move back to native
        % First, lets store the average orientated segmentations:
        disp('> Cleaning up files in thalamus orientated space');
        [segs,pfile,rpfile] = thalspacecleanup(tmpdir,newop,options);

        disp('> Creating native space results');
        thal2native(TS,tmsk,segs,pfile,rpfile,tmpdir,newop);
    end
else
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

        % Create an array for later

%         TS.INPUT=l;


        %% Lets copy the TPM and mask and coregister
        ttpm{1} = fullfile(output,'tmp','tpm.nii');
        tmsk{1} = fullfile(output,'tmp','mask.nii');
        copyfile(tpm{1},ttpm{1});

        %% Note in the original paper we used warps. This should probably work as well but be quicker

        disp('> Segment');
        thalseg(TS,ttpm,options);

        %% Now lets clean everything up, and move back to native
        % First, lets store the average orientated segmentations:
        disp('> Cleaning up files in thalamus orientated space');
        thalspacecleanup(tmpdir,newop,options);
    end
end
if ~options.tseg.thalamus.debug
    rmdir(tmpdir,'s');
end


disp(fullfile(['Finished',32,datestr(now)]));
end

%% --------------SUB-FUNCTIONS---------------------------------------------
%% (1) Create combined, normalised featuremaps

function TS = combinefmaps(l,r,newop,tmpdir,TS)

for i = 1:numel(l)

    N1 = nifti(l{i});N1i = N1.dat(:,:,:);N1i=N1i./(max(abs(N1i(:)))); %Rescale the data so all in the same range
    N2 = nifti(r{i});N2i = N2.dat(:,:,:);N2i=N2i./(max(abs(N2i(:)))); %Rescale the data so all in the same range

    % For naming
    [~,ft,~] = fileparts(l{i});
    dash=strfind(ft,'_');dish=strfind(ft,'_feature');
    fmaps = fullfile(newop,'featuremaps');
    if ~isfolder(fmaps),mkdir(fmaps);end

    newfile = fullfile(fmaps,[ft(1:dash(1)),'desc-thalamus-norm-combined',ft(dish:end),'.nii']);
    Z=zeros(N1.dat.dim);M=zeros(N1.dat.dim);
    TS.FM{i}=newfile;
    %Create mask, this allows use to take the average of any overlapping
    %voxels that may have marginally slightly different values. This is
    %probably ott as any difference will be extremely small, and probably
    %will only apply to a few voxels, but avoids any risk of introducing
    %an asymmetry bias, plus we want the mask anyway.

    M = (N1i~=0) + (N2i~=0);

    Z=(N1i+N2i)./M;
    N=N1;
    N.dat.fname=newfile;
    N.dat(:,:,:)=Z;
    create(N);

    % Make the featuremap
    newfile = fullfile(tmpdir,[ft(1:dash(1)),'desc-thalamus-norm-combined',ft(dish:end),'.nii']);
    TS.TMP{i} = newfile;copyfile(TS.FM{i},TS.TMP{i});

    newfile = fullfile(tmpdir,['r',ft(1:dash(1)),'desc-thalamus-norm-combined',ft(dish:end),'.nii']);
    TS.INPUT{i} = newfile;

    if i==1 %Make the mask on round 1
        newfile = fullfile(fmaps,[ft(1:dash(1)),'mask-thalamus-combined.nii']);
        TS.MASK{i} = newfile;
        N=N1;
        N.dat.fname=newfile;
        N.dat(:,:,:)=M;
        create(N);
        newfile = fullfile(tmpdir,[ft(1:dash(1)),'mask-thalamus-combined.nii']);
        TS.TMASK{i} = newfile;copyfile(TS.MASK{i},TS.TMASK{i});
        newfile = fullfile(tmpdir,['r',ft(1:dash(1)),'mask-thalamus-combined.nii']);
        TS.RMASK{i} = newfile;
    end
end
end
% -------------------------------------------------------------------------
%% (2) Rigidly move native data to HCP average template (also will impose
% small volume bounding box)

function native2thal(tmsk,TS,options) %Coreg step moving native data to thal-population (HCP population average thalamus)
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = tmsk(1);
matlabbatch{1}.spm.spatial.coreg.estwrite.source = TS.TMASK(1);
matlabbatch{1}.spm.spatial.coreg.estwrite.other = TS.TMP;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = options.tseg.thalamus.coreg.sep;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = options.tseg.thalamus.coreg.tol;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = options.tseg.thalamus.coreg.fwhm;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = options.tseg.thalamus.coreg.interp;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch)
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
    matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm = {fullfile([ttpm{1},',',num2str(i)])};
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
function [segs,pfile,rpfile] = thalspacecleanup(tmpdir,newop,options)
segs=cellstr(spm_select('FPListRec',tmpdir,'^c.*.nii'));segs(1:3)=[];
tmaps = fullfile(newop,'tseg-thalamus-ao');mkdir(tmaps);

N=nifti(segs{1});A=zeros(N.dat.dim);B=A;

for i=1:numel(segs)
    [~,bb,~]=fileparts(segs{i});
    dash=strfind(bb,'_feature');
    newfile=fullfile(tmaps,[bb(4:dash),'space-thalorient_tseg-c',num2str(i),'-',options.tseg.thalamus.labels{i},'.nii']);
    copyfile(segs{i},newfile);

    N=nifti(segs{i});Ni=N.dat(:,:,:);Ni(Ni<0.1)=0;
    A(Ni>B)=i;
    B(Ni>B)=Ni(Ni>B);
end

pfile=fullfile(tmpdir,[bb(4:dash),'space-thalorient_tseg-combined-nuclei.nii']);
N.dat.fname = pfile;
N.dat.dtype='FLOAT32';
N.dat(:,:,:) = A;
create(N);

newfile=fullfile(tmaps,[bb(4:dash),'space-thalorient_tseg-combined-nuclei.nii']);
copyfile(pfile,newfile);

rpfile=fullfile(tmpdir,['r',bb(4:dash),'space-thalorient_tseg-combined-nuclei.nii']);
end

% -------------------------------------------------------------------------
%% (5) Rigidly move thal data to native space and clean up

function thal2native(TS,tmsk,segs,pfile,rpfile,tmpdir,newop)
options = tractbox_defaults;
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = TS.MASK(1);
matlabbatch{1}.spm.spatial.coreg.estwrite.source = tmsk(1);
matlabbatch{1}.spm.spatial.coreg.estwrite.other = [segs;pfile];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = options.tseg.thalamus.coreg.sep;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = options.tseg.thalamus.coreg.tol;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = options.tseg.thalamus.coreg.fwhm;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);

segs=cellstr(spm_select('FPListRec',tmpdir,'^rc.*.nii'));
tmaps = fullfile(newop,'tseg-thalamus-native');mkdir(tmaps);

for i=1:numel(segs)
    [~,bb,~]=fileparts(segs{i});
    dash=strfind(bb,'_feature');
    newfile=fullfile(tmaps,[bb(5:dash),'space-native_tseg-c',num2str(i),'-',options.tseg.thalamus.labels{i},'.nii']);
    copyfile(segs{i},newfile);
end

newfile=fullfile(tmaps,[bb(5:dash),'space-native_tseg-combined-nuclei.nii']);
copyfile(rpfile,newfile);
end
