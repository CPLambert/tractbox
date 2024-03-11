function tractbox_fjoin(job)
%% function tractbox_tseg_thalamus(job)
% Create combined normalised feature maps. Doing this seperately is useful
% if there are space issues meaning you will have to manually move the data
% before segment, otherwise use tseg as it combines all of the necessary
% steps.
%_______________________________________________________________________
% Version History:
% Version 1.0, December 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

%% Get defaults, find data
clc
options = tractbox_defaults;
disp(fullfile(['Starting',32,datestr(now)]));
disp(options.tractbox.version);
if ~isfield(job,'jointdir')
    for k=1:numel(job.leftdir)
        disp(fullfile(['Combining feature maps: Subject',32,num2str(k),'/',num2str(numel(job.leftdir))]));
        l = cellstr(spm_select('FPListRec',job.leftdir{k},'feature.*.nii'));
        r = cellstr(spm_select('FPListRec',job.rightdir{k},'feature.*.nii'));

        %% Figure out some naming
        [aa,~,~]=fileparts(l{1});
        if contains(aa,'/'),dash=strfind(aa,'/');else,dash=strfind(aa,'\');end
        output = aa(1:dash(end-1));

        %% Create some directories for the data
        newop = fullfile(output,options.tseg.thalamus.folder);mkdir(newop); %Output


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
            if i==1 %Make the mask on round 1
                newfile = fullfile(fmaps,[ft(1:dash(1)),'mask-thalamus-combined.nii']);
                N=N1;
                N.dat.fname=newfile;
                N.dat(:,:,:)=M;
                create(N);
            end
        end
    end
end
