function tractbox_constrength(job)
%% function tractbox_constrength
% This creates some measure of connectivity strength based on tractography
% images and some parcellation. It assumes the input is appropriately
% scaled (note FSL fdt_paths is not appropriate if you are doing an ROI>1
% voxel because this is a summative image, and therefore meaningless).
% Advise voxel-wise tractography storing -omatrix2 output, and used tractbox
% for secondary processing steps (as it will also correct the slight
% discrepancy in voxel position caused by FSLs counting system compared to
% MATLAB/SPM)
%
% The fundamental problems using tractography as some pseudo-quantitative
% metric of interest is the distance-volume confound with tractography.
% Namely, probability of connectivity decreases naturally with distance
% (and is non-trivial to correct due to natural spatial variability in
% myeloarchitecture) and higher chances of connecting with bigger
% structures. Users should also note tractography is a dumb algorithm, if
% there is a disease process causing structural connectivity to fall in one
% region, then each tracking initiation will have to go somewhere, and therefore you
% may see artificial "increases" in connectivity that reflects this
% process. This is not necessarily a problem, but important to be aware of.
%
% Currently, the constrength code implements a step I have re-labelled as
% "Adaptive Ratio of Connectivity" (adcon), as it didn't seem to have a
% clearly designated term and it best reflects the process (in my mind).
% As far as I'm aware, the steps to calculate this metric was first described in:
%
% Aron AR, Behrens TE, Smith S, Frank MJ, Poldrack RA. Triangulating a
% cognitive control network using diffusion-weighted magnetic resonance imaging (MRI)
% and functional MRI. J Neurosci. 2007, 27(14):3743-52.
%
% There certainly are other approaches, and no one clear method. I will
% implement other metrics if requested
%_______________________________________________________________________
% Version History:
% Version 1.0, November 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------
options = tractbox_defaults;

disp(options.tractbox.version);
disp('Running Connectivity Strength: Calculate')

tbin=[];tractimg=[];cstype=[];csthr=[];parcin=[];outdir=[];
% Harvest inputs and allocate

for i=1:numel(job.csin)
    if isfield(job.csin{i},'tractboxfile')
        tbin = job.csin{i}.tractboxfile;
    end

    if isfield(job.csin{i},'seedname')
        seedname = regexprep(lower(job.csin{i}.seedname), ' +', '-');
    end

    if isfield(job.csin{i},'tractimg')
        tractimg = job.csin{i}.tractimg;
    end

    if isfield(job.csin{i},'cstype')
        cstype = job.csin{i}.cstype;
    end

    if isfield(job.csin{i},'csthr')
        csthr = job.csin{i}.csthr;
    end

    if isfield(job.csin{i},'parcin')
        parcin{end+1} = job.csin{i}.parcin;
    end

    if isfield(job.csin{i},'outdir')
        outdir = job.csin{i}.outdir;
    end
end

if isempty(tbin) %Using individual level con matricies is a special case, and slow (see later code)

    %% Fast version - Group level

    for i = 1:numel(parcin) %loop through all parcellations
        source = [];
        % Parse the parcin options
        for ii = 1:numel(parcin{i})
            fz=cell2mat(fields(parcin{i}{ii}));
            source.(fz) = parcin{i}{ii}.(fz);
        end

        N = nifti(source.parcimg);
        Ni = N.dat(:,:,:);

        % Pre-allocate and organise output

        root = [];

        if isfield(source,'parclabel')
            for ii = 1:numel(source.parclabel)
                root.(regexprep(lower(source.parclabel{ii}), ' +', '_')) = [];
            end
        else
            source.parcsub = unique(sort(Ni(Ni>0),'ascend')); %Loop through all

            for ii = 1:numel(labs)
                root.(fullfile(['roi_',num2str(source.parcsub(ii))])) = [];
            end
        end

        if numel(fields(root))~=numel(source.parcsub)
            disp('ERROR - Uneven parcellation index vs. labels input');
            disp('Stopping');
            break
        end


        % Now loop over all subjects (tracts) and parcellation ROIs to
        % calculate
        fx=fields(root);

        for ii = 1:numel(tractimg)
            T = nifti(tractimg{ii});
            Ti = T.dat(:,:,:);
            if strcmp(cstype,'adcon') %Adaptive Ratio of Connectivity
                thr=0;
                for k = 1:numel(source.parcsub)

                    mx = max(Ti(Ni == source.parcsub(k)));
                    if mx>thr
                        thr = mx;
                    end
                end
                thr = (thr*csthr);
                for k = 1:numel(source.parcsub)
                    ol = (Ti(Ni == source.parcsub(k)));
                    vol = sum((Ni(:) == source.parcsub(k)));
                    tot = sum(ol>thr);
                    root.(fx{k})(:,ii) = tot/vol;
                end
                root.threshold(:,ii) = thr;
            end

            if strcmp(cstype,'streamdensity') %Stream Density
                thr = csthr/100; %0 and 1 range
                for k = 1:numel(source.parcsub)

                    ol = (Ti(Ni == source.parcsub(k)));
                    vol = sum((Ni(:) == source.parcsub(k)));
                    tot = sum(ol>thr);
                    root.(fx{k})(:,ii) = tot/vol;
                end
                root.threshold(:,ii) = options.constrength.globalthr;
            end
        end

        filename = fullfile(outdir,['connectivity-strength_seed-',seedname,'_desc-parc_',(regexprep(lower(source.parcname), ' +', '-')),'.tsv']);
        spm_save(filename{1},root);
        results=tractbox_summarystats(root);
        filename = fullfile(outdir,['connectivity-strength_seed-',seedname,'_desc-parc_',(regexprep(lower(source.parcname), ' +', '-')),'_summary-stats.tsv']);
        spm_save(filename{1},results);

        % Parcellation details:
        root=[];
        root.index = source.parcsub;
        root.name = source.parclabel;
        filename = fullfile(outdir,['connectivity-strength_seed-',seedname,'_desc-parc-',(regexprep(lower(source.parcname), ' +', '-')),'_dseg.tsv']);
        spm_save(filename{1},root);
        jsonname = fullfile(outdir,['connectivity-strength_seed-',seedname,'_desc-parc-',(regexprep(lower(source.parcname), ' +', '-')),'.json']);
        if iscell(jsonname)
            jsonname=cell2mat(jsonname);
        end
        rootjson = tractbox_connectivitystrength_json(cstype,csthr,source,filename{1});
        spm_jsonwrite(jsonname,rootjson,struct('indent','  '));
    end
else
    disp('Not yet implemented!');
    return
end
end
