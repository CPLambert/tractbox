function tractbox_mat2tract(job)
%% function tractbox_mat2tract(job)
% A relatively flexible function for extracting and creating tracts from
% tractbox import .mat data
%_______________________________________________________________________
% Version History:
% Version 1.0, December 2021
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

options = tractbox_defaults;
mask=cell(1,2);

for i=1:numel(job.maketracttype)
    if isfield(job.maketracttype(i),'newtractfew')
        tracttype = options.tract.tracttype;
        thr = options.tract.thr;

        for ii = 1:numel(job.maketracttype(i).newtractfew)
            if isfield(job.maketracttype(i).newtractfew{ii},'tractboxfile1')
                load(job.maketracttype(i).newtractfew{ii}.tractboxfile1{1});
                output = root.output;
            end

            %% Sort out the options
            if isfield(job.maketracttype(i).newtractfew{ii},'threshold')
                thr = job.maketracttype(i).newtractfew{ii}.threshold;
            end

            if isfield(job.maketracttype(i).newtractfew{ii},'tracttype')
                tracttype = job.maketracttype(i).newtractfew{ii}.tracttype;
            end

            if isfield(job.maketracttype(i).newtractfew{ii},'maskimg')
                for iii = 1:numel(job.maketracttype(i).newtractfew{ii}.maskimg)
                    [~,mask{iii,1},~]=fileparts(job.maketracttype(i).newtractfew{ii}.maskimg{iii});
                    N=nifti(job.maketracttype(i).newtractfew{ii}.maskimg{iii});
                    mask{iii,2}=find(N.dat(:,:,:)>0);
                end
            end

            if isfield(job.maketracttype(i).newtractfew{ii},'outdir')
                output = job.maketracttype(i).newtractfew{ii}.outdir{1};
            end
        end
        if iscell(root.dot)
            dfile=root.dot{1};
        else
            dfile = root.dot;
        end
        % Get the data
        x0=tractbox_import_dot(dfile,thr);

        %% Now create the image
        Z=zeros(root.dim(1),root.dim(2),root.dim(3));A=Z;

        for k=1:size(mask,1)
            if isempty(mask{k,2})
                clist=root.seedcoord;
            else
                clist=(mask{k,2});
            end

            if strcmp(tracttype,'maxnorm') || strcmp(tracttype,'max')
                for iii = 1:numel(root.seedcoord)
                    if logical(sum(clist==root.seedcoord(iii)))
                        tmp=zeros(size(A));
                        tmp(root.tractcoord)=x0(:,iii);
                        tmp(tmp<thr)=0;
                        Z(tmp>A)=tmp(tmp>A);
                        A(tmp>A)=tmp(tmp>A);
                    end
                end

                if strcmp(tracttype,'maxnorm')
                    Z=Z./(max(Z(:)));
                end
            end

            if strcmp(tracttype,'mean') || strcmp(tracttype,'cuml') || strcmp(tracttype,'bin')
                for iii = 1:numel(root.seedcoord)
                    if logical(sum(clist==root.seedcoord(iii)))

                        tmp=zeros(size(A));
                        tmp(root.tractcoord)=x0(:,iii);
                        tmp(tmp<thr)=0;
                        Z=Z + tmp;
                    end
                end

                if strcmp(tracttype,'mean')
                    Z=Z./(numel(root.seedcoord));
                end

                if strcmp(tracttype,'bin')
                    Z(Z>0)=1;
                end
            end

            if isempty(mask{k,1})
                desc=char(strcat('desc-',tracttype,'-thr-',num2str(thr),'-region-all'));
            else
                desc=char(strcat('desc-',tracttype,'-thr-',num2str(thr),'-region-',mask{k,1}));
            end

            filename = fullfile(output,['sub-',root.subject,'_seed-',root.seed,'_',desc,'-tracts.nii']);
            N=nifti(root.seedspace);
            N.dat.fname=filename;
            N.dat.dtype='FLOAT32';
            N.dat(:,:,:)=Z;
            create(N);
        end
    end

    if isfield(job.maketracttype(i),'newtractmany')
        tracttype = options.tract.tracttype;
        thr = options.tract.thr;

        for ii = 1:numel(job.maketracttype(i).newtractmany)
            if isfield(job.maketracttype(i).newtractfew{ii},'tractboxfile')
                list=(job.maketracttype(i).newtractfew{ii}.tractboxfile);
            end

            %% Sort out the options
            if isfield(job.maketracttype(i).newtractfew{ii},'threshold')
                thr = job.maketracttype(i).newtractfew{ii}.threshold;
            end

            if isfield(job.maketracttype(i).newtractfew{ii},'tracttype')
                tracttype = job.maketracttype(i).newtractfew{ii}.tracttype;
            end
        end

        for ii=1:numel(list)
            load(list{ii});

            if iscell(root.dot)
                dfile=root.dot{1};
            else
                dfile = root.dot;
            end
            % Get the data
            x0=tractbox_import_dot(dfile,thr);

            %% Now create the image

            Z=zeros(root.dim(1),root.dim(2),root.dim(3));A=Z;

            for k=1:size(mask,1)
                if isempty(mask{k,2})
                    clist=root.seedcoord;
                else
                    clist=(mask{k,2});
                end

                if strcmp(tracttype,'maxnorm') || strcmp(tracttype,'max')
                    for iii = 1:numel(root.seedcoord)
                        if logical(sum(clist==root.seedcoord(iii)))
                            tmp=zeros(size(A));
                            tmp(root.tractcoord)=x0(:,iii);
                            tmp(tmp<thr)=0;
                            Z(tmp>A)=tmp(tmp>A);
                            A(tmp>A)=tmp(tmp>A);
                        end
                    end

                    if strcmp(tracttype,'maxnorm')
                        Z=Z./(max(Z(:)));
                    end
                end

                if strcmp(tracttype,'mean') || strcmp(tracttype,'cuml') || strcmp(tracttype,'bin')
                    for iii = 1:numel(root.seedcoord)
                        if logical(sum(clist==root.seedcoord(iii)))

                            tmp=zeros(size(A));
                            tmp(root.tractcoord)=x0(:,iii);
                            tmp(tmp<thr)=0;
                            Z=Z + tmp;
                        end
                    end

                    if strcmp(tracttype,'mean')
                        Z=Z./(numel(root.seedcoord));
                    end

                    if strcmp(tracttype,'bin')
                        Z(Z>0)=1;
                    end
                end

                if isempty(mask{k,1})
                    desc=char(strcat('desc-',tracttype,'-thr-',num2str(thr),'-region-all'));
                else
                    desc=char(strcat('desc-',tracttype,'-thr-',num2str(thr),'-region-',mask{k,1}));
                end

                filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_',desc,'-tracts.nii']);
                N=nifti(root.seedspace);
                N.dat.fname=filename;
                N.dat.dtype='FLOAT32';
                N.dat(:,:,:)=Z;
                create(N);
            end
        end
    end
end
end





