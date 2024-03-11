function tractbox_concomp(job)

%% Super quick and basic.. nothing clever here! Needs more time... sometime...
options=tractbox_defaults;
disp(options.tractbox.version);
disp('Running connectivity strength statistics');
if isfield(job.concompin{1},'twosample')
    for k = 1:numel(job.concompin{1}.twosample)
        if isfield(job.concompin{1}.twosample{k},'csmat')
            imp = job.concompin{1}.twosample{k}.csmat;
        end
        if isfield(job.concompin{1}.twosample{k},'cstest')
            tec = job.concompin{1}.twosample{k}.cstest;
        end
        if isfield(job.concompin{1}.twosample{k},'testname')
            testname = job.concompin{1}.twosample{k}.testname;
        end
    end

        root=[];
        %----- new code segment: start -------
        root1 = spm_load(imp{1});root2 = spm_load(imp{2});

        if isfield (root1,'threshold')
        root1=rmfield(root1,'threshold');
        end

        if isfield (root2,'threshold')
        root2=rmfield(root2,'threshold');
        end
        %----- new code segment: end -------
        fx=fields(root1);

        pval = zeros(numel(fx),1);

        for i = 1:numel(fx)
            if strcmp(tec,'ranksum')
                [pval(i),~]=ranksum(root1.(fx{i}),root2.(fx{i}));
            end

            if strcmp(tec,'ttest2')
                [~,pval(i)]=ttest2(root1.(fx{i}),root2.(fx{i}));
            end
        end

        fdr = pval_adjust(pval,'fdr'); %External code

        for i = 1:numel(fx)
            root.(fx{i}) = [pval(i) fdr(i)];
        end

        [aa,~,~] = fileparts(imp{1});
        newfile = fullfile(aa,[testname,'.tsv']);
        spm_save(newfile,root);

    elseif isfield(job.concompin{1},'correlate')
        imp = job.concompin{1}.correlate{1}.cscorin;
        x = job.concompin{1}.correlate{2}.metric;
        testname = job.concompin{1}.correlate{3}.metricname;
        testname = (regexprep(lower(testname), ' +', '-'));

        root=[];
        root1 = spm_load(imp{1});rmfield(root1,'threshold');
        fx=fields(root1);

        pval = zeros(numel(fx),1);
        rho = pval;
        for i = 1:numel(fx)
            [rho(i),pval(i)]=corr(root1.(fx{i})(:),x(:),'type','Spearman');
        end

        fdr = pval_adjust(pval,'fdr'); %External code

        for i = 1:numel(fx)
            root.(fx{i}) = [rho(i) pval(i) fdr(i)];
        end

        [aa,bb,~] = fileparts(imp{1});
        newfile = fullfile(aa,[bb,'_corr-',testname,'.tsv'])
        spm_save(newfile,root);

    end
end