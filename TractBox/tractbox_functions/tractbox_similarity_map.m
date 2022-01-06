function tractbox_similarity_map(dotimport,region,label,options)
%% function tractbox_similarity_map(dotimport,region,label,options)
% Import tractography matrix and then compare to a canonical target
% or tracts to see how similar each seed is. At the moment will output both
% absolute correlation and euclidean distance. Not sure what is best. Tried
% Earth Movers Distance but results did not look as good so removed.
%_______________________________________________________________________
% Version History:
% Version 1.0, January 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

% Import the tract data
load(dotimport,'root');
XO=tractbox_import_dot(root.dot,0);

for k = 1:numel(region)
    N = nifti(region{k});Ni=N.dat(:,:,:);
    target=Ni(root.tractcoord);

    if options.heatmap.correlation
        Cm = zeros(root.dim);
    end

    if options.heatmap.ed
        Em = zeros(root.dim);
    end

    for i = 1:numel(root.seedcoord)

        X=target';
        Y=XO(:,i)';
        Y=Y./(max(Y));

        if options.heatmap.correlation
            Cm(root.seedcoord(i)) = abs(corr(X',Y'));
        end

        if options.heatmap.ed
            Em(root.seedcoord(i)) = pdist2(X,Y,'euclidean');
        end
    end

    if options.heatmap.correlation
        filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_target-',label{k},'_heatmap-correlation.nii']);
        N.dat.fname=filename;
        N.dat(:,:,:) = Cm;
        create(N);
    end

    if options.heatmap.ed
        filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_target-',label{k},'_heatmap-euclidean.nii']);
        N.dat.fname=filename;
        N.dat(:,:,:) = Em;
        create(N);
    end
end
end


