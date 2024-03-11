function tractbox_create_featureimage(dotimport,options)
%% function tractbox_create_featureimage(root)
% This creates the feature images as described in:
%
% "Lambert C, Simon H, Colman J, Barrick TR. Defining thalamic nuclei and
% topographic connectivity gradients in vivo. Neuroimage. 2017
% Sep;158:466-479"
%
% Please cite this paper if used.
%
% In brief, this function uses a structured array storing all the necessary
% lookup data for a tractogram. It will calculate the correlation and
% euclidean distance matrices, then generate feature images (statistical
% moments) from the whole ED matric, and the first derivative (calculated
% here as the sobel edge gradient between the 27 connected neighbourhood
% for a voxel. The output will be _featuremap_[moment] where [moment] may
% be mean, variation, kurtosis, skew or median. Admitted these are a lot of
% images (10 per person, five for whole ED matric and five for connected gradients),
% but it is currently unclear how useful each of these are. If you wish to
% modify this to be more selective in the future, change the defaults file
%_______________________________________________________________________
% Version History:
% Version 1.0, January 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

if nargin <2
    options = tractbox_defaults;
end

load(dotimport,'root');
ksize=3;
XO=tractbox_import_dot(root.dot,options.featuremap.thr);
BCi=tractbox_fast_corr(XO,XO); clear XO;

if options.save.correlation
    disp('Saving correlation matrix');
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_tractbox-correlation']);
    save(filename,'BCi');
end

ED = squareform(tractbox_distance_ed(BCi,options.featuremap.edstep)); clear BCi;

if options.save.ed
    disp('Saving Euclidean Distance matrix');
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_tractbox-ed']);
    save(filename,'ED');
end

list = tractbox_find_neighbours(root.seedcoord,root.dim);
xN = tractbox_sobel_matrix(ED,list);

%% Let's make some images:
N = nifti(root.seedspace);

if options.featuremap.matrix.mean
    val = mean(ED,1);
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-mean.nii']);
    N.dat.fname=filename;
    N.dat.dtype='FLOAT32';
    N.dat(:,:,:) = Z;
    create(N);
    type = 'Euclidean Distance Matrix: Voxel-wise mean';
    sourcesize = size(val,1);
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
        filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-mean.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end

if options.featuremap.matrix.median
    val = median(ED,1);
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-median.nii']);
    N.dat.fname=filename;
    N.dat.dtype='FLOAT32';
    N.dat(:,:,:) = Z;
    create(N);
    type = 'Euclidean Distance Matrix: Voxel-wise median';
    sourcesize = size(val,1);
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
        filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-median.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end

if options.featuremap.matrix.var
    val = var(ED,1);
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-variance.nii']);
    N.dat.fname=filename;
    N.dat.dtype='FLOAT32';
    N.dat(:,:,:) = Z;
    create(N);
    type = 'Euclidean Distance Matrix: Voxel-wise variance';
    sourcesize = size(val,1);
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
        filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-variance.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end

if options.featuremap.matrix.kurtosis
    val = kurtosis(ED,1);
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-kurtosis.nii']);
    N.dat.fname=filename;
    N.dat.dtype='FLOAT32';
    N.dat(:,:,:) = Z;
    create(N);

    type = 'Euclidean Distance Matrix: Voxel-wise kurtosis';
    sourcesize = size(val,1);
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
        filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-kurtosis.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end

if options.featuremap.matrix.skew
    val = skewness(ED,1);
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-skewness.nii']);
    N.dat.fname=filename;
    N.dat.dtype='FLOAT32';
    N.dat(:,:,:) = Z;
    create(N);
    type = 'Euclidean Distance Matrix: Voxel-wise skewness';
    sourcesize = size(val,1);
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
        filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-skewness.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end

if options.featuremap.matrix.sum
    val = sum(ED,1);
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-total.nii']);
    N.dat.fname=filename;
    N.dat.dtype='FLOAT32';
    N.dat(:,:,:) = Z;

    create(N);
    type = 'Euclidean Distance Matrix: Voxel-wise total (summation of vector)';
    sourcesize = size(val,1);
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
        filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-global-total.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end


%% Featuremaps - Local gradients
if options.featuremap.gradient.mean
    val = zeros(size(xN));
    for ii = 1:numel(xN)
        val(ii)=mean(xN{ii}(~isnan(xN{ii})));
    end
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-mean.nii']);
    N.dat.fname=filename;
    N.dat(:,:,:) = Z;
    N.dat.dtype='FLOAT32';
    create(N);

    type = 'Regional Sobel Gradient of Euclidean Distance: Mean of 3D connected neighbourhood';
    sourcesize = ksize^3;
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-mean.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end

if options.featuremap.gradient.median
    val = zeros(size(xN));
    for ii = 1:numel(xN)
        val(ii)=median(xN{ii}(~isnan(xN{ii})));
    end
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-median.nii']);
    N.dat.fname=filename;
    N.dat.dtype='FLOAT32';
    N.dat(:,:,:) = Z;
    create(N);

     type = 'Regional Sobel Gradient of Euclidean Distance: Median of 3D connected neighbourhood';
    sourcesize = ksize^3;
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-median.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end

if options.featuremap.gradient.var
    val = zeros(size(xN));
    for ii = 1:numel(xN)
        val(ii)=var(xN{ii}(~isnan(xN{ii})));
    end
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-variance.nii']);
    N.dat.fname=filename;
    N.dat.dtype='FLOAT32';
    N.dat(:,:,:) = Z;
    create(N);

     type = 'Regional Sobel Gradient of Euclidean Distance: Variance of 3D connected neighbourhood';
    sourcesize = ksize^3;
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-variance.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end

if options.featuremap.gradient.kurtosis
    val = zeros(size(xN));
    for ii = 1:numel(xN)
        val(ii)=kurtosis(xN{ii}(~isnan(xN{ii})));
    end
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-kurtosis.nii']);
    N.dat.fname=filename;
    N.dat.dtype='FLOAT32';
    N.dat(:,:,:) = Z;
    create(N);

     type = 'Regional Sobel Gradient of Euclidean Distance: Kurtosis of 3D connected neighbourhood';
    sourcesize = ksize^3;
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-kurtosis.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end

if options.featuremap.gradient.skew
    val = zeros(size(xN));
    for ii = 1:numel(xN)
        val(ii)=skewness(xN{ii}(~isnan(xN{ii})));
    end
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-skewness.nii']);
    N.dat.fname=filename;
    N.dat.dtype='FLOAT32';
    N.dat(:,:,:) = Z;
    create(N);

     type = 'Regional Sobel Gradient of Euclidean Distance: Skewness of 3D connected neighbourhood';
    sourcesize = ksize^3;
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-skewness.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end

if options.featuremap.gradient.sum
    val = zeros(size(xN));
    for ii = 1:numel(xN)
        val(ii)=sum(xN{ii}(~isnan(xN{ii})));
    end
    Z=zeros(root.dim);
    Z(root.seedcoord)=val;
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-total.nii']);
    N.dat.fname=filename;
    N.dat.dtype='FLOAT32';
    N.dat(:,:,:) = Z;
    create(N);

     type = 'Regional Sobel Gradient of Euclidean Distance: Total summation of 3D connected neighbourhood';
    sourcesize = ksize^3;
    rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport);
    filename = fullfile(root.output,['sub-',root.subject,'_seed-',root.seed,'_threshold-',num2str(options.featuremap.thr),'_featuremap-gradient-total.json']);
    spm_jsonwrite(filename,rootjson,struct('indent','  '));
end
end
