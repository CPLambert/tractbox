function options = tractbox_defaults
%% Define the default options for the SPM toolbox: TractBox
%_______________________________________________________________________
% Version History:
% Version 1.1, December 2022 - Update with tseg added
% Version 1.0, December 2021
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------
options.tractbox.version = 'Tractbox Version 1.2, 7th December 2022';
options.fsl.path='/usr/local/fsl/bin'; %Easiest to define manually
[options.tractbox.path,~,~]=fileparts(which('tbx_cfg_tractbox'));

%% tractbox_import_dot
% Defaults for importing FSL dot file
options.dot.thr=0;%Default threshold
options.dot.step=100;%Default matrix building step

%% Build tracts
options.tract.thr = 0;
options.tract.tracttype ='maxnorm';

%% Featuremap tractogram threshold
options.featuremap.thr = 0; %Recommend to apply no thresholding/binarisation - You induce false features in the data
options.featuremap.edstep = 1;

%% Save individual correlation and ED matricies - Useful if doing other types of clustering. Default off
options.save.correlation = false;
options.save.ed = false;

%% Featuremaps - Whole matrix
options.featuremap.matrix.mean = true;
options.featuremap.matrix.median = true;
options.featuremap.matrix.var = true;
options.featuremap.matrix.kurtosis = true;
options.featuremap.matrix.skew = true;
options.featuremap.matrix.sum = true;

%% Featuremaps - Local gradients
options.featuremap.gradient.mean = true;
options.featuremap.gradient.median = true;
options.featuremap.gradient.var = true;
options.featuremap.gradient.kurtosis = true;
options.featuremap.gradient.skew = true;
options.featuremap.gradient.sum = true;

%% Heatmaps - I don't know which would be best - Could add more
options.heatmap.correlation = true;
options.heatmap.ed = true;

%% Map2ROI - I don't know which would be best - Could add more
options.map2roi.label = 'desc-map2roi-metrics';
options.map2roi.type = true; %Threshold the highest

%% CONNECTIVITY STRENGTH (NORMALISED STREAMLINE DENSITY)
options.constrength.globalthr = 0.001;% Anything below 0.01% probably unreliable (this is the same as PiCO of 5 using default 5k samples)

%% TRACTOGRAPHIC SEGMENTATION - THALAMUS
options.tseg.thalamus.tpm = '/Users/clambert/Desktop/4FF_TSEG_THALAMUS/TPM/rtseg-thalamus-tpm_9nuclei.nii';%fullfile(options.tractbox.path,'tractbox_tpm','tseg-thalamus-tpm_9nuclei.nii');
options.tseg.thalamus.mask = '/Users/clambert/Desktop/4FF_TSEG_THALAMUS/TPM/rtseg-thalamus-mask.nii';%fullfile(options.tractbox.path,'tractbox_tpm','tseg-thalamus-mask.nii');
options.tseg.thalamus.folder = 'spm-tsegment-thalamus';
options.tseg.thalamus.coreg.sep = [2 2];
options.tseg.thalamus.coreg.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
options.tseg.thalamus.coreg.fwhm = [2 2];
options.tseg.thalamus.coreg.interp = 1;

options.tseg.thalamus.seg.biasreg = 10;
options.tseg.thalamus.seg.biasfwhm = inf;
options.tseg.thalamus.seg.write = [0 0];
options.tseg.thalamus.seg.ngaus = 6;
options.tseg.thalamus.seg.native = [1 0];
options.tseg.thalamus.seg.warped = [0 0];
options.tseg.thalamus.seg.mrf = 0;
options.tseg.thalamus.seg.cleanup = 1;
options.tseg.thalamus.seg.reg = [0 0.001 0.5 0.05 0.2];
options.tseg.thalamus.seg.fwhm = 0;
options.tseg.thalamus.seg.samp = 1;
options.tseg.thalamus.seg.write = [0 0];

options.tseg.thalamus.labels{1} = 'vl';
options.tseg.thalamus.labels{2} = 'pul';
options.tseg.thalamus.labels{3} = 'gen';
options.tseg.thalamus.labels{4} = 'vp';
options.tseg.thalamus.labels{5} = 'postant';
options.tseg.thalamus.labels{6} = 'va';
options.tseg.thalamus.labels{7} = 'ant';
options.tseg.thalamus.labels{8} = 'md';
options.tseg.thalamus.labels{9} = 'ilam'; %may include edge voxels

% These are not written out by default
options.tseg.thalamus.labels{10} = 'midline';
options.tseg.thalamus.labels{11} = 'edges';
options.tseg.thalamus.labels{12} = 'nonthal';

options.tseg.thalamus.debug = false;

end

