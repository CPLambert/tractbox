function options = tractbox_defaults
%% Define the default options for the SPM toolbox: TractBox
%_______________________________________________________________________
% Version History:
% Version 1.0, December 2021
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------
options.fsl.path='/usr/local/fsl/bin'; %Easiest to define manually

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
end

