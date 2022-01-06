function cfg = tbx_cfg_tractbox
% Tractography SPM Toolbox (TractBox)
% Collection of useful functions for tractography data in SPM.
%_______________________________________________________________________
% Version History:
% Version 1.0, December 2021
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

if ~isdeployed
    addpath(fullfile(spm('Dir'),'toolbox','TractBox'));
end

%% GENERIC
% ---------------------------------------------------------------------
% Enter filename
% ---------------------------------------------------------------------
filename                 = cfg_entry;
filename.tag             = 'filename';
filename.name            = 'Output filename';
filename.help            = {'Define filename'};
filename.strtype         = 's';
filename.val             = {''};
filename.num             = [0 Inf];

% ---------------------------------------------------------------------
% Select output directory
% ---------------------------------------------------------------------
outdir              = cfg_files;
outdir.tag          = 'outdir';
outdir.name         = 'Output directory';
outdir.help         = {'Select output dirctory'};
outdir.filter       = 'dir';
outdir.ufilter      = '.*';
outdir.num          = [0 1];

% ---------------------------------------------------------------------
% Select tract_space_coords_for_fdt_matrix2 files
% ---------------------------------------------------------------------
tractboxfile1         = cfg_files;
tractboxfile1.tag     = 'tractboxfile1';
tractboxfile1.name    = 'TractBox import file';
tractboxfile1.help    = {'Select Tractbox Import'};
tractboxfile1.filter = 'mat';
tractboxfile1.ufilter = '.*';
tractboxfile1.num     = [0 1];

tractboxfile         = cfg_files;
tractboxfile.tag     = 'tractboxfile';
tractboxfile.name    = 'TractBox import files';
tractboxfile.help    = {'TractBox import files'};
tractboxfile.filter = 'mat';
tractboxfile.ufilter = '.*';
tractboxfile.num     = [0 inf];

%% IMPORT OPTIONS
% ---------------------------------------------------------------------
% Subject name
% ---------------------------------------------------------------------
subname                = cfg_entry;
subname.tag            = 'subname';
subname.name           = 'Subject name';
subname.help           = {'Name of the subject'};
subname.strtype        = 's';
subname.num            = [1 Inf];

% ---------------------------------------------------------------------
% Seed name
% ---------------------------------------------------------------------
seedname                = cfg_entry;
seedname.tag            = 'seedname';
seedname.name           = 'Seed name';
seedname.help           = {'Name of the seed'};
seedname.strtype        = 's';
seedname.num            = [1 Inf];

% ---------------------------------------------------------------------
% Select .dot files
% ---------------------------------------------------------------------
fsldot         = cfg_files;
fsldot.tag     = 'fsldot';
fsldot.name    = 'Select FSL .dot file';
fsldot.help    = {'Select FSL .dot file'};
fsldot.filter = 'dot';
fsldot.ufilter = '.*';
fsldot.num     = [1 Inf];

% ---------------------------------------------------------------------
% Select coords_for_fdt_matrix2 files
% ---------------------------------------------------------------------
fslcoord         = cfg_files;
fslcoord.tag     = 'fslcoord';
fslcoord.name    = 'Select FSL coords';
fslcoord.help    = {'Select coords_for_fdt_matrix2: Only required if not in the same folder as input .dot'};
fslcoord.ufilter = '.*';
fslcoord.num     = [0 inf];

% ---------------------------------------------------------------------
% Select tract_space_coords_for_fdt_matrix2 files
% ---------------------------------------------------------------------
fsltract         = cfg_files;
fsltract.tag     = 'fsltract';
fsltract.name    = 'Select FSL tract space coords';
fsltract.help    = {'Select tract_space_coords_for_fdt_matrix2: Only required if not in the same folder as input .dot'};
fsltract.ufilter = '.*';
fsltract.num     = [0 inf];

% ---------------------------------------------------------------------
% Select tract_space_coords_for_fdt_matrix2 files
% ---------------------------------------------------------------------
fslseedspace         = cfg_files;
fslseedspace.tag     = 'fslseedspace';
fslseedspace.name    = 'Seed Space example';
fslseedspace.help    = {'Select example nifti in seed space to base output on'};
fslseedspace.filter = 'nii';
fslseedspace.ufilter = '.*';
fslseedspace.num     = [1 Inf];

%% MAT2TRACT options
threshold            = cfg_entry;
threshold.tag        = 'threshold';
threshold.name       = 'Add PiCO threshold';
threshold.help       = {'Integer value'};
threshold.strtype    = 'i';
threshold.val        = {0};

tracttype              = cfg_menu;
tracttype.tag          = 'tracttype';
tracttype.name         = 'Type of tractography image to make';
tracttype.help         = {'There are several types of tractography image:' 
    'Normalised Maximum Projection: Biggest PICo per voxel over whole seed space, normalised to range 0 - 1 (default)'
    'Maximum Projection: Biggest PICo per voxel over whole seed space'
    'Mean: Voxel-wise average over all seeds per voxel'
    'Cumulative: Voxel-wise summation over all seeds (this is fdt_paths output in FSL)'
    'Binarised'};
tracttype.labels       = {'Normalised Maximum Projection','Maximum Projection','Mean','Cumulative','Binarised'};
tracttype.values       = {'maxnorm','max','mean','cuml','bin'};
tracttype.val           = {'maxnorm'};  

% ---------------------------------------------------------------------
% Subselect seeds from mask
% ---------------------------------------------------------------------
maskimg         = cfg_files;
maskimg.tag     = 'maskimg';
maskimg.name    = 'Mask image';
maskimg.help    = {'Use mask/ROI image to subselect seeds to make tract image'};
maskimg.filter = 'nii';
maskimg.ufilter = '.*';
maskimg.num     = [1 Inf];

% ---------------------------------------------------------------------
% Subselect seeds from file
% ---------------------------------------------------------------------
% maskcoord         = cfg_files;
% maskcoord.tag     = 'maskimg';
% maskcoord.name    = 'Mask coordinate list';
% maskcoord.help    = {'Use list of coordinates from csv/tsv, stored as nx3 (seed: x y z)'};
% maskcoord.filter = 'nii';
% maskcoord.ufilter = '.*';
% maskcoord.num     = [1 Inf];

savecoor              = cfg_menu;
savecoor.tag          = 'savecoor';
savecoor.name         = 'Save Correlation Matrix';
savecoor.help         = {'Save seed x seed correlation matrix (default off)'};
savecoor.labels       = {'Yes','No'};
savecoor.values       = {true,false};
savecoor.val           = {true}; 

saveed              = cfg_menu;
saveed.tag          = 'saveed';
saveed.name         = 'Save Euclidean Distance Matrix';
saveed.help         = {'Save seed x seed Euclidean Distance matrix (default off)'};
saveed.labels       = {'Yes','No'};
saveed.values       = {true,false};
saveed.val           = {true}; 

edstepsize            = cfg_entry;
edstepsize.tag        = 'edstepsize';
edstepsize.name       = 'ED Calculation - Step size';
edstepsize.help       = {'For large seed regions (>1000 voxels), it becomes quicker to calculate ED piecewise (i.e. breaking the matrix into blocks). Generally, I divide the total seed size by ~1500 to set step size'};
edstepsize.strtype    = 'i';
edstepsize.val        = {4};

% ---------------------------------------------------------------------
% Subselect seeds from mask
% ---------------------------------------------------------------------
targetimg         = cfg_files;
targetimg.tag     = 'targetimg';
targetimg.name    = 'Target image(s)';
targetimg.help    = {'Target tract image(s) for comparison against tractogram matrix. Must equal the number of TractBox import files'};
targetimg.filter = 'nii';
targetimg.ufilter = '.*';
targetimg.num     = [1 Inf];

% ---------------------------------------------------------------------
% Subject name
% ---------------------------------------------------------------------
targetname                = cfg_entry;
targetname.tag            = 'targetname';
targetname.name           = 'Target label';
targetname.help           = {'Label for output files, remove whitespace e.g. posterior-stn'};
targetname.strtype        = 's';
targetname.num            = [1 Inf];

%% SUBSECTION OPTIONS
% ---------------------------------------------------------------------
% 1. IMPORT DOT
% ---------------------------------------------------------------------
newdot              = cfg_repeat;
newdot.tag          = 'newdot';
newdot.name         = 'Add new FSL .dot';
newdot.values       = {subname seedname fsldot fslseedspace outdir fslcoord fsltract};
newdot.val          = {subname seedname fsldot fslseedspace};
newdot.help         = {'Import new FSL .dot'};

% ---------------------------------------------------------------------
% 1a. MAT2TRACT (few)
% ---------------------------------------------------------------------
newtractfew              = cfg_repeat;
newtractfew.tag          = 'newtractfew';
newtractfew.name         = 'Create tractography image - Subject wise';
newtractfew.values       = {tractboxfile1 threshold tracttype maskimg outdir};
newtractfew.val          = {tractboxfile1};
newtractfew.help         = {'Create tractography image - Default threshold 0, normalised maximum projection'};

%Feature input:
regiondata              = cfg_branch;
regiondata.tag          = 'regiondata';
regiondata.name         = 'Add new tract target';
% regiondata.values       = {targetimg targetname};
regiondata.val          = {targetimg targetname};
regiondata.help         = {'Add target tractography data for comparison'};

% ---------------------------------------------------------------------
% 1b. MAT2TRACT (many)
% ---------------------------------------------------------------------
newtractmany              = cfg_repeat;
newtractmany.tag          = 'newtractmany';
newtractmany.name         = 'Create tractography image - Group wise';
newtractmany.values       = {tractboxfile threshold tracttype};
newtractmany.val          = {tractboxfile};
newtractmany.help         = {'Create tractography images: Group wise (output will be input, default threshold 0, normalised maximum projection)'};

maketracttype     = cfg_choice;
maketracttype.tag = 'maketracttype';
maketracttype.name= 'MAT2TRACT';
maketracttype.values = {newtractfew newtractmany};
maketracttype.val =  {newtractfew};
maketracttype.help= {'Convert sparse matrix into an actual tractography image. Select few if you want to apply individual level options, select many if you want to apply same settings over many images'};

% Feature maps


featureinput              = cfg_repeat;
featureinput.tag          = 'featureinput';
featureinput.name         = 'Feature Maps Input';
featureinput.values       = {tractboxfile edstepsize savecoor saveed threshold};
featureinput.val          = {tractboxfile};
featureinput.help         = {'Create Feature Maps: Defaults: Output will be input, default threshold 0, Single pass ED calculation, Sobel-Gradient 26 Connected ED Neighbourhood, ED and Correlation matrix not saved'};

heatinput              = cfg_repeat;
heatinput.tag          = 'heatinput';
heatinput.name         = 'Heatmap Input';
heatinput.values       = {regiondata};
heatinput.val          = {regiondata};
heatinput.help         = {'Create HeatMaps (correlation and euclidean distance - change defaults to modify)'};

%% MAIN TOOLBOX OPTIONS
% ---------------------------------------------------------------------
% 1. IMPORT DOT TO TRACTBOX - Most of the functions then handle this array
% ---------------------------------------------------------------------
importfsldot     = cfg_exbranch;
importfsldot.tag = 'importfsldot';
importfsldot.name= 'Import FSL .dot to TractBox';
importfsldot.val = {newdot};
importfsldot.help= {'Import FSL .dot tractography to TractBox array. First four inputs are mandatory. If not specified, output will default to dot input'};
importfsldot.prog= @tractbox_importfsldot;

% ---------------------------------------------------------------------
% 2. IMPORT DOT TO TRACTBOX - Most of the functions then handle this array
% ---------------------------------------------------------------------
maketract     = cfg_exbranch;
maketract.tag = 'maketract';
maketract.name= 'MAT2TRACT';
maketract.val = {maketracttype};
maketract.help= {'Convert sparse matrix into an actual tractography image. Select few if you want to apply individual level options, select many if you want to apply same settings over many images'};
maketract.prog= @tractbox_mat2tract;

% ---------------------------------------------------------------------
% 3. CREATE FEATURE MAPS - Featuremaps from tracts
% ---------------------------------------------------------------------
featuremaps     = cfg_exbranch;
featuremaps.tag = 'featuremaps';
featuremaps.name= 'Create Feature Maps';
featuremaps.val = {featureinput};
featuremaps.help= {'Create Feature Maps from Tractography data as described in Lambert et al., 2017 Defining thalamic nuclei and topographic connectivity gradients in vivo. This moves away from viewing tractography as "tracts" per se but instead ways of sampling brain shape properties from specific spatial locations, and therefore escapes from issues surrounding thresholding, false positives etc.,'};
featuremaps.prog= @tractbox_featuremaps;

% ---------------------------------------------------------------------
% 3. CREATE FEATURE MAPS - Featuremaps from tracts
% ---------------------------------------------------------------------
heatmaps     = cfg_exbranch;
heatmaps.tag = 'heatmaps';
heatmaps.name= 'Create Heat Maps';
heatmaps.val = {tractboxfile heatinput};
heatmaps.help= {'Create voxel-wise heat maps from tractography matrix and target tract'};
heatmaps.prog= @tractbox_heatmaps;

%% MAIN SPM MENU: TRACTBOX
cfg                 = cfg_choice;
cfg.tag             = 'tractbox';
cfg.name            = 'TractBox';
cfg.values          = {importfsldot maketract featuremaps heatmaps};
end