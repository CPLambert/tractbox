function cfg = tbx_cfg_tractbox
% Tractography SPM Toolbox (TractBox)
% Collection of useful functions for tractography data in SPM.
%_______________________________________________________________________
% Version History:
% Version 1.2, December 2022
% - TSEG
% Version 1.1, August 2022
% - Add in thalamic segment tools
% Version 1.0, December 2021
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

if ~isdeployed
    addpath(fullfile(spm('Dir'),'toolbox','TractBox'));
end

options=tractbox_defaults;

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

% ---------------------------------------------------------------------
% BINARY ROIs
% ---------------------------------------------------------------------
roiimg         = cfg_files;
roiimg.tag     = 'roiimg';
roiimg.name    = 'Binary ROI image(s)';
roiimg.help    = {'Select some binary ROIs'};
roiimg.filter = 'nii';
roiimg.ufilter = '.*';
roiimg.num     = [1 Inf];

%% Measures options
metric            = cfg_entry;
metric.tag        = 'metric';
metric.name       = 'Add measures';
metric.help       = {'Add vector of measures'};
metric.strtype    = 'e';
metric.val        = {0};


% ---------------------------------------------------------------------
% Map images
% ---------------------------------------------------------------------
mapimg         = cfg_files;
mapimg.tag     = 'mapimg';
mapimg.name    = 'Map image(s)';
mapimg.help    = {'Select map images (either heat map or metric to ROI)'};
mapimg.filter = 'nii';
mapimg.ufilter = '.*';
mapimg.num     = [1 Inf];

%% Map threshold
mapthresh            = cfg_entry;
mapthresh.tag        = 'mapthresh';
mapthresh.name       = 'Add percentage threshold';
mapthresh.help       = {'Percentage threshold between 0 - 1: Defaults at top 10% of values (i.e. 0.9)'};
mapthresh.strtype    = 'e';
mapthresh.val        = {0.9};

threshtype              = cfg_menu;
threshtype.tag          = 'threshtype';
threshtype.name         = 'Threshold type';
threshtype.help         = {'Apply threshold to keep the highest values (default) or lowest values'};
threshtype.labels       = {'High','Low'};
threshtype.values       = {true,false};
threshtype.val           = {options.map2roi.type};

% ---------------------------------------------------------------------
% Select tract images
% ---------------------------------------------------------------------
tractimg         = cfg_files;
tractimg.tag     = 'tractimg';
tractimg.name    = 'Tractography images';
tractimg.help    = {'Select tractography images'};
tractimg.filter = 'nii';
tractimg.ufilter = '.*';
tractimg.num     = [1 Inf];

% ---------------------------------------------------------------------
% Select parcellation images
% ---------------------------------------------------------------------
cstype              = cfg_menu;
cstype.tag          = 'cstype';
cstype.name         = 'Connectivity Strength Method';
cstype.help         = {'There are different ways of calculating connectivity strength:'
    'Streamline Density: Set one global threshold for all individuals, then count the streamlines arriving at an ROI, as a ratio of ROI volume'
    'Normalised Streamline Density: As above but as a percentage of maximum connectivity strength to all ROIs of interest for each individual (i.e. attempting to factor head size confound)'};
cstype.labels       = {'Streamline Density'
    'Normalised Streamline Density'};
cstype.values       = {'streamdensity'
    'adcon'};
cstype.val           = {'adcon'};

% ---------------------------------------------------------------------
% Select parcellation images
% ---------------------------------------------------------------------
csthr            = cfg_entry;
csthr.tag        = 'csthr';
csthr.name       = 'Connectivity strength threshold';
csthr.help       = {'Many methods require a thresholding step. Specify here (as a percentage of the maximum PiCO of interest) or set to zero to disable. It has been set to 0.02 used for the Adaptive Ratio of Connectivity'};
csthr.strtype    = 'e';
csthr.val        = {0.02};

% ---------------------------------------------------------------------
% Select parcellation images
% ---------------------------------------------------------------------
parcimg         = cfg_files;
parcimg.tag     = 'parcimg';
parcimg.name    = 'Parcellation image';
parcimg.help    = {'Parcellation (or ROI) image to loop through and calculate connectivity strength'
    ''
    'NOTE: If only one image is selected, assumes this is all happening in group average space and will loop all input tracts through this,'
    'if more than one input, then it has to equal the number of tracts (i.e. individual subject space with unique roi per person.'
    'If you "repeat" these menu options, you can loop over different parcellation options (e.g. left, right, subcortical etc.,), using the'
    'parcellation subselect options below (otherwise will simply loop over all labels 1..N in the parcellation)'};
parcimg.filter = 'nii';
parcimg.ufilter = '.*';
parcimg.num     = [1 Inf];

% ---------------------------------------------------------------------
% Subselect parcellation labels
% ---------------------------------------------------------------------
parcsub            = cfg_entry;
parcsub.tag        = 'parcsub';
parcsub.name       = 'Sub-select parcellation';
parcsub.help       = {'Add vector of parcellation values to sub-select'};
parcsub.strtype    = 'e';
parcsub.val        = {[]};

% ---------------------------------------------------------------------
% Parcellation labels
% ---------------------------------------------------------------------
parclabel                = cfg_entry;
parclabel.tag            = 'parclabel';
parclabel.name           = 'Parcellation labels';
parclabel.help           = {'Labels for parcellations provided in ascending value order (if not sub-selecting), or 1:1 mapping with sub-select range'};
parclabel.strtype        = 's+';
parclabel.num            = [1 Inf];

% ---------------------------------------------------------------------
% Parcellation name
% ---------------------------------------------------------------------
parcname                = cfg_entry;
parcname.tag            = 'parcname';
parcname.name           = 'Parcellation name';
parcname.help           = {'Parcellation name for output'};
parcname.strtype        = 's';
parcname.num            = [0 Inf];

% ---------------------------------------------------------------------
% Parcellation name
% ---------------------------------------------------------------------
parcdesc               = cfg_entry;
parcdesc.tag            = 'parcdesc';
parcdesc.name           = 'Parcellation description';
parcdesc.help           = {'Parcellation description for .json'};
parcdesc.strtype        = 's';
parcdesc.num            = [0 Inf];

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

% ---------------------------------------------------------------------
% 2. FEATUREMAPS - Note fix the split version
% ---------------------------------------------------------------------


featureinput              = cfg_repeat;
featureinput.tag          = 'featureinput';
featureinput.name         = 'Feature Maps Input';
featureinput.values       = {tractboxfile edstepsize savecoor saveed threshold};
featureinput.val          = {tractboxfile};
featureinput.help         = {'Create Feature Maps: Defaults: Output will be input, default threshold 0, Single pass ED calculation, Sobel-Gradient 26 Connected ED Neighbourhood, ED and Correlation matrix not saved'};

% ---------------------------------------------------------------------
% 3. HEATMAPS
% ---------------------------------------------------------------------
heatinput              = cfg_repeat;
heatinput.tag          = 'heatinput';
heatinput.name         = 'Heatmap Input';
heatinput.values       = {regiondata};
heatinput.val          = {regiondata};
heatinput.help         = {'Create HeatMaps (correlation and euclidean distance - change defaults to modify)'};

% ---------------------------------------------------------------------
% 3. METRIC2ROI
% ---------------------------------------------------------------------
met2roiinput              = cfg_branch;
met2roiinput.tag          = 'met2roiinput';
met2roiinput.name         = 'Metric2ROI Input';
met2roiinput.val          = {metric roiimg targetname};
met2roiinput.help         = {'Map some measures (e.g. clinical score) to a binary ROI. For example, this can be used to create individual specific ROIS associated with clinical score, that when averaged map the region associated with top/bottom scores'};

% ---------------------------------------------------------------------
% 4. MAP2ROI
% ---------------------------------------------------------------------
map2roiinput              = cfg_branch;
map2roiinput.tag          = 'map2roiinput';
map2roiinput.name         = 'Metric2ROI Input';
map2roiinput.val          = {targetname mapthresh threshtype roiimg mapimg};
map2roiinput.help         = {'Compare a measurement map to a region of interest (e.g., distance between heatman and DBS ROI). Can do many ROIs and Maps, but will apply same target names, thresholds etc accross group'};

% ---------------------------------------------------------------------
% 3. CONNECTIVITY STRENGTH MAPS
% ---------------------------------------------------------------------
parcin              = cfg_repeat;
parcin.tag          = 'parcin';
parcin.name         = 'Parcellation input';
parcin.values       = {parcimg parcsub parclabel parcname parcdesc};
parcin.val          = {parcimg parcsub parclabel parcname parcdesc};
parcin.help         = {'One output will be generated per parcellation subselection. Repeat if >1 required'};

csin              = cfg_repeat;
csin.tag          = 'csin';
csin.name         = 'Connectivity strength inputs';
csin.values       = {seedname outdir tractimg cstype csthr parcin tractboxfile};
csin.val          = {seedname outdir tractimg cstype csthr parcin};
csin.help         = {'Import new FSL .dot'};

% ---------------------------------------------------------------------
% 3. CONNECTIVITY STRENGTH COMPARE
% ---------------------------------------------------------------------
csmat         = cfg_files;
csmat.tag     = 'csmat';
csmat.name    = 'Compare two matrices';
csmat.help    = {'TractBox connectivity strength files'};
csmat.filter = 'tsv';
csmat.ufilter = '.*';
csmat.num     = [0 2];

cstest              = cfg_menu;
cstest.tag          = 'cstest';
cstest.name         = 'Statistical test';
cstest.help         = {'Only the absolute basics have been implemented. Select one:'};
cstest.labels       = {'Two sample T-test'
    'Wilcoxon rank sum'};
cstest.values       = {'ttest2'
    'ranksum'};
cstest.val           = {'ranksum'};

% ---------------------------------------------------------------------
% Parcellation name
% ---------------------------------------------------------------------
testname               = cfg_entry;
testname.tag            = 'testname';
testname.name           = 'Name of contrast';
testname.help           = {'Name of contrast'};
testname.strtype        = 's';
testname.val        = {};

twosample              = cfg_repeat;
twosample.tag          = 'twosample';
twosample.name         = 'Two-sample comparison';
twosample.values       = {csmat cstest testname};
twosample.val          = {csmat cstest testname};
twosample.help         = {'Do some simple statistics on basic connectivity matrices'};

cscorin         = cfg_files;
cscorin.tag     = 'cscorin';
cscorin.name    = 'Input matrix';
cscorin.help    = {'TractBox connectivity strength file (1)'};
cscorin.filter = 'tsv';
cscorin.ufilter = '.*';
cscorin.num     = [0 1];

metricname            = cfg_entry;
metricname.tag        = 'metricname';
metricname.name       = 'Metric name';
metricname.help       = {'For file labelling'};
metricname.strtype    = 's';
metricname.val        = {};

correlate              = cfg_repeat;
correlate.tag          = 'correlate';
correlate.name         = 'Correlate';
correlate.values       = {cscorin metric metricname};
correlate.val          = {cscorin metric metricname};
correlate.help         = {'Row-wise Spearmans correlation vs some measure of interest'};

concompin              = cfg_repeat;
concompin.tag          = 'concompin';
concompin.name         = 'Connectivity strength: Statisics';
concompin.values       = {twosample correlate};
concompin.val          = {};
concompin.help         = {'Do some stats on basic connectivity matrix'};

%% TSEG options:
% ---------------------------------------------------------------------
% Select LEFT ROI data
% ---------------------------------------------------------------------
leftdir              = cfg_files;
leftdir.tag          = 'leftdir';
leftdir.name         = 'Featuremap directory: LEFT';
leftdir.help         = {'Select directory containing featuremaps for the left side. Note, this assumes you have run create feature maps and the directory contains all feature maps for the left side. If selecting multiple subjects, make sure it is the same ordering as the right'};
leftdir.filter       = 'dir';
leftdir.ufilter      = '.*';
leftdir.num          = [0 inf];

% ---------------------------------------------------------------------
% Select RIGHT ROI data
% ---------------------------------------------------------------------
rightdir              = cfg_files;
rightdir.tag          = 'rightdir';
rightdir.name         = 'Featuremap directory: RIGHT';
rightdir.help         = {'Select directory containing featuremaps for the left side. Note, this assumes you have run create feature maps and the directory contains all feature maps for the left side. If selecting multiple subjects, make sure it is the same ordering as the right'};
rightdir.filter       = 'dir';
rightdir.ufilter      = '.*';
rightdir.num          = [0 inf];

jointdir              = cfg_files;
jointdir.tag          = 'jointdir';
jointdir.name         = 'Featuremap directory: JOINT NORMALISED DATA';
jointdir.help         = {'Select directory containing the JOINT NORMALISED featuremaps if already calculated'};
jointdir.filter       = 'dir';
jointdir.ufilter      = '.*';
jointdir.num          = [0 inf];

% % ---------------------------------------------------------------------
% % Enter BIDS derivative, if variable is derived from others?
% % ---------------------------------------------------------------------
% makegrad          = cfg_menu;
% makegrad.tag      = 'makegrad';
% makegrad.name     = 'Calculate Gradients?';
% makegrad.help     = {'Once segmentation is complete, calculate gradients within each nuclei - Note this involves reading in the whole tractogram and can be very memory intensive'};
% makegrad.labels   = {'Yes','No'};
% makegrad.values   = {true,false};
% makegrad.val      = {true};

% ---------------------------------------------------------------------
% Enter BIDS derivative, if variable is derived from others?
% ---------------------------------------------------------------------
tsegver          = cfg_menu;
tsegver.tag      = 'tsegver';
tsegver.name     = 'Tseg version?';
tsegver.help     = {'Which structure?'};
tsegver.labels   = {'thalamus'};
tsegver.values   = {'thalamus'};
tsegver.val      = {'thalamus'};

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
% 2. MAT2TRACT - Make tract image from imported matrix
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
% 4. TRACTOGRAPHIC SEGMENTATION - Currently just thalamus
% ---------------------------------------------------------------------
fjoin     = cfg_exbranch;
fjoin.tag = 'fjoin';
fjoin.name= 'Combine feature maps';
fjoin.val = {leftdir,rightdir};
fjoin.help= {'Calculate combined normalised Feature Maps from Tractography data'};
fjoin.prog= @tractbox_fjoin;

% This is the version where you feed left and right directories seperately. Removed for now. 
% tseg     = cfg_exbranch;
% tseg.tag = 'tseg';
% tseg.name= 'Thalamus';
% tseg.val = {tsegver,leftdir,rightdir};
% tseg.help= {'Use pre-calculated combined normalised Feature Maps from Tractography data to run segmentation as described in Lambert et al., 2017 Defining thalamic nuclei and topographic connectivity gradients in vivo. Currently just thalamus implemented.'};
% tseg.prog= @tractbox_tseg_thalamus;

% ---------------------------------------------------------------------
% Select parcellation images
% ---------------------------------------------------------------------
thwarp         = cfg_files;
thwarp.tag     = 'thwarp';
thwarp.name    = 'Thalamic deformation field';
thwarp.help    = {'Deformation fields for thalamus'};
thwarp.filter = 'nii';
thwarp.ufilter = '.*';
thwarp.num     = [1 Inf];
% ---------------------------------------------------------------------
% 4. TRACTOGRAPHIC SEGMENTATION - Currently just thalamus
% ---------------------------------------------------------------------
thseg     = cfg_exbranch;
thseg.tag = 'thseg';
thseg.name= 'Thalamus Segment';
thseg.val = {jointdir,thwarp};
thseg.help= {'Use pre-calculated combined normalised Feature Maps from Tractography data to run segmentation as described in Lambert et al., 2017 Defining thalamic nuclei and topographic connectivity gradients in vivo. Currently just thalamus implemented.'};
thseg.prog= @tractbox_tseg2_thalamus;

% ---------------------------------------------------------------------
% 5. CREATE FEATURE MAPS - Featuremaps from tracts
% ---------------------------------------------------------------------
heatmaps     = cfg_exbranch;
heatmaps.tag = 'heatmaps';
heatmaps.name= 'Create Heat Maps';
heatmaps.val = {tractboxfile heatinput};
heatmaps.help= {'Create voxel-wise heat maps from tractography matrix and target tract'};
heatmaps.prog= @tractbox_heatmaps;

% ---------------------------------------------------------------------
% 6. MEASURE2ROI MAPS
% ---------------------------------------------------------------------
met2roi     = cfg_exbranch;
met2roi.tag = 'met2roi';
met2roi.name= 'Map Measures to ROI';
met2roi.val = {met2roiinput};
met2roi.help= {'Map a phenotyping measure to a binary ROI - Lots of uses e.g. mapping symptom severity to brain regions/lesions'};
met2roi.prog= @tractbox_met2roi;

% ---------------------------------------------------------------------
% 7. COMPARE MAPS
% ---------------------------------------------------------------------
map2roi     = cfg_exbranch;
map2roi.tag = 'map2roi';
map2roi.name= 'Map vs ROI analysis';
map2roi.val = {map2roiinput};
map2roi.help= {'This will compare a measurement map to a region of interest'
    ''
    'It is primarily designed to compare heat maps to a given volume of tissue activation, but the same principle could also be used with the map output from Measures2ROI'
    ''
    'It will calculate the centroid location, taking a set threshold for the input map (generally top 10%), and measure the distance between the two'
    ''
    'It will also extract some measures from within the centroid from the input map'
    ''
    'These will all be written out to a .tsv paired with a .json, generally in the input folder'};
map2roi.prog= @tractbox_map2roi;

% ---------------------------------------------------------------------
% 8. CONNECTIVITY STRENGTH: CALCULATE
% ---------------------------------------------------------------------
constrength     = cfg_exbranch;
constrength.tag = 'constrength';
constrength.name= 'Connectivity Strength: Calculate';
constrength.val = {csin};
constrength.help= {'This will calculate some metric of "connectivity strength" for tractography results and some ROI(s)'
    ''
    'Note your tractography images should be normalised to fall between 0 and 1 by dividing by the number of streamlines you used, otherwise you may get odd results with some settings' 
    ''
    'This whole brain (i.e. pattern of connectivity) or between limited ROIs'
    ''
    'You can either used previously calculated tract images (e.g., maximum projection images) or you can do this voxelwise for each individual via the tractbox import file (obviously much slower)'
    ''
    'Unless specified, it will simply loop through all the binary ROIs in an image'
    ''
    'If you specify the ROIs, it will output a matrix just for those. If you pair this with labels, then it will assign these as the variables'
    ''
    'Currently implemented connectivity strength metrics: Normalised strength of connectivity (see Lambert et al., 2013 for details'};
constrength.prog= @tractbox_constrength;

% ---------------------------------------------------------------------
% 9. CONNECTIVITY STRENGTH: COMPARE
% ---------------------------------------------------------------------
concomp     = cfg_exbranch;
concomp.tag = 'concomp';
concomp.name= 'Connectivity Strength: Statistics';
concomp.val = {concompin};
concomp.help= {'Very simple implementation - Take two connectivity matrices, and compare ROI for ROI accross population'
    'Will correct for multiple comparisons. Assumes the two matrices are same size (they should be!)'};
concomp.prog= @tractbox_concomp;


%% MAIN SPM MENU: TRACTBOX
cfg                 = cfg_choice;
cfg.tag             = 'tractbox';
cfg.name            = 'TractBox';
cfg.values          = {importfsldot maketract constrength concomp featuremaps fjoin thseg heatmaps met2roi map2roi};
end