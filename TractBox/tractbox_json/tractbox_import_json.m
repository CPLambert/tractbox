function root=tractbox_import_json
%% function root=tractbox_import_json
% Structured array for tractbox imported data - paired with JSON for BIDS
%_______________________________________________________________________
% Version History:
% Version 1.0, December 2021
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

root.subject = 'Subject Identifier (string)';
root.seed = 'Tractography Seed Identifier (string)';
root.output = 'Output path';
root.dot = 'Path to FSL .dot file (omatrix2)';
root.seedspace = 'Path to example nifti in seedspace';
root.dim = 'Original seed image dimensions';
root.mat = 'Original seed image mat from nifti header';
root.mat0 = 'Original seed image mat0 from nifti header';
root.seedcoord = 'Seed coordinates (x,y,z)';
root.tractcoord = 'Tract coordinates (x,y,z)';
end