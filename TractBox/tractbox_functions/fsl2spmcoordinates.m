function lkp = fsl2spmcoordinates(xyz,d)
%% function lkp = fsl2spmcoordinates(x,y,z,d)
% Convert FSL coordinates to MATLAB/SPM. Offset caused by differences in
% where counting starts - Note, this is now historic and has been corrected
% in FSL outputs.
%_______________________________________________________________________
% Version History:
% Version 1.1, December 2021 - Offsets removed
% Version 1.0, February 2012
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

lkp = double(xyz(:,1)+d(1)*((xyz(:,2)) + d(2)*(xyz(:,3)))); 
end