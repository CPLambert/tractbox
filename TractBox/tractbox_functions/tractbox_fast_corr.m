function C = tractbox_fast_corr(X,Y,chk)
%% function C = tractbox_fast_corr(X,Y,chk)
%--------------------------------------------------------------------------
% Faster than MATLAB inbuilt corr, useful for correlating tractogram
% matrices. Note needs to be dimension (seeds x tracts). If too big, see
% tractbox_massive_corr(X,Y). chk flag defaults to true - It will make sure
% the first dimension is the smallest and if not flip the matrix to avoid
% it killing your computer
%--------------------------------------------------------------------------
% C Lambert, J Ashburner - Wellcome Centre for Human Neuroimaging
% Version 1.1 - January 2022, simplified and added to tractbox
% Version 1.0 - May 2012 (originally DTI_fast_corr2)
%--------------------------------------------------------------------------

if nargin < 3,chk = true;end

Sx=size(X);Sy=size(Y);

if chk
    if Sx(2)<Sx(1);X=X';Sx=size(X);end
    if Sy(2)<Sy(1);Y=Y';Sy=size(Y);end
end

if ~isequal(Sy,Sx)
    disp('Error - Matrices need to be the same size');
    return
end

N=Sx(2);
XY   = X*Y'/N;
sX   = sum(X,2)/N;
sX2  = sum(X.^2,2)/N;
sY   = sum(Y,2)/N;
sY2  = sum(Y.^2,2)/N;
C    = (XY - sX*sY')./(sqrt(sX2-sX.^2)*sqrt(sY2-sY.^2)');
end