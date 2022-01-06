function xN = tractbox_sobel_matrix(datamatrix,list)
%% function xN = tractbox_sobel_matrix
% Apply sobel edge detection to 2D matrix generated from 3D image data
% (e.g. correlation, euclidean distance). Returns the gradient. Defaults to
% a 3x3x3 sobel kernel.
%
% In this version, return the 27 connected neighbourhood for each seed.
% Missing data will be NaN.
%_______________________________________________________________________
% Version History:
% Version 1.0, January 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

h = tractbox_sobel_kernel;
ym2=[h{1}(:) h{2}(:) h{3}(:)];
xN=cell(size(list,1),1);
for i=1:size(list,1)

    xm=[];
    xmt=list(i,:);
    tmp=xmt(xmt>0)';

    for iii=1:size(tmp,1)
        xm=[xm list(tmp(iii),:)']; %Each voxel needs a gradient value, so you need the surronding data for each one
    end

    xm(xm>0)=datamatrix(i,xm(xm>0));
    vec=nan(27,1);
    gd=sqrt(sum(((xm'*ym2).^2),2));
    vec(xmt>0)=gd;
    xN{i}=vec;
end
end