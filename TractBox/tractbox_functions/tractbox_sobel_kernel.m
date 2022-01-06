function h = tractbox_sobel_kernel
%% function h = tractbox_sobel_kernel
% Return 3D sobel kernel for MRI edge detection. Default 3 x 3. At some
% point generalise this for bigger kernel sizes
%_______________________________________________________________________
% Version History:
% Version 1.0, January 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

h{1}(:,:,1)=[-1,-2,-1;-2,-4,-2;-1,-2,-1];
h{1}(:,:,3)=h{1}(:,:,1).*-1;
h{2}=permute(h{1},[2,3,1]);
h{3}=permute(h{1},[3,2,1]);
end

             