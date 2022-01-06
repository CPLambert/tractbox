function list = tractbox_find_neighbours(seedcoord,dim)
%% function list = tractbox_find_neighbours(seedcoord,dim)
% This will search through a list of seed coordinates from image size dim
% and identify those in the immediate vicinity (i.e. 3D connected neighbourhood). 
% 
% By default it will map a neighbourhood size of 3, which corresponds to the three-dimensional
% 26-connected neighborhood (plus central voxel = 27) - In theory bigger neighbourhoods
% are possible but this has not been added yet. It will then search through the seedcoord to check if
% data for these "connections" exist (e.g. absent edges). It will return a
% coordinate list (seeds x 27) that maps to the corresponding
% coordinate given in seedcoord. This allows neighbourhood data from each seed 
% to be easily extracted from other sources (e.g. Euclidean distance squareform matrix)
%_______________________________________________________________________
% Version History:
% Version 1.0, January 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

options = tractbox_defaults;
ksize = 3;

%Prepare matrix:
list=zeros(size(seedcoord,1),ksize^3);

% How to step through matrix
step=median(1:1:ksize)-1;
jo=(step.*-1);oj=jo;Sj=jo;
opo=-1;
jump=(dim(1)*dim(2))*jo;
nudge=opo*dim(1);

% Work out all the positions
for i=1:ksize^3
    list(:,i)=seedcoord+jump+nudge+oj;
    oj=oj+1;

    if ~rem(i,ksize)
        oj=jo;
        opo=opo+1;
        nudge=opo*dim(1);
    end

    if ~rem(i,ksize^2) %Central voxel
        Sj=Sj+1;
        jump=(dim(1)*dim(2))*Sj;
        opo=-1;
        nudge=opo*dim(1);
    end
end

%Zero impossible values:
list(list<0)=0;

% Loop through all the seeds and only retain neighbours that exist in data
for j=1:size(seedcoord,1)

    if ~rem(j,1000) %If very big
        disp(char(strcat(num2str((j/size(seedcoord,1)*100),'% Complete'))));
    end

    if sum(logical(seedcoord(j)==list(:)))
        list(list==seedcoord(j))=j.*-1;
    end
end

list(list>=0)=0;
list=list.*-1;
end