function p = tractbox_distance_ed(C,steps)
%% function p = tractbox_distance_ed(C,steps)
%--------------------------------------------------------------------------
% Faster than MATLAB inbuilt pdist. Handles very large matrices by
% dividing matrix into blocks ("steps"), performing the calculation
% piecewise.
%_______________________________________________________________________
% Version History:
% Version 1.0, January 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

if nargin <2, steps = 1; end %Default to single pass calculation

Sq=size(C);
stepsize=round(Sq(1)/steps);
xa=cell(steps,1);xb=xa;

% Work out the coordinates and store seperate parts 
for iii=1:steps
    xa{iii}=((iii-1)*stepsize)+1;
    if iii~=steps
        xb{iii}=iii*stepsize;
    else
        xb{iii}=Sq(1);
    end
end

d = cell(steps,1);
s = zeros(steps,1);
for jjj=1:steps
    d{jjj} = sum(C(xa{jjj}:xb{jjj},:).^2,2); % Diagonal elements
    s(jjj) = size(C(xa{jjj}:xb{jjj},:),1); % Length of each block
end

% Prepare the ED vector (note, squareform(p) will convert to square,
% symmetric ED matrix

p = zeros(sum(s)*(sum(s)-1)/2,1);

for j1=1:steps
    for i1=j1:steps
        XX   = repmat(d{j1}',s(i1),1) +repmat(d{i1},1,s(j1)) - 2*(C(xa{i1}:xb{i1},:))*(C(xa{j1}:xb{j1},:))';
        ii    = (1:s(i1))'+sum(s(1:(i1-1)));
        jjj    = (1:s(j1)) +sum(s(1:(j1-1)));
        indp = repmat(ii,1,s(j1))+repmat((jjj-1)*(sum(s)-1)-jjj.*(jjj-1)/2,s(i1),1)-1;

        if i1==j1
            ind  = find(toeplitz([0 ones(1,s(i1)-1)],zeros(1,s(i1))));
            indp = indp(ind);
            p(indp) = sqrt(XX(ind));
        else
            p(indp) = sqrt(XX(:));
        end
    end
end
end