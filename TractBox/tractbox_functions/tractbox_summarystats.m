function results=tractbox_summarystats(root)
%% results=tractbox_summarystats(root)
% Very simple function to run a basic data summary for continious variables
% This will return the following array:
%
% ROI = [mean | std | median | iqr | normal | overlap]
%
% mean - Mean, omitting NaNs
% std - Standard Deviation
% median - Median
% iqr - Intraquartile Range
% normal - Anderson Darling Normality Test, 1 is not normal
% overlap - Population overlap
%--------------------------------------------------------------------------
% Version History:
% Version 1.0, November 2022
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------
warning('off');

fx=fields(root);

for i = 1:numel(fx)
output=[];
    input=root.(fx{i});

if sum(~isnan(input))>4
    output.normal = adtest(input); %1 is not normal
else
    output.normal = nan;
end
output.mean = mean(input,'omitnan'); %Use these if normal ==0
output.std = std(input,'omitnan');
output.median = median(input,'omitnan'); %normal == 1
output.iqr = iqr(input);
output.overlap = 100*(sum(input>0)./numel(input));

results.(fx{i}) = [output.mean output.std output.median output.iqr output.normal output.overlap];
end

end