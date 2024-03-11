function rootjson = tractbox_connectivitystrength_json(cstype,thresh,source,dseg)
if strcmp(cstype,'adcon')
    rootjson.Description = 'Stength of connectivity calculated using normalised streamline density. This version takes a set percentage (csthreshold) of the maximum value to all the ROIs provided as an adaptive step (this is to account for differences in global headsize) and binarises the data at this threshold (defined in .tsv). The number of streamlines arriving at an ROI above this threshold are expressed as a ratio of ROI volume';
    rootjson.csthreshold = thresh;
end

if strcmp(cstype,'streamdensity')
    rootjson.Description = 'Stength of connectivity calculated using streamline density. This version takes a global fixed threshold and binarises the data at this. The number of streamlines arriving at an ROI above this threshold are expressed as a ratio of ROI volume';
    rootjson.csthreshold = thresh;
end

rootjson.Atlas = source.parcimg{1};
rootjson.AtlasDescription = source.parcdesc; %Probably not correct for BIDs, got annoyed figuring it out!
rootjson.dseg = dseg;
end