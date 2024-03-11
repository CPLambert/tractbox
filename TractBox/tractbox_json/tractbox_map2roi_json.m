function root = tractbox_map2roi_json
root.roi_ap.Description='Region of interest centroid location in the anterior-posterior plane (world space)'; 
root.roi_ml.Description='Region of interest centroid location in the medial-lateral plane (world space)'; 
root.roi_si.Description='Region of interest centroid location in the superior-inferior plane (world space)'; 
root.map_thresh_level.Description='Percentage threshold between 0 - 1 (i.e. 1 = 100%)';
root.map_thresh_val.Description='Corresponding value in the map image';
root.map_thresh_type.Description='Type of thresholding: Either highest (>thresh) or lowest (<thresh) values in the map';
root.map_vol.Description='Volume of voxels in the map surviving threshold';
root.map_ed.Description='Euclidean distance between ROI and thresholded map in mm';
root.map_val_mean.Description='Mean value within thresholded map';
root.map_val_std.Description='Standard deviation of values within thresholded map';
root.map_val_max.Description='Maximum value within thresholded map';
root.map_roi_overlap.Description='Percentage overlap between the thresholded map and the ROI of interest i.e. 100*(thresholded map volume/(intersection volume between thresholded map and roi))';
end