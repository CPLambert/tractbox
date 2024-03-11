function tractbox_map2roi(job)

options = tractbox_defaults;
tmp = job.met2roiinput.targetname;
tmp(isspace(tmp))='-';
filename = fullfile([tmp,'_',options.map2roi.label]);

if numel(job.met2roiinput.roiimg)~=numel(job.met2roiinput.mapimg)
    disp('Error - Uneven inputs');
end

for i=1:numel(job.met2roiinput.mapimg)
    root = tractbox_map2roi_struct;
    rootjson = tractbox_map2roi_json;
    [aa,bb,~]=fileparts(job.met2roiinput.mapimg{i});
    if contains(bb,'sub')
        dash=strfind(bb,'_');
        pre=bb(1:dash(1));
    else
        [~,pre,~] = fileparts(job.met2roiinput.roiimg{i});
        pre=char(strcat(pre,'_'));
    end

    opfile=fullfile(aa,[pre,filename,'.tsv']);
    opjson=fullfile(aa,[pre,filename,'.json']);
    spm_jsonwrite(opjson,rootjson,struct('indent','  '));

    % All measurements voxel to world. All
    % coordinates are internally consistent with MATLAB centroid command where
    % x = AP, y = ML and z = SI (which requires some re-mapping when working
    % out world coordinates, but all are re-ordered to the convention given
    % unless otherwise specified)
    V=nifti(job.met2roiinput.roiimg{i});Vi=V.dat(:,:,:)>0;
    Affine = V.mat;
    vx2vol        = prod(sqrt(sum(Affine(1:3,1:3).^2)));
    S = regionprops(Vi,'Centroid');tmp=[];
    tmp.ap=S.Centroid(1);tmp.ml=S.Centroid(2);tmp.si=S.Centroid(3);
    X = (Affine(1:3,1:3)*[tmp.ml tmp.ap tmp.si]')+Affine(1:3,4);
    X = [X(2,:);X(1,:);X(3,:)];
    root.roi_ap=X(1);root.roi_ml=X(2);root.roi_si=X(3);

    H=nifti(job.met2roiinput.mapimg{i});Ho=H.dat(:,:,:);
    tmp=(Ho(:));
    tmp(tmp==0)=[];%Due to the way heatmaps and metric maps are made, exact zeros values should not really be present within the map area - Bit of an assumptoin, could always add the mask option later
    root.map_thresh_level=job.met2roiinput.mapthresh;
    root.map_thresh_val=job.met2roiinput.mapthresh*(max(tmp)-min(tmp));

    if job.met2roiinput.threshtype
        root.map_thresh_type={'top'};
        Hi=Ho>root.map_thresh_val;
    else
        root.map_thresh_type={'bottom'};
        Hi=Ho<root.map_thresh_val;
        Hi(Ho==0)=0;
    end

    S = regionprops(Hi,'Centroid');tmp=[];
    tmp.ap=S(1).Centroid(1);tmp.ml=S(1).Centroid(2);tmp.si=S(1).Centroid(3);%biggest is first
    Y = (Affine(1:3,1:3)*[tmp.ml tmp.ap tmp.si]')+Affine(1:3,4);
    Y = [Y(2,:);Y(1,:);Y(3,:)];

    root.map_ap=Y(1);
    root.map_ml=Y(2);
    root.map_si=Y(3);

    R  = bsxfun(@minus,X,Y);

    root.map_ed  = sqrt(sum(R.^2,1))';
    root.map_vol = sum(Hi(:))*vx2vol;

    root.map_val_mean=mean(Ho(Hi==1));
    root.map_val_std=std(Ho(Hi==1));
    root.map_val_max=max(Ho(Hi==1));
    root.map_roi_overlap = ((sum(sum(sum(Hi.*Vi))))/sum(Hi(:)==1))*100;
    spm_save(opfile,root);
end
end
