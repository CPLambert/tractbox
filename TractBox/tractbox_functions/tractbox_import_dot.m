function x0=tractbox_import_dot(dfile,thr)
%% function x0=tractbox_import_dot(dfile)
%Import a FSL dot file from probtrackx2 - Won't kill the memory as much or
%for as long as spconvert function.
%_________________________________________________________________________
% Version History:
% Version 1.1, December 2021, Code tidied up and formatted for tractbox
% Version 1.0, May 2012
%--------------------------------------------------------------------------
% C.Lambert - Wellcome Centre for Human Neuroimaging
%--------------------------------------------------------------------------

options = tractbox_defaults;

if nargin<2
    thr=options.dot.thr;
    if nargin<1
        disp('At least one input required');return
    end
end

if iscell(dfile),dfile=cell2mat(dfile);end

disp('Loading Data');fdt_matrix2=load(dfile);disp('DONE!');
disp(char(strcat('Importing Sparse Tractography Matrix from FSL .dot using PICo threshold=',num2str(thr))));

s=max(fdt_matrix2(:,1)); %Number of seed voxels
stran=max(fdt_matrix2(:,2)); %Number of brain voxels

%% Prepare matrices
x1   = sparse([],[],[],stran,s);
x0   = sparse([],[],[],stran,s);

%% Actual process
for i=1:s

    sub=fdt_matrix2((fdt_matrix2(:,1)==i),:); %Each seed
    sub=sub((sub(:,3)>=thr),:); %Only select voxels exceeding threshold

    if  ~isempty(sub)

        x1 = x1 + sparse(sub(:,2),i,sub(:,3),stran,s); %Fastest way to grow, in parts
        
        if ~rem(i,10),fprintf('.');end

        if ~rem(i,options.dot.step) || i==s
            x0   = x0 + x1; %Decant into bigger matrix every few iterations
            x1  = sparse([],[],[],stran,s);
            disp(char(strcat(num2str(i/s*100),'% done')));
        end
    end
end
end