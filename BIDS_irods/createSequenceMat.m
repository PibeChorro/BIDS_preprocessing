%% Create Sequence INFO:
fprintf('=========CREATING sequenceInfo.mat FILE:  \n\n')

%% General settings
Headers = {'Name','Description','Nr_vols','Modality'};

%% First session sequences
logs.ses(1).tbl = cell2table(cell(0,4),'VariableNames', Headers);

% Anatomical: add the new sequence to table
newSequence={"anat", "ADNI_192slices_64channel_me", {192}, "T1w"}; %#ok<*STRSCALR
logs.ses(1).tbl = [logs.ses(1).tbl; newSequence];

% functional: Add the new sequence to table
newSequence={"func", "ep2d_3mm_TE30_64channel_MB2_GRAPPA2_NOGAP_FIRST_SESSION", {[230,331]}, "bold"}; 
logs.ses(1).tbl = [logs.ses(1).tbl; newSequence];

% fieldmaps: Add the new sequence to table
newSequence={"fmap", "gre_field_mapping", {[]}, "bold"}; %#ok<*CLARRSTR> 
logs.ses(1).tbl = [logs.ses(1).tbl; newSequence];

%% Second session sequences

% Display the table
disp(logs.ses(1).tbl);

save('sequenceInfo.mat',"logs");

% Before:
%logs.ses(tmpses).sequenceDescriptions    = {}; % e.g. "AFNI_..."
%logs.ses(tmpses).sequenceNames           = {}; % e.g. "anat"
%logs.ses(tmpses).sequenceScanNrs         = {}; % e.g. "192"
%logs.ses(tmpses).sequenceModality        = {}; % e.g. "T1w" (ONLY IMPORTANT FOR ANATOMICAL SCANS)
