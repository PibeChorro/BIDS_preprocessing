%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT: run_BIDS_conversion
% This script will:
%     * Sort dicoms into folders:
%           --> using function: sortDicomsIntoFolders.m
%     * Create BIDS-conform folder structure (source and rawdata)
%     * Create a dataset.json file
%     * Create a _scan.tsv file
%     * Convert dicoms to niftis:
%           --> using function: DICOMconversion_BIDS
%     * Create functional json files.
%     * Modify fieldmap json files.
%     * Create a bids-ignore file
%
% TO-DO:
%    * Optimize the Sort Dicoms into folders script.
%    * Finish writing this script, save values, etc.
%    * Add conversion of SBref and fieldmaps in the DICOMconvert script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define what to do:
do.sortDicoms     = true;
do.rawdataDirs    = false;
do.datasetJson    = false;
do.scansTsv       = false;
do.dicom2nifti    = false;
do.funcJson       = false;
do.fmapJson       = false;
do.addBidsIgnore  = false;
do.addFprepIgnore = false;

%% Define parameters for all scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjects = 7;

params.rootDir    = '/Volumes/pgrassi/projects/TMS-fMRI-piloting'; % rootDir
params.formatSpec = '%03i'; % format specification for the subject. The formatSpec for runs is set as default to '%02i'
params.prefix     = 's';    % for subjects in dicoms
params.extensions = {'**.IMA','**.ima'}; % for images to find
params.anatModalitiesDirs   = {'T1w','tmsloc'};      % name of modalities directories in the dicoms folders
params.anatModalities       = {'T1w','T1w'};         % for BIDS conform naming _T1w (do not include tmsloc into the T1w-anat folder --> fmriprep will use all of them!)
params.anatAcquisition      = {'tmscoils','tmsloc'}; % one specific acquisition label for each modalities. NOTE: important for fmriprep ignore file (see below) TO-DO: ignore file is NOT working! 
params.anatDir    = 'anat'; % anatomical
params.funcDir    = 'func'; % functional
params.fmapDir    = 'fmap'; % fieldmaps
params.nDummies   = 5;      % number of dummy images

params.sourceDir = fullfile(params.rootDir, 'sourcedata'); % sourcedata
params.rawDir    = fullfile(params.rootDir, 'rawdata'); % rawdata


% Specific parameters for dicom to nifti conversion
params.software           = 'dcm2niix';      % software to convert with: 'SPM' or 'dcm2niix'
params.funcConversion     = false;    % functional conversion
params.anatConversion     = false;    % anatomical conversion
params.fieldmapConversion = false;    % fieldmap conversion
params.path2exe           = '/Applications/MRIcroGL.app/Contents/Resources/'; % if using dcm2niix

% Define subject names for rawdata (during conversion)
folders          = dir(fullfile(params.sourceDir,[params.prefix, '*'])); % Get subjects data from sourcedata
allSubNames      = {folders(:).name}; % all subNames
nSub             = length(folders);   % all folders in source
subNames         = {};                % pre-allocate
rawSubNames      = {};                % pre-allocate

if isempty(folders); error('The specified directory does not contain any folders starting with the specified subject prefix');
end

if isempty(subjects)
    for i=1:nSub
        rawSubNames{end+1} = ['sub-' num2str(i,params.formatSpec)]; %#ok<*SAGROW>
    end
else
    % Get the corresponding index of the subjects with the index provided
    % by subjects and use these to convert only their data.
    for i=1:length(subjects)
        tmpId     = contains(allSubNames,num2str(subjects(i),params.formatSpec));
        if tmpId
            rawSubNames{end+1} = ['sub-' num2str(subjects(i),params.formatSpec)];
        end
    end
end

for i = 1:length(subjects)
    subNames{i} = [params.prefix rawSubNames{i}(5:end)];
end

params.rawSubNames = rawSubNames;
params.subNames    = subNames;


%% SOURCEDATA: SORT DICOMS INTO FOLDERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do.sortDicoms
    sortDicomsIntoFolders(params.prefix, [], params);
end

%% RAWDATA: Make BIDS-conform folders in rawdata directory for all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the rawdata folders (for niftis) (unprocessed nifti data folders)
% using spm_mkdir (anatomical and functional folders)

if do.rawdataDirs
    spm_mkdir(params.rawDir, params.rawSubNames, params.anatDir);
    spm_mkdir(params.rawDir, params.rawSubNames, params.funcDir);
    %Uncomment to generate also special folders for the anatomical modalities:
    %spm_mkdir(params.rawDir, params.rawSubNames, params.anatDir, params.anatModalities);
end

%% DATASET JSON: Make json dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the dataset_description.json file
if do.datasetJson
    BIDS_dataset_json(params.rawDir);
end

%% SCAN TSV: Make *_scans.tsv for all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create *_scans.tsv file for your functional runs with task name and run for each subject
% The file contains run numbers the corresponding task
% TO-DO: include file name
% Define tsv columns
if do.scansTsv
    disp('===========================');
    disp('CREATING SCAN TSV FILE');
    disp('===========================');
    
    hdr         = {'filename', 'run' 'task'}; % header
    runs        = {'1','2','3','4','5','6'}'; % which functional runs
    tasks       = {'TMSlow','TMShigh'}; % which task
    
    
    % Enter values per subject of interest
    ss = 1; % subject index 1 from "subjects" array.
    vals(ss).hdr   = hdr;
    vals(ss).runs  = runs;
    vals(ss).tasks = {tasks{1},tasks{2},tasks{1},tasks{2},tasks{1},tasks{2}}'; % Low-high-low-high-low-high
    vals(ss).filename    = {[params.funcDir '/' params.rawSubNames{ss} '_task-' vals(ss).tasks{1} '_run-' runs{1} '_bold.nii.gz'];
        [params.funcDir '/' params.rawSubNames{ss} '_task-' vals(ss).tasks{2} '_run-' runs{2} '_bold.nii.gz'];
        [params.funcDir '/' params.rawSubNames{ss} '_task-' vals(ss).tasks{3} '_run-' runs{3} '_bold.nii.gz'];
        [params.funcDir '/' params.rawSubNames{ss} '_task-' vals(ss).tasks{4} '_run-' runs{4} '_bold.nii.gz'];
        [params.funcDir '/' params.rawSubNames{ss} '_task-' vals(ss).tasks{5} '_run-' runs{5} '_bold.nii.gz'];
        [params.funcDir '/' params.rawSubNames{ss} '_task-' vals(ss).tasks{6} '_run-' runs{6} '_bold.nii.gz']};
    
    for i = 1:length(subjects)
        currDir = fullfile(params.rawDir, params.rawSubNames{i}); % directory: rawdata/sub-*
        tmpfileName = [currDir, '/sub-' num2str(subjects(i),params.formatSpec) '_scans']; % filename: rawdata/sub-*/sub-*_scans
        tbl = cell2table(cat(2,vals(i).filename, vals(i).runs,vals(i).tasks));
        tbl.Properties.VariableNames = vals(i).hdr;
        writetable(tbl,tmpfileName,'Delimiter','\t');
        movefile([tmpfileName '.txt'],[tmpfileName '.tsv'],'f');
    end
end

%% CONVERT: DICOM TO NIFTI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform from dicom to niftis using the function DICOMconversion_BIDS
% Use either 'SPM' or 'dcm2niix'
if do.dicom2nifti
    disp('===========================');
    disp('RUNNING DICOM CONVERSION');
    disp('===========================');
    DICOMconversion_BIDS(rawSubNames, subNames, 'params', params, 'funcConversion', params.funcConversion, 'anatConversion', params.anatConversion, 'fieldmapConversion', params.fieldmapConversion)
end

%% FUNCTIONAL JSON: Add task and discarded images to func json files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After generating the functional json files: we want to add the TaskName
% to them. Instead of using the defined task names defined before, we will
% extract them from the name of the json file.
% We separate the json files into elements (separator = '_') and then again
% to get the task name after "task-*".
% DEFINE:
%        * position of "task" in func filename
% NOTE:
%        * modify accordingly if necessary
if do.funcJson
    disp('===========================');
    disp('MODIFING FUNCTIONAL JSON FILES');
    disp('===========================');
    params.positionOfTaskInName = 2;
    params.NumberOfVolumesDiscardedByUser = 0; % important: are you discarding dummies before including files in the dataset? If so, how many?
    
    for i = 1:length(subjects) % for-loop across subjects
        currDir = fullfile(params.rawDir, params.rawSubNames{i}, params.funcDir); % directory: rawdata/sub-*
        tmpjson = dir([currDir '/*.json']);
        for ii = 1:length(tmpjson) % for-loop across json files in functional folder
            tmpstr = regexp(tmpjson(ii).name,'_','split'); % get all elements separated by '_'
            tmpstr = regexp(tmpstr{params.positionOfTaskInName},'-','split'); %separate the "task-TASKNAME" in two elements
            try
                tmpvals = jsonread([tmpjson(ii).folder '/' tmpjson(ii).name]); % read json file
                tmpvals.TaskName = tmpstr{2}; % add taskname
                tmpvals.NumberOfVolumesDiscardedByUser = params.NumberOfVolumesDiscardedByUser;
                jsonwrite([tmpjson(ii).folder '/' tmpjson(ii).name], tmpvals, 'prettyPrint', true) % write json file again.
            catch
                warning('%s\n%s\n%s\n%s', ...
                    'Reading/Writing the JSON file seems to have failed.', ...
                    'Make sure that the following library is in the matlab/octave path:', ...
                    'https://github.com/gllmflndn/JSONio');
            end
        end
    end
end

%% FMAP JSON: Add "IntendedFor" and B0FieldIdentifier metadata for the Fieldmap scans.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This will add an "IntendedFor" field in the fieldmap Json file.
% See: https://bids-specification.readthedocs.io/en/latest/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#using-intendedfor-metadata
if do.fmapJson
    disp('===========================');
    disp('MODIFING FIELDMAP JSON FILES');
    disp('===========================');
    
    for i = 1:length(subjects) % for-loop across subjects
        
        % Intended For Functional scans:
        funcDir  = fullfile(params.rawDir, params.rawSubNames{i}, params.funcDir); % directory: rawdata/sub-*/fund
        tmpfunc  = dir([funcDir '/*.gz']);
        tmpfunc  = join([repmat("func/",6,1),string({tmpfunc(1:length(tmpfunc)).name}')],''); % add the folder func/ to the run_names
        
        currDir = fullfile(params.rawDir, params.rawSubNames{i}, params.fmapDir); % directory: rawdata/sub-*/fmap
        tmpjson = dir([currDir '/sub*.json']);       
        if isempty(tmpjson)
            warning('\nNo fieldmap found in folder. Skipping step\n');
        elseif length(tmpjson) ~= 3
            warning('\nWrong number of fieldmaps in folder. Skipping step. \n Expected number of fieldmaps: 3 \n Found fieldmaps: %s\n', num2str(length(tmpjson)));
        else
            for ii = 1:length(tmpjson) % for-loop across json files in fieldmap folder
                try
                    tmpvals             = jsonread([tmpjson(ii).folder '/' tmpjson(ii).name]); % read json file
                    tmpvals.IntendedFor = tmpfunc; % add taskname
                    jsonwrite([tmpjson(ii).folder '/' tmpjson(ii).name], tmpvals, 'prettyPrint', true) % write json file again.
                catch
                    warning('%s\n%s\n%s\n%s', ...
                        'Reading/Writing the JSON file seems to have failed.', ...
                        'Make sure that the following library is in the matlab/octave path:', ...
                        'https://github.com/gllmflndn/JSONio');
                end
            end
            
        end
    end
end


%% CREATE A BIDSIGNORE FILE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ignore specific folders (such as tmsloc) by creating a .bidsignore file
if do.addBidsIgnore
    tmpstr = ["notes/", "fmriprep_BIDS_filter.json", "**/*_scans*", "tmsloc/"];
    tmpf   = fopen([params.rawDir '/.bidsignore'],'w');
    fprintf(tmpf,'%s\n', tmpstr);
    tmpf   = fclose(tmpf);
end


%% SAVE: parameters and further infos.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO-DO



%% FMRI-PREP-IGNORE JSON: Generate a json file to ignore certain images in fmri prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do.addFprepIgnore
    BIDs_filter.t1w.datatype = 'anat';
    BIDs_filter.t1w.suffix   = params.anatModalities{1};  % this corresponds to T1w
    BIDs_filter.t1w.acquisition = '64ch';                 % this corresponds to 1mm
    jsonwrite([params.rawDir '/fmriprep_BIDS_filter.json'], BIDs_filter, 'prettyPrint', true);
end

%% FMRI-PREP: Now run fmri using Docker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO-DO






