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
%    * consider renaming files AFTER moving some folders into the "excluded" file.
%    * Acquistion of anatomical images NOTE: important for fmriprep ignore file (see below) TO-DO: ignore file is NOT working!
%
% USAGE NOTES:
%    * Make sure not include more than one T1w (e.g. TMS-localizer sequences)
%      into the 'anat' folder --> as fmriprep will use all of them (averaging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define which steps to do in this script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do.sortDicoms     = true; % Sort dicoms into folder in "sourcedata" (see function: sortDicomsIntoFolders.m)
do.excludeRuns    = true; % Put some functional runs into an "exclude" folder in sourcedata. Name of remaining runs stays the same
do.rawdataDirs    = true; % Create and pre-allocate new directories in "rawData"
do.datasetJson    = true;  % Create the dataset_description.json file (see function: BIDS_dataset_json.m);
do.scansTsv       = true;  % Create scans TSV file.
do.dicom2nifti    = true;  % Transform Dicoms into Niftis (from sourcedata to rawdata, see function: DICOMconversion_BIDS.m)
do.funcJson       = true;  % Add task and discarded images to func json files
do.fmapJson       = true;  % Add an "IntendedFor" field in the fieldmap Json file.
do.addBidsIgnore  = true;  % Add a BIDS-ignore file
do.addFprepIgnore = true;  % Add an fprepIgnore file
do.save           = true;  % Save the Workspace and command outputs in at rootDir/code/Dicom2Bids

if do.save
    % Everything is going to be save in a "log"-structure with sub-structures
    % consisting of "params", "subj_params" and "do".
    log.savefile      = ['BIDS_conversion_workspace' date];
    diary(log.savefile)
    diary on
end
disp('Running the following steps'); disp(do);



%% Define parameters for all scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjects = [5,8,11]; %#ok<*NBRAK> % Array with subject number (for Dicoms AND niftis). If EMPTY: if will get ALL subjects with param.prefix and RENAME them to sub-001, sub-002, etc!
n_subjects = length(subjects);

% Prefix of subjects in sourcedata folder
params.prefix     = 's';       % Prefix used in the subjects in dicoms ('PW', 's', 'p', etc)

% Directories
params.rootDir    = '/Volumes/pgrassi/projects/TMS-fMRI-piloting'; % Directory of the project. '/Volumes/DATA2/BIDS_test';
params.sourceDir  = fullfile(params.rootDir, 'sourcedata'); % Dir: sourcedata. DEFAULT: rootDir/sourcedata
params.rawDir     = fullfile(params.rootDir, 'rawdata');    % Dir: rawdata. DEFAULT:  rootDir/rawdata
params.anatDirs         = {'anat','findtms'};     % anatomical. DEFAULT: 'anat'. In case you have different forms anatomical scans in the dicoms folders, use this to import then to rawdata. NOTE: names given by sortDicomsIntoFolders.m
params.anatModalities   = {'T1w','T1w'};          % for BIDS conform naming _T1w. DEFAULT: 'T1w'.
params.anatAcquisition  = {'tmscoils','findtms'}; % one specific acquisition label for each modalities. DEFAULT: '1mm'.
params.funcDir    = 'func';     % functional. DEFAULT: 'func'
params.funcrefDir = 'func_ref'; % functional SB reference scans (MB sequences). DEFAULT: 'func_ref'
params.fmapDir    = 'fmap';     % fieldmaps. DEFAULT: 'fmap'
params.excludeDir = 'excluded'; % excluded functional runs. DEFAULT: 'excluded'

% Other specificiations
params.formatSpec    = '%03i';    % Format specification for the subject (e.g., 'sub-001').
params.formatSpecRun = '%03i'; % Format specification for the runs (e.g., 'run-01').
params.extensions    = {'**.IMA','**.ima'}; % for dicom images to find
params.runQuestions  = true; % true or false. With false: no questions asked! With true: some questions included in some crucial steps (e.g. nDummies, subject names, etc)

% Check if array 'subjects' is empty and make sure this is not an error
if isempty(subjects) && params.runQuestions
    warning('No specific subject entered. This will get ALL subjects with the defined prefix and rename them to sub-001, sub-002, etc');
    wait4confirmation = true;
    while wait4confirmation; tmp = input('Is this correct? y or n :', 's');
        if tmp == 'y'; wait4confirmation = false;
        elseif tmp == 'n'; error('Please define your subjects');
        end
    end
end

%% Define parameters for some scripts (what to do, etc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific parameters for the sort Dicoms into Folders
if do.sortDicoms
    params.steps2sort       = [true true true]; % Logical array: Three entries. Default = 1 1 1. STEP 1 = MAKE FOLDERS (01,...) , STEP 2 = RE-NAME FOLDER ('anat', ...), STEP 3 = RE-NAME SUBFOLDERS ('run-01', ...)
    params.dir2sequenceInfo = fullfile(params.rootDir, 'sequenceInfo.mat'); % dir to sequenceInfo DEFAULT: "sequenceInfo.mat"
    params.mkModalityDirs   = false; % Logical. creates modality specific sub-directories. DEFAULT = false.
end

% Specific parameters for dicom to nifti conversion
if do.dicom2nifti
    params.software           = 'dcm2niix'; % software to convert with: 'SPM' or 'dcm2niix'
    params.anatConversion     = true;      % anatomical conversion
    params.funcConversion     = true;       % functional conversion
    params.fieldmapConversion = true;       % fieldmap conversion
    params.SBref2nii          = true;       % conversion of single band reference images (only for multiband)
    params.path2exe           = '/Applications/MRIcroGL.app/Contents/Resources/'; % if using dcm2niix
    params.nDummies           = 5;          % number of dummy images. DEFAULT: none. If using SPM: dummy volumes are discarded, if using dcm2niix: no volume is discarded
end

%% Get subjects folders from source data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define subject names for rawdata (during conversion)
folders          = dir(fullfile(params.sourceDir,[params.prefix, '*'])); % Get subjects data from sourcedata
allSubNames      = {folders(:).name}; % all subNames (in sourcedata)
nSub             = length(folders);   % all folders in source ( = number of subjects)
subNames         = {};                % pre-allocate
rawSubNames      = {};                % pre-allocate

if isempty(folders); error('The specified directory does not contain any folders starting with the specified subject prefix'); end

if isempty(subjects) % If subjects is EMPTY, then go through all folders
    disp('NO SPECIFIC SUBJECT PROVIDED. RUNNING ALL SUBJECTS AND CREATING RAWDATA NAME AS SUB-001, SUB-002, etc');
    for i=1:nSub
        rawSubNames{end+1} = ['sub-' num2str(i,params.formatSpec)]; %#ok<*SAGROW> % make rawSubNames
    end
else
    % Get the corresponding index of the subjects with the index provided
    % by subjects and use these to convert only their data.
    disp('SUBJECTS PROVIDED. RUNNING SUBJECTS AND CREATING RAWDATA NAME AS SUB-001, SUB-002, etc');
    for i=1:length(subjects)
        tmpId     = contains(allSubNames,join(['s', num2str(subjects(i),params.formatSpec)]));
        if any(tmpId) % check if any tmpId matches
            rawSubNames{end+1} = ['sub-' num2str(subjects(i),params.formatSpec)];
        end
    end
end

if length(rawSubNames) ~= length(subjects); error('(Some) subjects entered were not found. Check values'); end

% Get subset of subjects. Important assuming that they are called: "PREFIX" +
% "NUMBER" (in formatSpec 00X). E.g.: s001, s002, etc.
for i = 1:length(subjects); subNames{i} = [params.prefix rawSubNames{i}(5:end)]; end
params.rawSubNames = rawSubNames;
params.subNames    = subNames;


%% Define subject specific parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdr         = {'filename', 'run' 'task'}; % header
tasks       = {'TMSlow','TMShigh'}; % which task

% Enter values per subject of interest
ss = 1; % subject index 1 from "subjects" array.
subj_params(ss).hdr         = hdr;
subj_params(ss).runs        = {{'2'},{'3'},{'4'},{'5'},{'6'},{'7'}}'; % which functional runs;
subj_params(ss).run2exclude = 1; % which runs to exclude (will be put into a separate folder
subj_params(ss).tasks       = {tasks{1},tasks{2},tasks{1},tasks{2},tasks{1},tasks{2}}'; % Low-high-low-high-low-high

ss = ss + 1;
subj_params(ss).hdr         = hdr;
subj_params(ss).runs        = {{'1'},{'2'},{'3'},{'5'},{'6'},{'7'}}'; % which functional runs;
subj_params(ss).run2exclude = 4; % which runs to exclude (will be put into a separate folder
subj_params(ss).tasks       = {tasks{2},tasks{1},tasks{2},tasks{1},tasks{2},tasks{1}}'; % Low-high-low-high-low-high

ss = ss + 1;
subj_params(ss).hdr         = hdr;
subj_params(ss).runs        = {{'1'},{'2'},{'3'},{'4'},{'5'},{'6'}}'; % which functional runs;
subj_params(ss).run2exclude = 7; % which runs to exclude (will be put into a separate folder
subj_params(ss).tasks       = {tasks{1},tasks{2},tasks{1},tasks{2},tasks{1},tasks{2}}'; % Low-high-low-high-low-high


for ss = 1:n_subjects
    for i = 1:length(subj_params(ss).runs)
        runs_filenames{i} = [params.funcDir '/' params.rawSubNames{ss} '_task-' subj_params(ss).tasks{1} '_run-' subj_params(ss).runs{i}{1} '_bold.nii.gz'];
    end
    subj_params(ss).filename = runs_filenames';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXECUTE                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SOURCEDATA: SORT DICOMS INTO FOLDERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort dicoms into folders (with folder names). Needs a sequenceInfo.mat
% file to know what sequence name belongs to what folder, number of scans
% and modalities. See sortDicomsIntoFolders for more info.
if do.sortDicoms
    sortDicomsIntoFolders(params.prefix, [], params);
end

%% PUT SOME RUNS INTO EXCLUDE FOLDER IN SOURCEDATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT: this will put some runs into an extra folder called for
% exclusion. This information needs to be given via the subj_params
% structure.
if do.excludeRuns
    disp('===========================');
    disp('EXCLUDING RUNS');
    disp('===========================');
    
    for i = 1:length(subjects) % for-loop across subjects
        if subj_params(i).run2exclude % If subject has runs to exclude
            currDir = fullfile(params.sourceDir, params.subNames{i}, params.funcDir); % directory: sourcedata/prefix*
            folders  = dir(fullfile(currDir,'run-*')); % get all runs
            
            if isempty(folders); error('No runs inside the func folder, but runs asked to be excluded: cannot continue');end % sanity check: there need to be runs
            spm_mkdir(params.sourceDir, params.subNames{i}, params.excludeDir); % Create excluded folder
            for tmprun = 1:length(folders) % Go throw runs
                for tmprun2exclude = 1:length(subj_params(ss).run2exclude) % and compare to the run to exclude
                    if contains(folders(tmprun).name,num2str(subj_params(i).run2exclude(tmprun2exclude), '%02i'))
                        movefile(fullfile(folders(tmprun).folder, folders(tmprun).name), fullfile(params.sourceDir, params.subNames{i}, params.excludeDir));
                        if params.SBref2nii % If SB reference present: also take those to the 'exclude' folder
                            movefile(fullfile(params.sourceDir, params.subNames{i}, params.funcrefDir,folders(tmprun).name), fullfile(params.sourceDir, params.subNames{i}, params.excludeDir, 'SBrefs'));
                        end
                    end
                end
            end
            % TO-DO: consider renaming files AFTER moving some folders into the
            % "excluded" file. OR NOT (because of the TSV file)
        end
    end
end


%% RAWDATA: Make BIDS-conform folders in rawdata directory for all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the rawdata folders (for niftis) (unprocessed nifti data folders)
% using spm_mkdir (anatomical, functional folders and fieldmap folders)

if do.rawdataDirs
    disp('===========================');
    disp('CREATING RAWDIRS');
    disp('===========================');
    spm_mkdir(params.rawDir, params.rawSubNames, params.anatDirs);
    spm_mkdir(params.rawDir, params.rawSubNames, params.funcDir);
    spm_mkdir(params.rawDir, params.rawSubNames, params.fmapDir);
    %Uncomment to generate also special folders for the anatomical modalities:
    %spm_mkdir(params.rawDir, params.rawSubNames, params.anatDir, params.anatModalities);
end

%% DATASET JSON: Make json dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the dataset_description.json file
if do.datasetJson
    disp('===========================');
    disp('CREATING DATASET JSON FILE');
    disp('===========================');
    BIDS_dataset_json(params.rawDir);
end

%% SCAN TSV: Make *_scans.tsv for all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create *_scans.tsv file for your functional runs with task name and run for each subject
% The file contains run numbers the corresponding task
% Define tsv columns
if do.scansTsv
    disp('===========================');
    disp('CREATING SCAN TSV FILE');
    disp('===========================');
    
    for i = 1:length(subjects)
        currDir = fullfile(params.rawDir, params.rawSubNames{i}); % directory: rawdata/sub-*
        tmpfileName = [currDir, '/sub-' num2str(subjects(i),params.formatSpec) '_scans']; % filename: rawdata/sub-*/sub-*_scans
        tbl = cell2table(cat(2,subj_params(i).filename, subj_params(i).runs,subj_params(i).tasks));
        tbl.Properties.VariableNames = subj_params(i).hdr;
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
% Add an "IntendedFor" field in the fieldmap Json file.
% NOTE: Implemented for gre_fieldmap (in which one has 3 images).
% TO-DO: Implement pepolar (and image with other polarity!)
% See: https://bids-specification.readthedocs.io/en/latest/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#using-intendedfor-metadata
if do.fmapJson
    disp('===========================');
    disp('MODIFING FIELDMAP JSON FILES');
    disp('===========================');
    
    for i = 1:length(subjects) % for-loop across subjects
        
        % Intended For Functional scans:
        funcDir  = fullfile(params.rawDir, params.rawSubNames{i}, params.funcDir); % directory: rawdata/sub-*/func
        tmpfunc  = dir([funcDir '/*.gz']);
        tmpfunc  = join([repmat(join([params.funcDir "/"],''),length(tmpfunc),1),string({tmpfunc(1:length(tmpfunc)).name}')],''); % add the folder "func/" to the run_names
        
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
    disp('===========================');
    disp('ADDING BIDS IGNORE FILE');
    disp('===========================');
    tmpstr = ["notes/", "fmriprep_BIDS_filter.json", "**/*_scans*", "findtms/", "excluded/"];
    tmpf   = fopen([params.rawDir '/.bidsignore'],'w');
    fprintf(tmpf,'%s\n', tmpstr);
    tmpf   = fclose(tmpf);
end


%% FMRI-PREP-IGNORE JSON: Generate a json file to ignore certain images in fmri prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do.addFprepIgnore
    disp('===========================');
    disp('ADDING FMRI-PREP IGNORE FILE');
    disp('===========================');
    BIDs_filter.t1w.datatype = 'anat';
    BIDs_filter.t1w.suffix   = params.anatModalities{1};  % this corresponds to T1w
    BIDs_filter.t1w.acquisition = '64ch';                 % this corresponds to 1mm
    jsonwrite([params.rawDir '/fmriprep_BIDS_filter.json'], BIDs_filter, 'prettyPrint', true);
end

%% FMRI-PREP: Now run fmri-prep using Docker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO-DO


%% SAVE: parameters and further infos.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do.save
    log.params = params;
    log.params.subjects = subjects;
    log.do   = do;
    log.date = date;
    log.savefolder = fullfile(params.rootDir, 'Code','Dicom2Bids');
    save([log.savefolder '/' log.savefile], 'log');
    diary off
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



