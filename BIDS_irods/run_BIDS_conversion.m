%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT: run_BIDS_conversion
%
% This script starts with unsorted DICOM images and transforms them into
% NIFTI format with a BIDS-compatible structure. Thus, it will create json
% sidecar files, rename niftis, etc.
%
% DEPENDENCIES:
%    * SPM (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/):
%               for DICOM 2 nifti conversion (optional) and file control
%    * dcm2niix (https://github.com/rordenlab/dcm2niix):
%               for DICOM 2 nifti conversion (optional) and JSON file generation
%    * JSONio (https://github.com/gllmflndn/JSONio):
%               for reading and writing json files in matlab
%    * and secondary matlab functions. Specially:
%         - sortDicomIntoFolders.m (STEP 1)
%         - DICOMconversion_BIDS.m (STEP 6)
%
% IMPORTANT: please read the help of sortDicomsIntoFolder.m and DICOMSconversion_BIDS.m
%            to run the pipeline without problems.
%
% GENERAL NOTES ABOUT THE USE:
% The script consists of a PREPARATION and an EXECUTION part. Each
% consisting of different steps. The user NEEDS to define parameters in the
% PREPARATION PART.
%
% PREPARATION PART:
%    * STEP 1: Define WHAT TO DO (in structure: "do")
%    * STEP 2: Define PARAMETERS for all steps (in structure: "params")
%    * STEP 3: Define PARAMETERS for specific steps (in structure: "params")
%    * STEP 4: Define SUBJECT SPECIFIC PARAMETERS (in structure: "subj_params")
%
% EXECUTION PART:
%    * STEP 1: do.sortDicoms
%              Sort dicoms into folders: --> using function: sortDicomsIntoFolders.m
%    * STEP 2: do.excludeRuns
%              Exclude runs in a subject-based manner (e.g., exclude run 3 from s001, and run 5 from s006).
%              Individualised definition is necessary (using structure "subj_params")
%    * STEP 3: do.rawDataDirs
%              Creates 'rawdata' directories.
%    * STEP 4: do.datasetJson
%              Creates a dataset_description.json file
%    * STEP 5: do.scanTsv
%              Creates a scan TSV file (IMPORTANT: it includes ONLY the non-excluded
%              functional scans from STEP 2). Individualised definition is necessary
%              (using structure "subj_params")
%    * STEP 6: do.dicom2nifti
%              Convert dicoms to niftis: --> using function: DICOMconversion_BIDS
%              Uses either dcm2niix OR SPM. If using SPM dummy images are
%              discarded. If using dcm2niix dummy images are NOT discarded.
%              Importantly: individual dummy images are NOT supported (for
%              now).
%    * STEP 7: do.funcJson
%              Create functional json files AFTER dicom2nifti conversion
%    * STEP 8: do.fmapJson
%              Modify fieldmap json files AFTER dicom2nifti conversion
%    * EXTRA STEPS HAPPEN HERE
%    * STEP 9: do.save
%              Save workspace, and diary. (final step)
%
% EXTRAS STEPS: Add "ignore" files.
%    *  do.addBidsIgnore  = true; % Add a BIDS-ignore file
%    *  do.addFprepIgnore = true; % Add an fprepIgnore file
%
% USAGE NOTES:
%    * Main steps are STEP 1 (sort dicom into folders) and STEP 6
%    (transform dicoms into niftis).
%    * Several steps REQUIRE specific parameters, make sure these are
%    defined adequately.
%    * The subj_params ARE needed for: STEP 2 (run exclusion) and STEP 5 (scan-tsv
%    creation)
%
%    * Make sure not include more than one T1w (e.g. TMS-localizer sequences)
%      into the 'anat' folder --> as fmriprep will use all of them (averaging)
%    * PEPOLAR images for fieldmap correction should be called: fmap_pepolar
%
% TO-DO:
%    * Optimize the Sort Dicoms into folders script.
%    * Finish writing this script, save values, etc.
%    * consider renaming files AFTER moving some folders into the "excluded" file.
%    * Acquistion of anatomical images NOTE: important for fmriprep ignore file (see below) TO-DO: ignore file is NOT working!
%    * Consider making the dummy images subject and/or run dependent (as
%    sometimes differences CAN occur).
%
% Original PRG: 11/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean up and check for dependencies
clear do logs params subj_params;

if which('spm'); disp('SPM is in path'); else; error('SPM is not in path'); end
if which('jsonread'); disp('JSONio is in path'); else; error('JSONio is not in path'); end
if which('sortNiftisIntoRaw.m'); disp('sortNiftisIntoRaw.m is in path'); else; error('sortNiftisIntoRaw.m is not in path'); end

%% Preparation STEP 1: Define which steps to do in this script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main steps in EXECUTION:
do.rawdataDirs    = true;  % Create and pre-allocate new directories in "rawData"
do.scansTsv       = true;  % Create scans TSV file.
do.loadSeq        = true;  % load sequenceInfo or create a new one.
do.sortNiftis     = true;  % Sort niftis into folders in "rawdata" (see function: sortNiftisIntoRaw.m)
do.datasetJson    = true; % Create the dataset_description.json file (see function: BIDS_dataset_json.m);
do.funcJson       = true; % Add task and discarded images to func json files
do.fmapJson       = true;  % Add an "IntendedFor" field in the fieldmap Json file.
do.save           = true; % Save the Workspace and command outputs in at rootDir/code/Dicom2Bids

% Extras:
do.addBidsIgnore  = true; % Add a BIDS-ignore file
do.addFprepIgnore = true; % Add an fprepIgnore file

%% Define main directories (rootDir and save Dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.rootDir    = '/Volumes/bartels_data/TMS-WM-on/WM2_fMRI'; % Directory of the project. '/Volumes/DATA2/BIDS_test';
params.saveDir    = 'Dicom2Bids'; % rootDir/code/Dicom2Bids;

% Now display what to do, start the diary if wanted.
if do.save
    % Everything is going to be save in a "logs"-structure with sub-structures
    % consisting of "params", "subj_params" and "do".
    logs.savefile      = ['BIDS_conversion_workspace_' date];
    if ~isfolder(fullfile(params.rootDir, 'Code', params.saveDir))
        mkdir(fullfile(params.rootDir, 'Code', params.saveDir));
    end
    diary(fullfile(params.rootDir, 'Code', params.saveDir, logs.savefile));
    diary on
end
disp('Running the following steps'); disp(do);

%% Preparation STEP 2: Define parameters for all scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    subjects = [3,4,5]; %#ok<*NBRAK> % Array with subject number (for Dicoms AND niftis). If EMPTY: if will get ALL subjects with param.prefix and RENAME them to sub-001, sub-002, etc!
    n_subjects = length(subjects);

    % Prefix of subjects in sourcedata folder
    params.prefix      = 'sub-'; % Prefix used in the subject folder ('pwm', 'PW', 's', 'p', etc)
    params.formatSpecSource  = '%03i'; % Format specification for the subject from source data (e.g., 'prefix001').

    % Further directories
    params.sourceDir  = fullfile(params.rootDir, 'sourcedata'); % Dir: sourcedata. DEFAULT: rootDir/sourcedata
    params.rawDir     = fullfile(params.rootDir, 'rawdata');    % Dir: rawdata. DEFAULT:  rootDir/rawdata
    params.sesDirs    = {'ses-01','ses-02'}; % Cell - name of each session
    params.anatDirs         = {{'anat'},{'anat','findtms'}}; % Cell, entry per sesssion. anatomical. DEFAULT: 'anat'. In case you have different forms anatomical scans in the dicoms folders, use this to import them to rawdata. NOTE: names given by sortNiftisIntoFolders.m
    params.anatModalities   = {{'T1w'},{'T1w','T1w'}};  % Cell, entry per session. for BIDS conform naming _T1w. DEFAULT: 'T1w'.
    params.anatAcquisition  = {{'1mm'},{'tmscoils','findtms'}};  % Cell, entry per session. One specific acquisition label for each modalities. DEFAULT: '1mm'.
    params.funcDir    = 'func';     % functional. DEFAULT: 'func'
    params.funcrefDir = 'func_ref'; % functional SB reference scans (MB sequences). DEFAULT: 'func_ref'
    params.fmapDir    = 'fmap';     % fieldmaps. DEFAULT: 'fmap'
    params.excludeDir = 'excluded'; % excluded functional runs. DEFAULT: 'excluded'
    params.sequenceInfoName = 'sequenceInfo.m';

    % Other specificiations
    params.formatSpec    = '%03i'; % Format specification for the subject (e.g., 'sub-001').
    params.formatSpecRun = '%03i'; % Format specification for the runs (e.g., 'run-01').
    params.extensions    = {'*.nii.gz','*.json'}; % for nifti images to find
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

    %%  Preparation STEP 3: Define parameters for some steps (what to do in which function, etc)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Specific parameters for the sort Niftis into Folders
    params.nDummies = 5; % number of dummy images. DEFAULT: none.

    %% Preparation STEP 4: Get subjects folders from source data
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
            tmpId     = contains(allSubNames,join([params.prefix, num2str(subjects(i),params.formatSpecSource)]));
            if any(tmpId) % check if any tmpId matches
                rawSubNames{end+1} = ['sub-' num2str(subjects(i),params.formatSpec)];
            end
        end
    end

    if length(rawSubNames) ~= length(subjects); error('(Some) subjects entered were not found. Check values'); end

    % Get subset of subjects. Important assuming that they are called: "PREFIX" +
    % "NUMBER" (in formatSpec 00X). E.g.: s001, s002, sub-010, etc.
    for i = 1:length(subjects); subNames{i} = [params.prefix num2str(subjects(i),params.formatSpecSource)]; end
    params.rawSubNames = rawSubNames; % in rawdata
    params.subNames    = subNames;    % in sourcedata

    %% Preparation STEP 5: Define subject specific parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Important: this definitions are REQUIRED for STEP 2 (exclude runs) and
    % STEP 5 (create tsv scan file).
    hdr         = {'filename', 'run' 'task'}; % header
    ses1_tasks  = {'WM2','ROIloc'};       % which task in session 1
    ses2_tasks  = {'TMSlow','TMShigh'};   % which task in session 2

    % Enter values per subject of interest
    %sub-003
    ss = 1; % subject INDEX from "subjects" array, rawSubNames and subNames; sub-003
    subj_params(ss).rawname     = params.rawSubNames{ss};
    subj_params(ss).sourcename  = params.subNames{ss};
    subj_params(ss).ses2run     = 1; % index from session in param.sesDirs % 1 == session 1, 2 == session 2
    subj_params(ss).hdr         = hdr;
    subj_params(ss).ses(1).runs     = {{'4'},{'5'},{'6'},{'7'},{'8'},{'9'},{'10'},{'11'}}'; % which functional runs to trasnfer and rename.
    subj_params(ss).ses(1).nrImages = [230, 230, repmat(331,1,6)]; % If known: add this to sanity check the number of images (not necessary).
    subj_params(ss).ses(1).run2exclude = []; % which runs to exclude (will be put into a separate folder). Leave empty is no run should be excluded.
    subj_params(ss).ses(1).tasks       = {ses1_tasks{2},ses1_tasks{2},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1}}';
    params.ses2run{ss}          = subj_params(ss).ses2run; % add to params for other scripts

    % sub-004
    ss = ss+1; % subject INDEX from "subjects" array, rawSubNames and subNames; 
    subj_params(ss).rawname     = params.rawSubNames{ss};
    subj_params(ss).sourcename  = params.subNames{ss};
    subj_params(ss).ses2run     = 1; % index from session in param.sesDirs % 1 == session 1, 2 == session 2
    subj_params(ss).hdr         = hdr;
    subj_params(ss).ses(1).runs     = {{'4'},{'5'},{'6'},{'7'},{'8'},{'9'},{'10'},{'11'}}'; % which functional runs to trasnfer and rename.
    subj_params(ss).ses(1).nrImages = [230, 230, repmat(331,1,6)]; % If known: add this to sanity check the number of images (not necessary).
    subj_params(ss).ses(1).run2exclude = []; % which runs to exclude (will be put into a separate folder). Leave empty is no run should be excluded.
    subj_params(ss).ses(1).tasks       = {ses1_tasks{2},ses1_tasks{2},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1}}';
    params.ses2run{ss}          = subj_params(ss).ses2run; % add to params for other scripts

    % sub-005
    ss = ss+1; % subject INDEX from "subjects" array, rawSubNames and subNames; 
    subj_params(ss).rawname     = params.rawSubNames{ss};
    subj_params(ss).sourcename  = params.subNames{ss};
    subj_params(ss).ses2run     = 1; % index from session in param.sesDirs % 1 == session 1, 2 == session 2
    subj_params(ss).hdr         = hdr;
    subj_params(ss).ses(1).runs     = {{'4'},{'5'},{'6'},{'7'},{'9'},{'10'},{'11'},{'12'}}'; % which functional runs to trasnfer and rename.
    subj_params(ss).ses(1).nrImages = [230, 230, repmat(331,1,6)]; % If known: add this to sanity check the number of images (not necessary).
    subj_params(ss).ses(1).run2exclude = []; % which runs to exclude (will be put into a separate folder). Leave empty is no run should be excluded.
    subj_params(ss).ses(1).tasks       = {ses1_tasks{2},ses1_tasks{2},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1}}';
    params.ses2run{ss}          = subj_params(ss).ses2run; % add to params for other scripts

    if length(subj_params) ~= n_subjects; error('Wrong number of subjects specified/wrong number of subjects defined'); end

    for ss = 1:n_subjects
        for sesid = 1:length(subj_params(ss).ses2run)
            sesid2 = subj_params(ss).ses2run(sesid);
            % Sanity check: control that the included/excluded runs have no
            % intersection.
            tmpruns = cellfun(@str2num, cat(1,subj_params(ss).ses(sesid2).runs{:})');
            if any(ismember(subj_params(ss).ses(sesid2).run2exclude, tmpruns)); error('Subject %s run inclusion/exclusion definition is not consistent',subj_params(ss).sourcename); end

            for i = 1:length(subj_params(ss).ses(sesid2).runs)
                runs_filenames{i} = [params.funcDir '/' params.rawSubNames{ss} '_' params.sesDirs{sesid2} '_task-' subj_params(ss).ses(sesid2).tasks{i} '_run-' sprintf('%03d',str2num(subj_params(ss).ses(sesid2).runs{i}{1})) '_bold.nii.gz']; %#ok<*ST2NM>
            end
            subj_params(ss).ses(sesid2).filename = runs_filenames';
        end
    end

catch ME
    disp('Preparation ended with errors');
    if do.save
        logs.error = ME;
        logs.savefolder = fullfile(params.rootDir, 'code', params.saveDir);
        save([logs.savefolder '/' logs.savefile '_error'], 'logs');
        diary off
        disp('Saving...');
    end
    rethrow(ME);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXECUTE                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try % To execute steps, else catch ME
    %% Run STEP 0: RAWDATA: Make BIDS-conform folders in rawdata directory for all subjects
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the rawdata folders (for niftis) (unprocessed nifti data folders)
    % using spm_mkdir (anatomical, functional folders and fieldmap folders)
    % separately per session
    if do.rawdataDirs
        disp('===========================');
        disp('CREATING RAWDIRS');
        disp('===========================');
        spm_mkdir(params.rawDir, params.rawSubNames, params.sesDirs, params.funcDir);
        spm_mkdir(params.rawDir, params.rawSubNames, params.sesDirs, params.fmapDir);
        for i = 1:length(params.sesDirs)
            spm_mkdir(params.rawDir, params.rawSubNames, params.sesDirs{i}, params.anatDirs{i});
        end
        %Uncomment to generate also special folders for the anatomical modalities:
        %spm_mkdir(params.rawDir, params.rawSubNames, params.sesDirs, params.anatDir, params.anatModalities);
    end

    %% Run STEP 1: SCAN TSV: Make *_scans.tsv for all subjects
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create *_scans.tsv file for your functional runs with task name and run for each subject
    % The file contains run numbers the corresponding task.
    % IMPORTANT: it will ONLY include non-excluded runs! (from STEP 2)
    % Define tsv columns
    if do.scansTsv
        disp('===========================');
        disp('CREATING SCAN TSV FILE');
        disp('===========================');

        for i = 1:length(subjects)
            for sesid = 1:length(subj_params(ss).ses2run)
                sesid2 = subj_params(ss).ses2run(sesid);
                currDir = fullfile(params.rawDir, params.rawSubNames{i}, params.sesDirs{sesid2}); % directory: rawdata/sub-*
                tmpfileName = [currDir, '/sub-' num2str(subjects(i),params.formatSpec), '_' params.sesDirs{sesid2}, '_scans']; % filename: rawdata/sub-*/sub-*_ses-*_scans
                tbl = cell2table(cat(2,subj_params(i).ses(sesid2).filename, subj_params(i).ses(sesid2).runs,subj_params(i).ses(sesid2).tasks));
                tbl.Properties.VariableNames = subj_params(i).hdr;
                writetable(tbl,tmpfileName,'Delimiter','\t');
                movefile([tmpfileName '.txt'],[tmpfileName '.tsv'],'f');
            end
        end
    end

    %% Run STEP 2: FROM SOURCEDATA INTO RAWDATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort Niftis into rawdata, also put some runs into an extra folder
    % called "exclusion". This information needs to be given via the
    % subj_params structure.
    if do.sortNiftis
        %Check if there already exists a sequence description mat file
        if exist([params.sequenceInfoName],'file')
            disp('Loading sequenceInfo.mat from directory');
            load(params.sequenceInfoName,'logs');
        else
            disp('No sequenceInfo! Will be created by running createSequenceMat.m');
            try disp('Creating new sequenceInfo!'); createSequenceMat;
            catch; error('createSequenceMat.m is not in file. No sequence info is created.');
            end
        end
        sortNiftisIntoRaw(params.prefix, logs, params);
    end

    %% Run STEP 3: DATASET JSON: Make json dataset
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the dataset_description.json file
    if do.datasetJson
        disp('===========================');
        disp('CREATING DATASET JSON FILE');
        disp('===========================');
        BIDS_dataset_json(params.rawDir);
    end

    %% Run STEP 4: FUNCTIONAL JSON: Add task and discarded images to func json files
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
            for sesid = 1:length(subj_params(ss).ses2run)
                sesid2 = subj_params(ss).ses2run(sesid);

                currDir = fullfile(params.rawDir, params.rawSubNames{i}, params.sesDirs{sesid2}, params.funcDir); % directory: rawdata/sub-*
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
    end

    %% Run STEP 5: FMAP JSON: Add "IntendedFor" and B0FieldIdentifier metadata for the Fieldmap scans.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add an "IntendedFor" field in the fieldmap Json file.
    % NOTE: Implemented for gre_fieldmap (in which one has 3 images).
    % TO-DO: Implement pepolar (and image with other polarity!)
    % See: https://bids-specification.readthedocs.io/en/latest/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#using-intendedfor-metadata
    % See for sessions: https://neurostars.org/t/intendedfor-field-needs-to-point-to-an-existing-file/1884/2
    if do.fmapJson
        disp('===========================');
        disp('MODIFING FIELDMAP JSON FILES');
        disp('===========================');
        for i = 1:length(subjects) % for-loop across subjects
            for sesid = 1:length(subj_params(ss).ses2run)
                sesid2 = subj_params(ss).ses2run(sesid);
                % Intended For Functional scans:
                funcDir  = fullfile(params.rawDir, params.rawSubNames{i}, params.sesDirs{sesid2}, params.funcDir); % directory: rawdata/sub-*/func
                tmpfunc  = dir([funcDir '/' params.prefix '*.gz']);
                tmpfunc  = join([repmat(join([params.sesDirs{sesid2} "/" params.funcDir "/"],''),length(tmpfunc),1),string({tmpfunc(1:length(tmpfunc)).name}')],''); % add the folder "func/" to the run_names

                currDir = fullfile(params.rawDir, params.rawSubNames{i}, params.sesDirs{sesid2}, params.fmapDir); % directory: rawdata/sub-*/fmap
                tmpjson = dir([currDir '/sub*.json']);
                if isempty(tmpjson)
                    warning('\nNo fieldmap found in folder: %s. Skipping step\n', currDir);
                elseif length(tmpjson) ~= 3
                    warning('\nWrong number of fieldmaps in folder: %s. Skipping step. \n Expected number of fieldmaps: 3 \n Found fieldmaps: %s\n', currDir, num2str(length(tmpjson)));
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
    end

    %% Run EXTRA 1: CREATE A BIDSIGNORE FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ignore specific folders (such as tmsloc) by creating a .bidsignore file
    if do.addBidsIgnore
        disp('===========================');
        disp('ADDING BIDS IGNORE FILE (.bidsignore file)');
        disp('===========================');
        tmpstr = ["**/fmriprep_BIDS_filter_ses1.json", "**/fmriprep_BIDS_filter_ses2.json", "**/*_scans*", "**/*findtms*", "**/*excluded*", "**/legacy"];
        fprintf(['BIDS ignore file includes:  %s'], join(tmpstr,', ')); % display what it includes
        tmpf   = fopen([params.rawDir '/.bidsignore'],'w');
        fprintf(tmpf,'%s\n', tmpstr);
        tmpf   = fclose(tmpf);
    end


    %% Run EXTRA 2: FMRI-PREP-IGNORE JSON: Generate a json file to ignore certain images in fmri prep
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do.addFprepIgnore

        %fMRIPrep Ignore for session 1
        %fMRIPrep Ignore for session 2
        disp('===========================');
        disp('CREATING FMRI-PREP FILTER FILE');
        disp('===========================');
        BIDs_filter.t1w.datatype = "anat";
        BIDs_filter.t1w.suffix   = "T1w";     % this corresponds to T1w
        BIDs_filter.t1w.acquisition = "1mm";  % this corresponds to 1mm
        BIDs_filter.t1w.session = "01";  % this corresponds to 1mm
        BIDs_filter.bold.datatype = "func";
        BIDs_filter.bold.session  = "01";
        BIDs_filter.bold.suffix   = "bold";
        BIDs_filter.fmap.datatype = "fmap";
        BIDs_filter.fmap.session  = "01";
        jsonwrite([params.rawDir '/fmriprep_BIDS_filter_ses1.json'], BIDs_filter, 'prettyPrint', true);


        %fMRIPrep Ignore for session 2
        disp('===========================');
        disp('CREATING FMRI-PREP FILTER FILE');
        disp('===========================');
        BIDs_filter.t1w.datatype = "anat";
        BIDs_filter.t1w.suffix   = "T1w";     % this corresponds to T1w
        BIDs_filter.t1w.acquisition = "1mm";  % this corresponds to 1mm
        BIDs_filter.t1w.session = "01";       % this corresponds to 1mm
        BIDs_filter.bold.datatype = "func";
        BIDs_filter.bold.session  = "02";
        BIDs_filter.bold.suffix   = "bold";
        BIDs_filter.fmap.datatype = "fmap";
        BIDs_filter.fmap.session  = "02";
        jsonwrite([params.rawDir '/fmriprep_BIDS_filter_ses2.json'], BIDs_filter, 'prettyPrint', true);

    end

    %% Run FINAL STEP: SAVE: parameters and further infos.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do.save
        logs.params = params;
        logs.params.subjects = subjects;
        logs.do   = do;
        logs.date = date; %#ok<*DATE>
        logs.savefolder = fullfile(params.rootDir, 'Code',params.saveDir);

        if ~isfolder(logs.savefolder); spm_mkdir(logs.savefolder); end % make dir if necessary

        save([logs.savefolder '/' logs.savefile], 'logs');
        diary off
    end

catch ME
    disp('Execution ended with errors');
    if do.save
        logs.error = ME;
        logs.savefolder = fullfile(params.rootDir, 'Code', params.saveDir);
        save([logs.savefolder '/' logs.savefile '_error'], 'logs');
        diary off
        disp('Saving...');
    end
    rethrow(ME);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



