%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT: run_BIDS_conversion
%
% This script starts with unsorted DICOM images and transforms them into
% NIFTI format with a BIDS-compatible structure. Thus it will create json
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
if which('sortDicomsIntoFolders'); disp('sortDicomsIntoFolders is in path'); else; error('sortDicomsIntoFolders is not in path'); end
if which('DICOMconversion_BIDS'); disp('DICOMconversion_BIDS is in path'); else; error('DICOMconversion_BIDS is not in path'); end

%% STEP 1: Define which steps to do in this script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main steps:
do.sortDicoms     = false; % Sort dicoms into folder in "sourcedata" (see function: sortDicomsIntoFolders.m)
do.excludeRuns    = false;  % Put some functional runs into an "exclude" folder in sourcedata. Name of remaining runs stays the same
do.rawdataDirs    = false;  % Create and pre-allocate new directories in "rawData"
do.datasetJson    = false;  % Create the dataset_description.json file (see function: BIDS_dataset_json.m);
do.scansTsv       = true;  % Create scans TSV file.
do.dicom2nifti    = false;  % Transform Dicoms into Niftis (from sourcedata to rawdata, see function: DICOMconversion_BIDS.m)
do.funcJson       = false; % Add task and discarded images to func json files
do.fmapJson       = false; % Add an "IntendedFor" field in the fieldmap Json file.
do.save           = false; % Save the Workspace and command outputs in at rootDir/code/Dicom2Bids

% Extras:
do.addBidsIgnore  = false; % Add a BIDS-ignore file
do.addFprepIgnore = false; % Add an fprepIgnore file

%% Define main directories (rootDir and save Dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.rootDir    = '/Volumes/bartels_data/pgrassi/TMS-fMRI-WM2'; % Directory of the project. '/Volumes/DATA2/BIDS_test';
params.saveDir    = 'Dicom2Bids'; % rootDir/Code/Dicom2Bids;

% Now display what to do, start the diary if wanted.
if do.save
    % Everything is going to be save in a "logs"-structure with sub-structures
    % consisting of "params", "subj_params" and "do".
    logs.savefile      = ['BIDS_conversion_workspace_' date]; %#ok<DATE>
    diary(fullfile(params.rootDir, 'Code', params.saveDir, logs.savefile));
    diary on
end
disp('Running the following steps'); disp(do);

%% STEP 2: Define parameters for all scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    subjects = [1]; %#ok<*NBRAK> % Array with subject number (for Dicoms AND niftis). If EMPTY: if will get ALL subjects with param.prefix and RENAME them to sub-001, sub-002, etc!
    n_subjects = length(subjects);

    % Prefix of subjects in sourcedata folder
    params.prefix      = 'pwm'; % Prefix used in the subjects in dicoms ('PW', 's', 'p', etc)
    params.formatSpecSource  = '%02i'; % Format specification for the subject from source data (e.g., 'prefix001').

    % Further directories
    params.sourceDir  = fullfile(params.rootDir, 'sourcedata'); % Dir: sourcedata. DEFAULT: rootDir/sourcedata
    params.rawDir     = fullfile(params.rootDir, 'rawdata');    % Dir: rawdata. DEFAULT:  rootDir/rawdata
    params.sesDirs    = {'ses-01','ses-02'};
    params.anatDirs         = {{'anat'},{'anat','findtms'}}; % anatomical. DEFAULT: 'anat'. In case you have different forms anatomical scans in the dicoms folders, use this to import them to rawdata. NOTE: names given by sortDicomsIntoFolders.m
    params.anatModalities   = {{'T1w'},{'T1w','T1w'}};  % for BIDS conform naming _T1w. DEFAULT: 'T1w'.
    params.anatAcquisition  = {{'1mm'},{'tmscoils','findtms'}};  % one specific acquisition label for each modalities. DEFAULT: '1mm'.
    params.funcDir    = 'func';     % functional. DEFAULT: 'func'
    params.funcrefDir = 'func_ref'; % functional SB reference scans (MB sequences). DEFAULT: 'func_ref'
    params.fmapDir    = 'fmap';     % fieldmaps. DEFAULT: 'fmap'
    params.excludeDir = 'excluded'; % excluded functional runs. DEFAULT: 'excluded'

    % Other specificiations
    params.formatSpec    = '%03i'; % Format specification for the subject (e.g., 'sub-001').
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

    %% STEP 3: Define parameters for some steps (what to do in which function, etc)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Specific parameters for the sort Dicoms into Folders
    if do.sortDicoms
        params.steps2sort       = [true true true]; % Logical array: Three entries. Default = 1 1 1. STEP 1 = MAKE FOLDERS (01,...) , STEP 2 = RE-NAME FOLDER ('anat', ...), STEP 3 = RE-NAME SUBFOLDERS ('run-01', ...)
        params.dir2sequenceInfo = []; %fullfile(params.rootDir, 'sequenceInfo.mat'); % dir to sequenceInfo DEFAULT: "sequenceInfo.mat"
        params.mkModalityDirs   = false; % Logical. creates modality specific sub-directories. DEFAULT = false.
    end

    % Specific parameters for dicom to nifti conversion
    if do.dicom2nifti
        params.software           = 'dcm2niix'; % software to convert with: 'SPM' or 'dcm2niix'
        params.anatConversion     = true;       % anatomical conversion
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
            tmpId     = contains(allSubNames,join([params.prefix, num2str(subjects(i),params.formatSpecSource)]));
            if any(tmpId) % check if any tmpId matches
                rawSubNames{end+1} = ['sub-' num2str(subjects(i),params.formatSpec)];
            end
        end
    end

    if length(rawSubNames) ~= length(subjects); error('(Some) subjects entered were not found. Check values'); end

    % Get subset of subjects. Important assuming that they are called: "PREFIX" +
    % "NUMBER" (in formatSpec 00X). E.g.: s001, s002, etc.
    for i = 1:length(subjects); subNames{i} = [params.prefix num2str(subjects(i),params.formatSpecSource)]; end
    params.rawSubNames = rawSubNames; % in rawdata
    params.subNames    = subNames;    % in sourcedata

    %% STEP 4: Define subject specific parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Important: this definitions are REQUIRED for STEP 2 (exclude runs) and
    % STEP 5 (create tsv scan file).
    hdr         = {'filename', 'run' 'task'}; % header
    ses1_tasks  = {'WM2','ROIloc'}; %{'TMSlow','TMShigh'};       % which task

    % Enter values per subject of interest
    ss = 1; % subject index from "subjects" array, rawSubNames and subNames;
    subj_params(ss).rawname     = params.rawSubNames{ss};
    subj_params(ss).sourcename  = params.subNames{ss};
    subj_params(ss).ses2run     = 1; % index from session in param.sesDirs
    subj_params(ss).hdr         = hdr;
    subj_params(ss).ses(1).runs        = {{'3'},{'4'},{'5'},{'6'},{'7'},{'8'},{'9'},{'10'}}'; % which functional runs;
    subj_params(ss).ses(1).run2exclude = [1, 2]; % which runs to exclude (will be put into a separate folder)
    subj_params(ss).ses(1).tasks       = {ses1_tasks{1},ses1_tasks{1},ses1_tasks{2},ses1_tasks{1},ses1_tasks{1},ses1_tasks{1},ses1_tasks{2},ses1_tasks{2}}';
    params.ses2run{ss}          = subj_params(ss).ses2run; % add to params for other scripts

    if length(subj_params) ~= n_subjects; error('Wrong number of subjects specified/wrong number of subjects defined'); end

    for ss = 1:n_subjects
        for sesid = 1:length(subj_params(ss).ses2run)
            sesid2 = subj_params(ss).ses2run(sesid);
            % Sanity check: control that the included/excluded runs have no
            % intersection. TO-DO alternative: first define runs to INCLUDE, and
            % THEN define the excluded runs based on that variable, then there is
            % no room for mistakes.
            tmpruns = cellfun(@str2num, cat(1,subj_params(ss).ses(sesid2).runs{:})');
            if any(ismember(subj_params(ss).ses(sesid2).run2exclude, tmpruns)); error('Subject %s run inclusion/exclusion definition is not consistent',subj_params(ss).sourcename); end

            for i = 1:length(subj_params(ss).ses(sesid2).runs)
                runs_filenames{i} = [params.funcDir '/' params.rawSubNames{ss} '_' params.sesDirs{sesid2} '_task-' subj_params(ss).ses(sesid2).tasks{i} '_run-' sprintf('%03d',str2num(subj_params(ss).ses(sesid2).runs{i}{1})) '_bold.nii.gz'];
            end
            subj_params(ss).ses(sesid2).filename = runs_filenames';
        end
    end

catch ME
    disp('Preparation ended with errors');
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
% EXECUTE                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try % To execute steps, else catch ME
    %% STEP 1: SOURCEDATA: SORT DICOMS INTO FOLDERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort dicoms into folders (with folder names). Needs a sequenceInfo.mat
    % file to know what sequence name belongs to what folder, number of scans
    % and modalities. See sortDicomsIntoFolders for more info.
    if do.sortDicoms
        sortDicomsIntoFolders(params.prefix, [], params);
    end

    %% STEP 2: PUT SOME RUNS INTO EXCLUDE FOLDER IN SOURCEDATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IMPORTANT: this will put some runs into an extra folder called for
    % exclusion. This information needs to be given via the subj_params
    % structure.
    if do.excludeRuns
        disp('===========================');
        disp('EXCLUDING RUNS');
        disp('===========================');

        for i = 1:length(subjects) % for-loop across subjects
            for sesid = 1:length(subj_params(ss).ses2run)
                sesid2 = subj_params(ss).ses2run(sesid);

                if subj_params(i).ses(sesid2).run2exclude % If subject has runs to exclude
                    currDir  = fullfile(params.sourceDir, params.subNames{i}, params.sesDirs{sesid2}, params.funcDir); % directory: sourcedata/prefix*
                    folders  = dir(fullfile(currDir,'*run-*')); % get all runs

                    if isempty(folders); error('No runs inside the func folder, but runs asked to be excluded: cannot continue');end % sanity check: there need to be runs
                    spm_mkdir(params.sourceDir, params.subNames{i},params.sesDirs{sesid2},params.excludeDir); % Create excluded folder
                    for tmprun = 1:length(folders) % Go throw runs
                        for tmprun2exclude = 1:length(subj_params(ss).ses(sesid2).run2exclude) % and compare to the run to exclude
                            if contains(folders(tmprun).name,['run-' num2str(subj_params(i).ses(sesid2).run2exclude(tmprun2exclude), params.formatSpecRun)])
                                movefile(fullfile(folders(tmprun).folder, folders(tmprun).name), fullfile(params.sourceDir, params.subNames{i},params.sesDirs{sesid2}, params.excludeDir));

                                %if params.SBref2nii % If SB reference present: also take those to the 'exclude' folder
                                %    movefile(fullfile(params.sourceDir, params.subNames{i},params.sesDirs{sesid2},params.funcrefDir,folders(tmprun).name), fullfile(params.sourceDir, params.subNames{i}, params.sesDirs{sesid2}, params.excludeDir, 'SBrefs'));
                                %end

                            end
                        end
                    end
                    % TO-DO: consider renaming files AFTER moving some folders into the
                    % "excluded" file. OR NOT (because of the TSV file)
                end
            end
        end
    end

    %% STEP 3: RAWDATA: Make BIDS-conform folders in rawdata directory for all subjects
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

    %% STEP 4: DATASET JSON: Make json dataset
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the dataset_description.json file
    if do.datasetJson
        disp('===========================');
        disp('CREATING DATASET JSON FILE');
        disp('===========================');
        BIDS_dataset_json(params.rawDir);
    end

    %% STEP 5: SCAN TSV: Make *_scans.tsv for all subjects
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

    %% STEP 6: CONVERT: DICOM TO NIFTI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transform from dicom to niftis using the function DICOMconversion_BIDS
    % Use either 'SPM' or 'dcm2niix'
    if do.dicom2nifti
        disp('===========================');
        disp('RUNNING DICOM CONVERSION');
        disp('===========================');
        DICOMconversion_BIDS(params.rawSubNames, params.subNames, params.ses2run, 'params', params, 'funcConversion', params.funcConversion, 'anatConversion', params.anatConversion, 'fieldmapConversion', params.fieldmapConversion)
    end

    %% STEP 7: FUNCTIONAL JSON: Add task and discarded images to func json files
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

    %% STEP 8: FMAP JSON: Add "IntendedFor" and B0FieldIdentifier metadata for the Fieldmap scans.
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
            for sesid = 1:length(subj_params(ss).ses2run)
                sesid2 = subj_params(ss).ses2run(sesid);
                % Intended For Functional scans:
                funcDir  = fullfile(params.rawDir, params.rawSubNames{i}, params.sesDirs{sesid2}, params.funcDir); % directory: rawdata/sub-*/func
                tmpfunc  = dir([funcDir '/*.gz']);
                tmpfunc  = join([repmat(join([params.funcDir "/"],''),length(tmpfunc),1),string({tmpfunc(1:length(tmpfunc)).name}')],''); % add the folder "func/" to the run_names

                currDir = fullfile(params.rawDir, params.rawSubNames{i}, params.sesDirs{sesid2}, params.fmapDir); % directory: rawdata/sub-*/fmap
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
    end

    %% EXTRA 1: CREATE A BIDSIGNORE FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ignore specific folders (such as tmsloc) by creating a .bidsignore file
    if do.addBidsIgnore
        disp('===========================');
        disp('ADDING BIDS IGNORE FILE (.bidignore file)');
        disp('===========================');
        tmpstr = ["notes/", "fmriprep_BIDS_filter.json", "**/*_scans*", "findtms/", "excluded/"];
        fprintf(['BIDS ignore file includes:  %s'], join(tmpstr,', ')); % display what it includes
        tmpf   = fopen([params.rawDir '/.bidsignore'],'w');
        fprintf(tmpf,'%s\n', tmpstr);
        tmpf   = fclose(tmpf);
    end


    %% EXTRA 2: FMRI-PREP-IGNORE JSON: Generate a json file to ignore certain images in fmri prep
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do.addFprepIgnore
        disp('===========================');
        disp('ADDING FMRI-PREP IGNORE FILE (JSON file)');
        disp('===========================');
        BIDs_filter.t1w.datatype = 'anat';
        BIDs_filter.t1w.suffix   = params.anatModalities{1};  % this corresponds to T1w
        BIDs_filter.t1w.acquisition = '64ch';                 % this corresponds to 1mm
        jsonwrite([params.rawDir '/fmriprep_BIDS_filter.json'], BIDs_filter, 'prettyPrint', true);
    end

    %% EXTRA 3: FMRI-PREP: Now run fmri-prep using Docker
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TO-DO: include fmri-prep here or not?


    %% STEP 9: SAVE: parameters and further infos.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do.save
        logs.params = params;
        logs.params.subjects = subjects;
        logs.do   = do;
        logs.date = date;
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



