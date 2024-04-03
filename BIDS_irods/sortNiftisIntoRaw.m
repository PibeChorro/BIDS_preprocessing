function sortNiftisIntoRaw(prefix, logs, params)


%% SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 1)
    clc;
    warning('No settings given!')
    help sortNiftisIntoRaw
    return
end

fprintf('\n====SORTING NIFTIS INTO RAWFOLDERS====\n');

%% Get values or set defaults from params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get/enter rootDir
if isfield(params, 'rootDir'); rootDir = params.rootDir; else
    fprintf(['====ROOTDIR: No default. Please select your project folder \n (ideally it should contain a folder named "sourcedata")\n\n'])
    rootDir = uigetdir(pwd, 'Select Project Folder [rootDir]');
end
if rootDir == 0; error('No folder was selected. Re-start'); end

% Check source and raw dirs
if isfield(params, 'sourceDir') && isfield(params, 'rawDir')
    sourceDir = params.sourceDir;
    rawDir    = params.rawDir;
else
    fprintf(['====SOURCE AND RAW DIRS: No dirs specified. Using defaults name: sourcedata and rawdata.\n\n']);
    sourceDir = fullfile(rootDir, 'sourcedata');  % sourcedata
    rawDir    = fullfile(rootDir, 'rawdata');     % rawdata
end

% Check sessions
if isfield(params, 'sesDirs')
    sesDirs = params.sesDirs;
else
    fprintf(['====SESSION DIRS: No session dir specified. Using defaults name: ses-01 \n\n']);
    sesDirs = {'ses-01'};  % sourcedata
end

% Get/enter "ses2run" per subject
if isfield(params, 'ses2run'); ses2run = params.ses2run; else
    ses2run = {[1]};
    fprintf('====SESSION to run: Missing. Using a default session of 1)\n\n');
end

% Extension we care about
if isfield(params, 'extensions'); extensions = params.extensions; else; extensions = {'*.nii.gz','*.json'}; % extension you care about
    fprintf(['====EXTENSIONS: No extensions specified. Using default extensions: .nii.gz and .json \n\n']);
end

% Run some specific checks?
if isfield(params, 'runQuestions'); runQuestions = params.runQuestions; else; runQuestions = true; % ask questions or not, run some checks
    fprintf(['====RUNNING CHECKS: Missing. Default: true. Running checks \n\n']);
end

% Anat
if isfield(params, 'anatModalities') && isfield(params, 'anatAcquisition') && isfield(params, 'anatDirs') % specify the anatomical modalities to use and acquisition values
    anatModalities      = params.anatModalities;
    anatDirs            = params.anatDirs;
    anatAcquisition     = params.anatAcquisition;
else
    fprintf(['====ANATMODALITIES & ACQUISITION & DIRS: No anatModalities/acquisition/dirs specified. Using default: anat, T1w and 1 mm for the first session.\n\n']);
    anatModalities      = {{'T1w'}};
    anatDirs            = {{'anat'}};
    anatAcquisition     = {{'1mm'}};
end

% Functional
if isfield(params, 'funcDir'); funcDir = params.funcDir; else; funcDir = 'func';
    fprintf(['====FUNC DIR: No dir specified. Using default name: func.\n\n']);
end %#ok<*NBRAK>

%if isfield(params, 'funcrefDir'); funcrefDir = params.funcrefDir; else; funcrefDir = 'func_ref';
%    fprintf(['====FUNC REF DIR: No dir specified. Using default name: func_ref.\n\n']);
%end
%if isfield(params, 'SBref2nii'); SBref2nii = params.SBref2nii; else; SBref2nii  = false; % NOT including the Single Band reference
%    fprintf(['====SB reference for multiband images: Using Default. NOT INCLUDING ANY SINGLE BAND REFERENCE. \n']);
%end

% Fieldmap
if isfield(params, 'fmapDir'); fmapDir = params.fmapDir; else; fmapDir = 'fmap';
    fprintf(['====FMAP DIR: No dir specified. Using default name: fmap\n\n']);
end

% Format specification for RUNS
if isfield(params, 'formatSpecRun'); formatSpecRun = params.formatSpecRun; else; formatSpecRun  = '%03i'; % specify how you name your RUNS
    fprintf(['====FORMATSPEC: No formatSpectRuns specified. Using default formatSpec for RUNS %03i. \n\n']);
end

%% Get subjects in dir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE: prefix of subject folders. Make sure that the prefix is unique for
% your subject folders and you do not have any other files starting with
% this prefix.
if isempty(prefix)
    prefix = input (['Please specify the prefix of your participant data.\n' ...
        '(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');
    fprintf('\n\n');
end

sub = dir(fullfile(sourceDir,[prefix '*'])); % This will give you also the number of subjects = number of folders.
if isempty(sub); error('The specified directory does not contain any folders starting with the specified prefix'); end

% Sanity check: check if folder and no hidden folder
areFolders = all([sub(:).isdir]);
areHidden  = any(contains({sub(:).name},'.'));
if ~areFolders; error('Your prefix seems not to be unique for folders')
elseif areHidden; error('It looks like you did not select a valid prefix, because a hidden folder was selected.')
end

% If specific subjects provided: sort only those!
if isfield(params, 'rawSubNames')
    tmpsub = sub(1); % just pre-allocate one subject
    for tmps1 = 1:length(params.rawSubNames)
        for tmps2 = 1:length(sub)
            if strcmp(sub(tmps2).name,params.rawSubNames{tmps1}); tmpsub(tmps1) = sub(tmps2); end
        end
    end
    sub = tmpsub; % overwrite the sub structure with only those sub
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BRING NIFTIS FROM SOURCEDATA INTO RAW DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f = 1:length(sub) % Outer loop for the subjects
    fprintf('Sorting subject: %s \n', sub(f).name); % display.
    for s2r = 1:length(ses2run{f}); tmp_ses2run{s2r} = sesDirs{ses2run{f}(s2r)}; end %#ok<AGROW> % get the right session to run in format "ses-01"

    for tmpses = 1:length(tmp_ses2run) % Loop for the session
        %% Set folders
        subjectDir = fullfile(sourceDir,sub(f).name, tmp_ses2run{tmpses},'DICOM_NIFTI/'); % specific for session
        subjrawDir = fullfile(rawDir,sub(f).name,tmp_ses2run{tmpses});

        %% Read TSV file, get task names, runs and number of images
        if isfile(fullfile(subjrawDir,[sub(f).name '_' tmp_ses2run{tmpses} '_scans.tsv']))
            movefile(fullfile(subjrawDir,[sub(f).name '_' tmp_ses2run{tmpses} '_scans.tsv']),fullfile(subjrawDir,[sub(f).name '_' tmp_ses2run{tmpses} '_scans.txt']),'f');
            [tmpTSV]= readtable(fullfile(subjrawDir,[sub(f).name '_' tmp_ses2run{tmpses} '_scans.txt']));
            movefile(fullfile(subjrawDir,[sub(f).name '_' tmp_ses2run{tmpses} '_scans.txt']),fullfile(subjrawDir,[sub(f).name '_' tmp_ses2run{tmpses} '_scans.tsv']),'f');
        else
            error('No *_scans.tsv file found');
        end

        % Get task names:
        taskNames = tmpTSV.task;
        runs = tmpTSV.run;

        %% Get session we are in
        if strcmp(tmp_ses2run(tmpses),'ses-01')
            tmpsesID = 1;
        elseif strcmp(tmp_ses2run(tmpses),'ses-02')
            tmpsesID = 2;
        end

        %% Get images and jsons from source dir
        tmpImgSource  = dir(fullfile(subjectDir,['0' extensions{1}])); % niftis (note that the niftis and json start with '0', we do did to prevent collecting hidden files or so.
        tmpJsonSource = dir(fullfile(subjectDir,['0' extensions{2}])); % json

        % Extract filenames without extensions from tmpImgsource and tmpJsonSource
        imgFilenames  = cellfun(@(x) strtok(x, '.'), {tmpImgSource.name}, 'UniformOutput', false);
        jsonFilenames = cellfun(@(x) strtok(x, '.'), {tmpJsonSource.name}, 'UniformOutput', false);
        unmatchedFilenames = setdiff(imgFilenames, jsonFilenames); % assuming same order and not extra files.
        if runQuestions
            if ~isempty(unmatchedFilenames)
                error('Not every .nii file has a matching .json file or there is a mismatch');
            end
        end

        %% Loop through images in tmpsource and copy/rename files to rawdata
        for tmpScanNr = 1:length(tmpImgSource)
            tmpFile   = fullfile(tmpImgSource(tmpScanNr).folder,tmpImgSource(tmpScanNr).name);
            tmpJsFile = fullfile(tmpJsonSource(tmpScanNr).folder,tmpJsonSource(tmpScanNr).name);

            % ANATOMY: Find, copy and rename anatomical scan
            if contains(tmpFile,char(logs.ses(tmpsesID).tbl{logs.ses(tmpsesID).tbl{:,"Name"}=="anat","Description"}))
                destFile = fullfile(subjrawDir,anatDirs{tmpsesID}{1}, [sub(f).name '_' tmp_ses2run{tmpses} '_acq-' anatAcquisition{tmpsesID}{1}  '_' anatModalities{tmpsesID}{1} '.nii.gz']);
                copyfile(tmpFile,destFile);

                destJsFile = fullfile(subjrawDir,anatDirs{tmpsesID}{1}, [sub(f).name '_' tmp_ses2run{tmpses} '_acq-' anatAcquisition{tmpsesID}{1} '_' anatModalities{tmpsesID}{1} '.json']);
                copyfile(tmpJsFile,destJsFile);
            end

            % FUNCTIONAL: Find, copy and rename functional scans
            if contains(tmpFile,char(logs.ses(tmpsesID).tbl{logs.ses(tmpsesID).tbl{:,"Name"}=="func","Description"}))
                % If it is a functional run: get run number
                tmpRunNr = regexp(tmpImgSource(tmpScanNr).name, '\d+', 'match', 'once');
                tmpRunId = find(runs == str2double(tmpRunNr)); % compare run number to the ones we are interested in.

                if ~isempty(tmpRunId) % If the run number is one of the runs we are interested in: copy and rename niftis and json
                    destFile   = fullfile(subjrawDir,funcDir, [sub(f).name '_' tmp_ses2run{tmpses} '_task-' taskNames{tmpRunId} '_run-' num2str(runs(tmpRunId),formatSpecRun) '_bold.nii.gz']);
                    copyfile(tmpFile,destFile);

                    destJsFile = fullfile(subjrawDir,funcDir, [sub(f).name '_' tmp_ses2run{tmpses} '_task-' taskNames{tmpRunId} '_run-' num2str(runs(tmpRunId),formatSpecRun) '_bold.json']);
                    copyfile(tmpJsFile,destJsFile);
                end
            end

            % FIELDMAP: Find, copy and rename fieldmap scans
            if contains(tmpFile,char(logs.ses(tmpsesID).tbl{logs.ses(tmpsesID).tbl{:,"Name"}=="fmap","Description"}))
                if contains(tmpImgSource(tmpScanNr).name,'_e1.nii.gz') % magnitude 1
                    destFile   = fullfile(subjrawDir,fmapDir, [sub(f).name '_' tmp_ses2run{tmpses} '_magnitude1.nii.gz']);
                    copyfile(tmpFile,destFile);

                    destJsFile = fullfile(subjrawDir,fmapDir, [sub(f).name '_' tmp_ses2run{tmpses} '_magnitude1.json']);
                    copyfile(tmpJsFile,destJsFile);
                elseif contains(tmpImgSource(tmpScanNr).name,'_e2.nii.gz') % magnitude 2
                    destFile   = fullfile(subjrawDir,fmapDir, [sub(f).name '_' tmp_ses2run{tmpses} '_magnitude2.nii.gz']);
                    copyfile(tmpFile,destFile);

                    destJsFile = fullfile(subjrawDir,fmapDir, [sub(f).name '_' tmp_ses2run{tmpses} '_magnitude2.json']);
                    copyfile(tmpJsFile,destJsFile);
                elseif contains(tmpImgSource(tmpScanNr).name,'_e2_ph.nii.gz') % phase difference
                    destFile   = fullfile(subjrawDir,fmapDir, [sub(f).name '_' tmp_ses2run{tmpses} '_phasediff.nii.gz']);
                    copyfile(tmpFile,destFile);

                    destJsFile = fullfile(subjrawDir,fmapDir, [sub(f).name '_' tmp_ses2run{tmpses} '_phasediff.json']);
                    copyfile(tmpJsFile,destJsFile);
                else
                    error('Fmap file not recognized! Check the files!');
                end
            end
        end

    end % session is complete
    fprintf('....Subject: %s  is complete\n', sub(f).name); % display.
end %for-loop across subjects
end % end SortNiftisIntoRaw


