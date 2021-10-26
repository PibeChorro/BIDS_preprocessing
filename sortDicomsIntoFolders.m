function sortDicomsIntoFolders( prefix, log, params )
%% Sort Dicoms into Folders
% The following script sorts MRI output data (*.IMA/*.ima) into folders
% in "sourcedata" dir, while generating a folder for each sequence
% (exceptions: fieldmaps and localizer runs if named accordingly).
% The script sorts all dicoms into separate folders for each sequence.
% Afterwards it reads out one dicom for every sequence and checks its
% sequence description. If the description is unknown it asks for a proper
% name (like 'func', 'anat', 'localizer', etc.) and for the number of
% images per sequence.
%
% IMPORTANT:
%   The number of images must be correct. Every non-Anatomical folder containing more or
%   less images are stored in an 'error' folder.
%   Thus, it is recommended to create a "sequenceInfo" mat file BEFOREHAND with the structure
%   "log" containing the fields: sequenceDescriptions, sequenceNames,
%   sequenceScanNrs and sequenceModality (important for anatomical scans!)
%  
%    The location of the "sequenceInfo.mat" file CAN be provided with the
%    "params" structure. Alternatively one can provide the "log" structure
%    directly when calling this function.
%
%
% IMPORTANT: About fieldmaps and localizer runs!
%   Per convention please name fieldmaps sequences "fmap" and localizer sequences
%   "loc". These are NOT going to be separated into different
%   "runs-folders". Any sequence name CONTAINING either "loc" of "fmap"
%   will be put together! (NOT sequenceDescriptions as read from DICOM headers, BUT sequenceNames)
%
% Example of "log" structure in sequenceInfo:
%       log.sequenceDescriptions{2}: 'localizer 64-channel' --> from hdr
%       log.sequenceNames{2}: 'loc' --> name given.
%       log.sequenceScanNrs{2}: [13, 3] --> nr. of images.
%       log.sequenceModality{2}: 'T1w' --> modality of scan (e.g. 'BOLD',
%       'T1w', 'T2w', etc.)
%
%
% INPUT:
%   prefix (optional): string with prefix of your subjects
%   log (optional): structure with sequence information. ELSE: try to load,
%           or get sequence info from user.
%   params (optional): structure with fields
%           rootDir: root directory (project directory)
%           sourceDir: sourcedata dir (where the dicoms are. Usually: "/sourcedata" for BIDS-conform format)
%           formatSpecRun: format specification for RUNS (default: %0.2i).
%                       E.g. run-01, -02, etc.
%           extentions: dicom extentions to look for (*.ima and *.IMA)
%
% EXAMPLE USAGE:
%   sortDicomsIntoFolders('subPrefix','s', log, params)
%
% SUGGESTION:
%   Run 1 subject, save the sequenceInfo file and then re-run with the rest
%   of the subjects!
%
% Original: JB 07/2019
% Modified: VP 03/2021
% Modified: PG 10/2021
%
% In case you want to "unsort" a Subject (e.g. for debugging this script),
% here an useful snippet:
% files = dir(pwd);
% for tmpdir = 1:length(files)
%    try
%        movefile(string(fullfile(files(tmpdir).name,'*.IMA')),files(tmpdir).folder);
%    catch
%    end
% end

%% SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 1)
    clc;
    warning('No settings given!')
    help sortDicomsIntoFolders
    return
end


fprintf('\n====SORTING DICOMS INTO FOLDERS====\n');

%% Check the input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% What step to do?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(params, 'steps2sort'); steps2sort = params.steps2sort; 
    disp(join(['Running steps: ' string(params.steps2sort)])); else
    fprintf('====STEPS: Missing. Using default: true true true\n\n');
    steps2sort = [true true true];
end



% Get Directories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get/enter rootDir
if isfield(params, 'rootDir'); rootDir = params.rootDir; else
    fprintf('====ROOTDIR: Missing. No default. Please select your project folder \n (ideally it should contain a folder named "sourcedata")\n\n')
    rootDir = uigetdir(pwd, 'Select Project Folder [rootDir]');
end
if rootDir == 0; error('No folder was selected. Re-start'); end

% Get/enter sourceDir
if isfield(params, 'sourceDir'); sourceDir = params.sourceDir; else
    fprintf('====SOURCEDIR: Missing. No default. Please specify the directory that contains your subjectfolders (ideally it should be named "sourcedata")\n\n')
    sourceDir = uigetdir(rootDir, 'Select your sourcedata');
end
if sourceDir == 0; error('Error: no folder was selected'); end

% Make modalities dirs?
if isfield(params, 'mkModalityDirs'); mkModalityDirs = params.mkModalityDirs; else
    fprintf('====Make Modality Dirs: Missing. Using no modality dirs (default).')
    mkModalityDirs = false;
end

% Structure containing information about sequences (create/load)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if LOG is provided and if it has the necessary fields.
% ELSE: make a list for all series descriptions (e.g. myFuncSequence, MB_MPI_dwiSequence)
% and save it in a .mat file (called 'sequenceInfo'). In case you already used this script on a
% subset of your data try to load the saved .mat file

if ~isfield(log, 'sequenceDescriptions') || ~isfield(log, 'sequenceNames') || ~isfield(log, 'sequenceScanNrs') || ~isfield(log, 'sequenceModality')
    
    % Check if sequenceInfo name is provided, else use default name.
    if isfield(params, 'dir2sequenceInfo'); sequenceInfoName = params.dir2sequenceInfo; else
        fprintf('====No dir to sequenceInfo provided. Using default: sequenceInfo.mat')
        sequenceInfoName = 'sequenceInfo.mat';
    end
    
    % Check if there already exists a sequence description mat file
    if exist([sequenceInfoName],'file')
        disp('Loading sequenceInfo.mat from directory');
        load(sequenceInfoName,'log');
    else
        disp('No sequenceInfo! Values will need to be entered manually');
        log.sequenceDescriptions    = {}; % e.g. "AFNI_..."
        log.sequenceNames           = {}; % e.g. "anat"
        log.sequenceScanNrs         = {}; % e.g. "192"
        log.sequenceModality        = {}; % e.g. "T1w" (ONLY IMPORTANT FOR ANATOMICAL SCANS)
    end
else
    disp('SequenceInfo provided when running function');
end

% Get subjects in dir
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

if isempty(sub)
    error('The specified directory does not contain any folders starting with the specified prefix');
end

% Sanity check: check if folder and no hidden folder
areFolders = all([sub(:).isdir]);
areHidden  = any(contains({sub(:).name},'.'));

if ~areFolders % check if files are all folders
    error('Your prefix seems not to be unique for folders')
elseif areHidden % check if the prefix was correctly set and no hidden files are selected
    error('It looks like you did not select a valid prefix, because a hidden folder was selected.')
end

% Extra specficiations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format specification for RUNS
if isfield(params, 'formatSpecRun'); formatSpecRun = params.formatSpecRun; else; formatSpecRun  = '%02i'; % specify how you name your RUNS
    fprintf(['====FORMATSPECRUN: No formatSpect specified. Using default formatSpecRun for RUNS %02i. \n\n']);  %#ok<*NBRAK>
end

% extension to convert
if isfield(params, 'extensions'); extensions = params.extensions; else
    extensions = {'**.IMA','**.ima'}; % extension you care about
    fprintf(['====EXTENSIONS: No extensions specified. Using default extensions: .IMA and .ima \n\n']);
end


%% SORT DICOMS IN FOLDERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f = 1:length(sub) % Outer loop for the subjects
    fprintf('Sorting subject: %s \n', sub(f).name); % display.
    subjectDir = fullfile (sourceDir,sub(f).name);
    
    %% STEP 1: Move images from separate sequences in different folders
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fill an array with all relevent data files
    
    if steps2sort(1) == true
        % Get all images
        dicoms = [];
        for ext = 1:length(extensions)
            dicoms = [dicoms; spm_select('FPList', subjectDir, extensions{ext})]; % this looks for the *IMA or *ima, if no image is found, then this step was either performed or no images are present.
        end
        
        % Sanity check: check if the dicoms array is empty.
        % If so it could mean that this sorting procedure already has been
        % performed on the folder or it is empty
        if isempty(dicoms)
            folderContent = dir(subjectDir);
            if all(contains({folderContent(:).name}, '.'))
                warning('Folder %s does not contain any image \n', subjectDir)
                continue; % next iteration (subject)
            elseif sum(isfolder(arrayfun(@(x) fullfile(sourceDir, sub(f).name, x.name), folderContent, 'UniformOutput', false))) > 3 % Check if you have any folders (apart from '.' and '..')
                warning(['It seems this folder was sorted before in any way (it contains sub-folders)\n'...
                    'You better go check this one out before you proceed.\n\n']);
                continue; % next iteration (subject)
            end
        else
            fprintf('\n Step 1:....Sorting dicoms into folders (01,02,03, ...). \n Total number of files to sort: %s \n', num2str(length(dicoms))); % display.
            % Now look for different runs and make a new folder for each one
            for tmpFile  = 1:length(dicoms) % go through all relevant files
                fileName = regexp(dicoms(tmpFile,:),'\.','split'); % returns list of splitted parts of filename:
                seqNum   = str2double(fileName{4}); % 4th part = sequence number
                if ~isfolder(fullfile(subjectDir, num2str(seqNum,formatSpecRun))) % if not existing, make folder
                    mkdir(fullfile(subjectDir, num2str(seqNum,formatSpecRun)))
                end
                movefile(string(dicoms(tmpFile,:)),...
                    fullfile(subjectDir,num2str(seqNum,formatSpecRun))); % move file into sequence folder
            end
        end
        
    else
        fprintf('\n Skipping step 1 \n'); % display.
    end
    
    %% STEP 2: Move the folders containing images into the corresponding sequence folders
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % At this point your subject folder should look like this:
    % 01
    % 02
    % ...
    % each subfolder contains the DICOMs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if steps2sort(2) == true
        subDirs = dir(subjectDir); % get the folders you just made within one subject
        fprintf('\n Step 2:....Renaming folders to "sequenceNames". Total number of folders: %s \n', num2str(length(subDirs))); % display.
        % Iterate over all folders and if it is a numeric (from the
        % sequencing before, check the seriesDescription in the header of
        % the first dicom and move to the corresponding sequenceName
        for sd = 1:length(subDirs)
            [status,~] = str2num(subDirs(sd).name); % check if the folder is not any kind of hidden folder but one of the previously created (--> looks for number! e.g., 01, 02, etc.)
            if status
                currentDir = fullfile(subjectDir,subDirs(sd).name);
                % read in the first dicom header
                dicoms = [];
                for ext = 1:length(extensions)
                    dicoms = [dicoms; spm_select('FPList', currentDir, extensions{ext})]; %#ok<*AGROW>
                end
                
                firstDicom        = spm_dicom_headers(dicoms(1,:));  % get header
                seriesDescription = firstDicom{1}.SeriesDescription; % read out the 'SeriesDescription' field
                sequenceIndex     = find (strcmp(log.sequenceDescriptions, seriesDescription));
                
                % If there is a new sequence, we add it to our list
                if isempty(sequenceIndex)
                    fprintf (['A new sequence description was found.\n'...
                        'Please assign the correct data type to the new series description.\n'...
                        'Use meaningful names, ideally in line with BIDS naming,\n'...
                        'e.g. anat, func, dwi, etc.']);
                    log.sequenceDescriptions{end+1} = seriesDescription;
                    log.sequenceNames{end+1}        = input (['\n' log.sequenceDescriptions{end} ': '],'s');
                    imageNumber                     = input(['\n\nPlease asign a number of scans to this sequence.\n\n'...
                        '!!!!!!IMPORTANT!!!!!\n\n'...
                        'Make sure you asign the correct number of scans (for EPI). Every folder containing more or less scans will be stored as "error-run-XX"\n\n'...
                        'If more than one number of scans are possible for a specific sequence, separate the numbers by a SPACE.\n'],'s');
                    log.sequenceScanNrs{end+1}      = str2num(imageNumber); %#ok<*ST2NM>
                    sequenceIndex                   = length(log.sequenceNames);
                    log.sequenceModality{end+1}     = input(['\n\n Please assign the modality (i.e. BOLD, T1w, T2w or T2star): '],'s');
                end
                
                if ~isfolder(fullfile(subjectDir,log.sequenceNames{sequenceIndex}))
                    mkdir (fullfile(subjectDir,log.sequenceNames{sequenceIndex}));
                end
                
                % Check for anatomical scans: built in heuristic that if an "anat" or "struct" folder exists just rename the folder to anat
                % NOT controlling for the number of images.
                if (contains(log.sequenceNames{sequenceIndex},'anat')||contains(log.sequenceNames{sequenceIndex},'struct')) % if anatomical
                    if mkModalityDirs % make modality specific dirs?
                        movefile(string(fullfile(currentDir,'*')), fullfile(subjectDir,'anat', log.sequenceModality{sequenceIndex}));
                        rmdir(fullfile(currentDir));
                    else
                        movefile(string(fullfile(currentDir,'*')), fullfile(subjectDir,'anat'));
                        rmdir(fullfile(currentDir));
                    end
                else % move files to the corresponding folders
                    movefile(string(fullfile(currentDir)),...
                        fullfile(subjectDir,log.sequenceNames{sequenceIndex}));
                end
            end
        end
    else
        fprintf('\n Skipping step 2 \n'); % display.
    end
    
    
    %% STEP 3: Go into the sequence folders and rename the containing folders to run-<runNr>
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iterate over the sequence namefolders and rename the contained
    % folders to run-<runNr>.
    % TO-DO: wenn die falsche Anzahl an scans in der folder gefunden werden dann gibts ein Problem.
    if steps2sort(3) == true
        fprintf('\n Step 3:....Changing sub-folder names to run-<runNr> or collapsing folders \n'); % display.
        for i = 1:length(log.sequenceNames) % go throw the different sequence types.
            % get all folders within the sequencefolder
            subDirs    = dir(fullfile(subjectDir,log.sequenceNames{i}));
            runCounter = 1; % set a counter that determines the run number
            
            if contains(log.sequenceNames{i},"fmap") || contains(log.sequenceNames{i},"loc") % fieldmaps and localizers are thrown together
                together = true;
            else
                together = false;
            end
            
            for sd = 1:length(subDirs)
                [status,~]     = str2num(subDirs(sd).name); % Here we avoid hidden folders like '..'
                if status
                    currentDir = fullfile(subjectDir, log.sequenceNames{i}, subDirs(sd).name);
                    
                    % Count all dicoms in the folder
                    dicoms = [];
                    for ext = 1:length(extensions)
                        dicoms = [dicoms; spm_select('FPList', currentDir, extensions{ext})];
                    end
                    nrOfScans = length(dicoms(:,1));
                    
                    % only if the number of dicoms in the folder is equal to
                    % the number asigned above, the folder accepted as a 'run'
                    % folder
                    if together % put files into the sequenceNames folder; ignore the wanted number of scans
                        movefile(string(fullfile(currentDir, '*.IMA')), fullfile(subjectDir,log.sequenceNames{i}));
                        rmdir(fullfile(currentDir));
                    else
                        if (ismember(nrOfScans, log.sequenceScanNrs{i}))
                            movefile(string(currentDir), fullfile(subjectDir,log.sequenceNames{i},['run-' num2str(runCounter,formatSpecRun)]))
                            runCounter = runCounter+1;
                        else
                            movefile(string(currentDir), fullfile(subjectDir,log.sequenceNames{i},['error-run-' num2str(runCounter,formatSpecRun)]))
                            runCounter = runCounter+1;
                        end
                    end
                end
            end
        end % rename to run-<runNr>
    else
        fprintf('\n Skipping step 3 \n'); %
    end
    fprintf('....Subject: %s  is complete\n', sub(f).name); % display.
end %for-loop across subjects

fprintf('\n Sorting Dicoms into folders complete! \n'); %
% Save the log in case you want to rerun the script on more subjects
% but don't want to reenter all the sequence information
save(sequenceInfoName,'log')

end
