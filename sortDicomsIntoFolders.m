%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% The following script sorts MRI output data (*.IMA/*.ima) into folders
% while generating a folder for each sequence.
% The script sorts all dicoms into seperate folders for each sequence.
% Afterwards it reads out one dicom for every sequence and checks its
% sequence description. If the descriptio is unknown it asks for a proper
% name (like 'func', 'anat', 'localizer', etc.) and for the number of
% images per sequence
% !!!IMPORTANT!!! The number of images must be correct. Every folder
% containing more or less images are stored in an 'error' folder.
%
% author: Jasper Bischofberger <jasper.bischofberger@posteo.de>
% modified: Vincent Plikat <vincent.plikat@student.uni-tuebingen.de>
%
% USAGE:
% sortDicomsIntoFolders
%
% Original: JB 07/2019
% Modified: VP 03/2021
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

clear all; close all;

%% Define important details of your file structure and location
workingDir = pwd;
%-------------------------------------------------------------------------%
% DEFINE path and get relevant filenames.
%-------------------------------------------------------------------------%
fprintf(['Please specify the directory that contains your subjectfolders'...
    '(ideally it should be named "sourcedata")\n\n']);

sourceDir = uigetdir(homedir, 'Select your sourcedata');

if sourceDir == 0
    error('No folder was selected --> I terminate the script')
end

%% the structure that contains the information about sequences
% make a list for all series descriptions (e.g. myFuncSequence, MB_MPI_dwiSequence)
% and save it in a .mat file. In case you already used this script on a
% subset of your data try to load the saved .mat file

sequenceInfoName = 'sequenceInfo';

% check if there already exists a sequence description mat file
if exist([sequenceInfoName '.mat'],'file')
    load(sequenceInfoName)
else
    log.sequenceDescriptions    = {};
    log.sequenceNames           = {};
    log.sequenceScanNrs         = {};
end


% specify format for folder numeration 
formatSpec = '%02i';

%-------------------------------------------------------------------------%
% DEFINE prefix of subject folders. Make sure that the prefix is unique for
% your subject folders and you do not have any other files starting with
% this prefix.
%-------------------------------------------------------------------------%

prefix = input (['Please specify the prefix of your participant data.\n' ...
    '(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');

fprintf('\n\n');

sub = dir(fullfile(sourceDir,[prefix '*'])); % This will give you also the number of subjects = number of folders.

if isempty(sub)
    error('The specified directory does not contain any folders starting with the specified prefix');
end

%-------------------------------------------------------------------------%
% DEFINE file extensions you want sort
extensions = {'**.IMA','**.ima'}; % extension you care about
% IMPORTANT: It seems like dir() (at least on MacOS) is not case sensitive, but 
% spm_select('FPList',...) is, so we use spm_select
%-------------------------------------------------------------------------%

%% check if the selected files are only folders and that there are no hidden folders selected

areFolders = all([sub(:).isdir]);
areHidden  = any(contains({sub(:).name},'.'));
% check if files are all folders
if ~areFolders
    error('Your prefix seems not to be unique for folders')
% check if the prefix was correctly set and no hidden files are selected
elseif areHidden
    error('It looks like you did not select a prefix, because a hidden folder was selected.')
end
%% The actual sorting of your files

% Outer loop for the subjects
for f = 1:length(sub)
    subjectDir = fullfile (sourceDir,sub(f).name);
    %% Move images from seperate sequences in different folders
    % fill an array with all relevent data files
    dicoms = [];
    for ext = 1:length(extensions)
        dicoms = [dicoms; spm_select('FPList', subjectDir, extensions{ext})];
    end
    
    % check if the dicoms array is empty. 
    % If so it could mean that this sorting procedure already has been
    % performed on the folder or it is empty
    if isempty(dicoms)
        folderContent = dir(subjectDir);
        if all(contains({folderContent(:).name}, '.'))
            warning ('Folder %s seems to be empty\n', subjectDir)
            continue;
        else
            warning (['It seems this folder was sorted before in any way\n'...
                'You better go check this one out before you proceed.\n\n'])
        end
    else
        %-----------------------------------------------------------------%
        % Look for different sequences and establish a folder each. 
        %-----------------------------------------------------------------%
        for tmpFile = 1:length(dicoms) % go through all relevant files
            fileName = regexp(dicoms(tmpFile,:),'\.','split'); % returns list of splitted parts of filename: 
            seqNum = str2double(fileName{4}); % 4th part = sequence number
            if ~isfolder(fullfile(subjectDir, num2str(seqNum,formatSpec))) % if not existing, make folder
                mkdir(fullfile(subjectDir, num2str(seqNum,formatSpec)))
            end
            status = movefile(string(dicoms(tmpFile,:)),...
                fullfile(subjectDir,num2str(seqNum,formatSpec))); % move file into sequence folder
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % At this point your subject folder should look like this:
    % 01 
    % 02
    % ...
    % each subfolder contains the DICOMs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subDirs = dir(subjectDir);     % get the folders you just made within one subject
    
    %% Move the folders containing images into the corresponding sequence folders
    % Iterate over all folders and if it is a numeric (from the
    % sequencing before, check the seriesDescription in the header of
    % the first dicom and move to the corresponding sequenceName

    for sd = 1:length(subDirs)
        [status,num] = str2num(subDirs(sd).name);       % check if the folder is not any kind of hidden folder but one of the previously created (--> looks for number! e.g., 01, 02, etc.)
        if status
            currentDir = fullfile(subjectDir,subDirs(sd).name);
            % read in the first dicom header
            dicoms = [];
            for ext = 1:length(extensions)
                dicoms = [dicoms; spm_select('FPList', currentDir, extensions{ext})];
            end
            firstDicom = spm_dicom_headers(dicoms(1,:));
            % read out the 'SeriesDescription' field
            seriesDescription = firstDicom{1}.SeriesDescription;
            sequenceIndex = find (strcmp(log.sequenceDescriptions,seriesDescription));
            
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
                    'Make sure you asign the correct number of scans. Every folder containing more or less scans will be stored as "error-run-XX"\n\n'...
                    'If more than one number of scans are possible for a specific sequence, seperate the numbers by a SPACE.\n'],'s');
                log.sequenceScanNrs{end+1}      = str2num(imageNumber);
                sequenceIndex                   = length(log.sequenceNames);
            end

            if ~isfolder(fullfile(subjectDir,log.sequenceNames{sequenceIndex}))
                mkdir (fullfile(subjectDir,log.sequenceNames{sequenceIndex}));
            end
            
            % Built in heuristic that if an "anat" or "struct" folder exists,
            % just rename the folder to anat
            if (contains(log.sequenceNames{sequenceIndex},'anat')||contains(log.sequenceNames{sequenceIndex},'struct'))
                nrOfScans = length(dicoms(:,1));
                anatomicalModality = input('You assigned the anatomical sequence. Please assign the modality (i.e. T1w, T2w or T2star): ','s');
                % check if the number of images is right
                if ismember(nrOfScans, log.sequenceScanNrs{sequenceIndex})
                    movefile(string(fullfile(currentDir,'*')), fullfile(subjectDir,'anat',anatomicalModality))
                else
                    movefile(string(fullfile(currentDir,'*')), fullfile(subjectDir,'anat',['error-' anatomicalModality]))
                end
                rmdir(fullfile(currentDir))
            else
                movefile(string(fullfile(currentDir)),...
                    fullfile(subjectDir,log.sequenceNames{sequenceIndex}));
            end
        end
    end

    %% Go into the sequence folders and rename the containing folders 
    % Iterate over the sequence namefolders and rename the contained
    % folders to run-<runNr>. 
    for i = 1:length(log.sequenceNames)
        % get all folders within the sequencefolder
        subDirs = dir(fullfile(subjectDir,log.sequenceNames{i}));
        runCounter = 1;     % set a counter that determines the run number
        for sd = 1:length(subDirs)
            [status,num] = str2num(subDirs(sd).name);   % here we avoid hidden folders like '..'
            if status
                currentDir = fullfile(subjectDir, log.sequenceNames{i}, subDirs(sd).name);
                
                %count all dicoms in the folder
                dicoms = [];
                for ext = 1:length(extensions)
                    dicoms = [dicoms; spm_select('FPList', currentDir, extensions{ext})];
                end
                nrOfScans = length(dicoms(:,1));
                
                % only if the number of dicoms in the folder is equal to
                % the number asigned above, the folder accepted as a 'run'
                % folder
                if (ismember(nrOfScans, log.sequenceScanNrs{i}))
                    movefile(string(currentDir), fullfile(subjectDir,log.sequenceNames{i},['run-' num2str(runCounter,formatSpec)]))
                    runCounter = runCounter+1;
                else
                    movefile(string(currentDir), fullfile(subjectDir,log.sequenceNames{i},['error-run-' num2str(runCounter,formatSpec)]))
                end
            end
        end
    end
    % End of subject
end
% Safe the log in case you want to rerun the script on more subjects
% but don't want to reenter all the sequence information
save(sequenceInfoName,'log')
