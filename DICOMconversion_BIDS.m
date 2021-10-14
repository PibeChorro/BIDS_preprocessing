function DICOMconversion_BIDS( varargin )
%% DICOM to NIFTI conversion using spm12
%  This script will convert sorted DICOM folders in "sourcedata" into a
%  BIDS-conform "rawdata" with data in NIFTI-format. The conversion can
%  either be done with SPM OR dcm2niix (default and suggested).
%
%  INPUTS:
%       The script need specific parameters, if these are not provided, these
%       will either be asked for or the default values are being use.
%   
%       rawSubNames (required): cell containing the wanted subjects name (for rawdata folder, e.g. "sub-001")
%       subNames (required): cell containing the subjects you want to
%               convert FROM (from sourcedata folder, e.g. "PW001")
%       softwareFlag (optional): 'dcm2niix' OR 'SPM' (dcm2niix)
%       funcConversion (optional): 'true' or 'false' (true)
%       anatConversion (optional): 'true' or 'false' (true)
%       fieldmapConversion (optional): 'true' or 'false' (true)
%       params (optional): parameter-defining structure with fields:
%           rootDir (optional) ELSE get dir using "uigetdir"
%           sourceDir and rawDir (optional) ELSE use defaults 'sourcedata'
%                    and 'rawdata'
%           funcDir (optional) ELSE use default: func
%           anatDir (optional) ELSE use default: anat
%           fmapDir (optional) ELSE use default: fmap
%           anatModalities, anatAcquisition, anatModalitiesDir (optional)
%                    ELSE use defaults: T1w, T2w, T2star, with 1 mm
%           formatSpec (optional) ELSE use default: '%.02i' (for runs, e.g., run-01)
%           extensions (optional) ELSE use default: 'ima' and 'IMA'
%
%       SOFTWARE SPECIFIC PARAMETERS:
%           FOR dcm2niix:
%               path2exe (optional) ELSE get with 'uigetdir' OR press 0 to
%                   use defaults (depending on operation system)
%           FOR SPM:
%               nDummies (optional - no default) ELSE ask for number of
%               dummies (in func scans) and ask for confirmation.
%
%  OUTPUT: rawdata folder(s) in BIDS-format
%         sub-<SUBJ_NR> (e.g., sub-001)
%                --> func
%                       --> sub-001-task_<TASKNAME>_run-<RUN_NR>_bold.json
%                       --> sub-001-task_<TASKNAME>_run-<RUN_NR>_bold.nii.gz
%                       ...
%                --> anat
%                       --> sub-001_acq->ACQUISITION_LABEL>_<MODALITY>.json
%                       --> sub-001_acq->ACQUISITION_LABEL>_<MODALITY>.gz
%                       ...
%                --> fmap
%                       --> sub-001_<FIELDMAP_LABEL>.json
%                       --> sub-001_<FIELDMAP_LABEL>.nii.gz
%                ...
%        ...
%                            
%
%  IMPORTANT: 
%            NO DUMMIES ARE EXTRACTED WHEN USING DCM2NIIX! 
%            DUMMIES ARE EXTRACTED WHEN USING SPM!   
%
% VP 10/2020 (adapted from JBs preprocessingSPM12)
% PG 10/2021 (modified, included parsing, dcm2niix, etc)
%
% TO-DO List
%         * create a log file to save warning and error messages
%         * implement flag "overwrite", delete old files and do new
%         conversion?

%% SETTINGS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 1)
    clc;
    warning('No settings given!')
    help DICOMconversion_BIDS
    return
end

%% Parse the input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addRequired(p, 'rawSubNames', @(x) iscell(x)); % cell containing the subjects you want to convert.
addRequired(p, 'subNames', @(x) iscell(x)); % cell containing the subjects you want to convert.
addParameter(p, 'softwareFlag', 'dcm2niix', @(x) ischar(x)); % Select: 'SPM' or 'dcm2niix' % Dicom transformation using SPM or dcm2niix
addParameter(p, 'funcConversion', true, @(x) islogical(x)); % do functional conversions
addParameter(p, 'anatConversion', true, @(x) islogical(x)); % do anatomical conversions
addParameter(p, 'fieldmapConversion', true, @(x) islogical(x)); % do fieldmap conversions
addParameter(p, 'params', [], @(x) isstruct(x)); % Check if you have a settings input
parse(p, varargin{:});
do = p.Results;

%% Get SubNames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawSubNames = do.rawSubNames;
subNames    = do.subNames;

%% Get values or set defaults from params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = do.params;

% Get/enter rootDir
if isfield(params, 'rootDir'); rootDir = params.rootDir; else
    fprintf(['====ROOTDIR: No default. Please select your project folder \n (ideally it should contain a folder named "sourcedata")\n\n'])
    rootDir = uigetdir(pwd, 'Select Project Folder [rootDir]');
end

% Check rootDir
if rootDir == 0; error('No folder was selected. Re-start'); end

% source and raw dirs
if isfield(params, 'sourceDir') && isfield(params, 'rawDir')
    fprintf(['====SOURCE AND RAW DIRS: No dir specified. Using defaults name: sourcedata and rawdata.\n\n']);
    sourceDir = params.sourceDir;
    rawDir    = params.rawDir;
else
    sourceDir = fullfile(rootDir, 'sourcedata');  % sourcedata
    rawDir    = fullfile(rootDir, 'rawdata');     % rawdata
end

if isfield(params, 'funcDir'); funcDir = params.funcDir; else; funcDir = 'func';
    fprintf(['====FUNC DIR: No dir specified. Using default name: func.\n\n']); 
end %#ok<*NBRAK>

if isfield(params, 'anatDir'); anatDir = params.anatDir; else; anatDir = 'anat';
    fprintf(['====ANAT DIR: No dir specified. Using default name: anat\n\n']); 
end

if isfield(params, 'fmapDir'); fmapDir = params.fmapDir; else; fmapDir = 'fmap';
    fprintf(['====FMAP DIR: No dir specified. Using default name: fmap\n\n']); 
end

% anatomical modalities to use
if isfield(params, 'anatModalities') && isfield(params, 'anatAcquisition') && isfield(params, 'anatModalitiesDirs') % specify the anatomical modalities to use and acquisition values
    anatModalities      = params.anatModalities;
    anatModalitiesDirs  = params.anatModalitiesDirs;
    anatAcquisition     = params.anatAcquisition;
else
    anatModalities      = {'T1w','T2w','T2star'};
    anatAcquisition     = {'1mm','1mm','1mm'};
    anatModalitiesDirs  = {'T1w','T2w','T2star'};
    fprintf(['====ANATMODALITIES & ACQUISITION & DIRS: No anatModalities/acquisition/dirs specified. Using default anatModalities: T1wm T2w and T2star and 1mm.\n\n']);
end

% format specification for RUNS
if isfield(params, 'formatSpec'); formatSpec = params.formatSpec; else; formatSpec  = '%02i'; % specify how you name your RUNS
    fprintf(['====FORMATSPEC: No formatSpect specified. Using default formatSpec for RUNS %02i. \n\n']); 
end

% extension to convert
if isfield(params, 'extensions'); extensions = params.extensions; else
    extensions = {'**.IMA','**.ima'}; % extension you care about
    fprintf(['====EXTENSIONS: No extensions specified. Using default extensions: .IMA and .ima \n\n']);
end

%% Software specific settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dcm2nixx: path 2 exe
if strcmp(do.softwareFlag,'dcm2niix')
    if isfield(params, 'path2exe')
        path2exe = params.path2exe;
    else
        fprintf(['====PATH2EXE: Enter path to dcm2niix or cancel for default values. \n\n'])
        path2exe = uigetdir(pwd, 'Select path to dcm2niix or press 0 (to use defaults)');
        if path2exe == 0
            warning('\n\nUsing default paths for you operating system. Probably wrong!\n\n');
            if ispc
                error('No dcm2niix default path implemented for windows yet'); 
            elseif ismac
                path2exe = '/Applications/MRIcroGL.app/Contents/Resources/'; % Path to dcm2niix exe in Mac OS
                disp(['Using path: ' path2exe]);
            else
                error('No dcm2niix default path implemented for linux yet');
            end
        end
    end
end

% SPM: nDummies.
if do.funcConversion == 1
    wait4confirmation = true;
    if isfield(params, 'nDummies')
        nDummies = params.nDummies;
        disp(['====nDUMMIES. If running SPM, using specified values of: ' num2str(nDummies)]);
        while wait4confirmation
            tmp = input('IS THIS CORRECT? (y/n) [n]:','s');
            if strcmp(tmp,'y')
                wait4confirmation = false;
            end
        end
    else
        while wait4confirmation
            nDummies = input(['\nPlease specify the number of dummy images (important for SPM convertion).\n'...
                '!!!IMPORTANT!!!\n'...
                'Please give the right number of images, because those will be skipped in the conversion process!\n\n']);
            disp(['Value entered: ' num2str(nDummies)]);
            tmp = input('IS THIS CORRECT? (y/n) [n]:','s');
            if strcmp(tmp,'y')
                wait4confirmation = false;
            end
        end
    end
end



%% Start conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ss = 1:length(rawSubNames) % For all subjects do each ...
    
    sourceSubDir       = fullfile(sourceDir,subNames{ss}); % dicoms
    sourceSubFuncDir   = fullfile(sourceSubDir,funcDir); % functional dicoms Dirs
    rawSubDir          = fullfile(rawDir,rawSubNames{ss}); % niftis
    rawSubFuncDir      = fullfile(rawSubDir,funcDir); % functional niftis Dirs
    
    %% STEP 1: Conversion from DICOM to NIfTI of ANATOMICAL image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do.anatConversion
        for mod = 1:length(anatModalities)
            try
                % Without any specific folder for anatomical modalities in rawdata:
                currDir = fullfile(sourceSubDir,anatDir, anatModalitiesDirs{mod}); % sourcedata dir for anat modality
                destDir = fullfile(rawSubDir, anatDir);   % rawdata dir for anat modality
                
                % Uncomment for specific subfolder for each anatomical modality in rawdata:
                % currDir = fullfile(sourceSubDir,anatDir, anatModalities{mod}); % sourcedata dir for anat modality
                % destDir = fullfile(rawSubDir, anatDir, anatModalities{mod});   % rawdata dir for anat modality
                
                folderContent = dir(fullfile(currDir,'run*')); % check if you have different runs, these will be put together in the same modality folder
                if isempty(folderContent) % no runs in this modality, set folderContent to 1
                    folderContent = 1;
                end
                
                for i = 1:length(folderContent) % loop through all folders found in this modality, if no run is there the just do this once.
                    if isnumeric(folderContent)
                        tmpDir = currDir;
                    else
                        tmpDir = fullfile(currDir, folderContent(i).name);
                    end
                    dirfiles = [];
                    for ext = 1:length(extensions)
                        dirfiles = [dirfiles; spm_select('FPList', tmpDir, extensions{ext})];
                    end
                    
                    if isempty(dirfiles)
                        error('FOLDER CONTENT IS EMPTY - PROBABLY WRONG PATH.')
                    else
                        fprintf('=> importing anatomical image to %s\n', destDir);
                        if strcmp(do.softwareFlag,'SPM')
                            % specify spm options
                            matlabbatch{1}.spm.util.import.dicom.data               = cellstr(dirfiles);
                            matlabbatch{1}.spm.util.import.dicom.root               = 'flat';
                            matlabbatch{1}.spm.util.import.dicom.outdir             = cellstr(destDir);
                            matlabbatch{1}.spm.util.import.dicom.protfilter         = '.*';
                            matlabbatch{1}.spm.util.import.dicom.convopts.format    = 'nii';
                            matlabbatch{1}.spm.util.import.dicom.convopts.meta      = 0;
                            matlabbatch{1}.spm.util.import.dicom.convopts.icedims   = 0;
                            spm_jobman('run', matlabbatch);
                            
                        elseif strcmp(do.softwareFlag,'dcm2niix')
                            if isnumeric(folderContent) % in case there was only one anat image in this modality
                                str4dcm2niix = sprintf([fullfile(path2exe, 'dcm2niix'), ' -b y -v 0 -z y -o', ' %s', ' -f', ' %s', ' %s'],destDir, ['/' rawSubNames{ss} '_acq-' anatAcquisition{mod} '_' anatModalities{mod}], tmpDir);
                                unix(str4dcm2niix);
                            else % in case there are more than one run in this modality
                                str4dcm2niix = sprintf([fullfile(path2exe, 'dcm2niix'), ' -b y -v 0 -z y -o', ' %s', ' -f', ' %s', ' %s'],destDir, ['/' rawSubNames{ss} '_acq-' anatAcquisition{mod} '_run-' num2str(i,formatSpec) '_' anatModalities{mod}], tmpDir);
                                unix(str4dcm2niix);
                            end
                        end % software
                    end  % if there are files
                end % for folderContent
                
                if strcmp(do.softwareFlag,'SPM') % change the names in case you use SPM
                    if folderContent == 1
                        anatImg     = spm_select('FPList',destDir,'.nii');
                        [stat, mes] = movefile(anatImg, fullfile(destDir,[rawSubNames{ss} '_acq-' anatAcquisition{mod} '_' anatModalities{mod} '.nii']));
                    else
                        anatImg = spm_select('FPList',destDir,'.nii');
                        for j = 1:length(anatImg)
                            [stat, mes] = movefile(anatImg(j,:), fullfile(destDir,[rawSubNames{ss} '_acq-' anatAcquisition{mod} '_run-' num2str(j,formatSpec) '_' anatModalities{mod} '.nii']));
                        end
                    end
                    if ~stat
                        warning(mes)
                    end
                end
                
            catch ME
                warning('Anatomical dicom in %s could not be converted.\n',tmpDir)
                fprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                    ME.stack(1).name, ME.stack(1).line, ME.message);
            end
            
            if strcmp(do.softwareFlag,'SPM')
                % after creating the anatomical NIfTIs we create its
                % corresponding .json file. Note we manually create only
                % ONE json file per file per modality.
                BIDS_anatT1w_json(destDir, dirfiles(1,:),[rawSubNames{ss} '_acq-' anatAcquisition{mod} '_' anatModalities{mod} '.nii.gz']);
            end
        end % for anatomical modalities
    end % for anatomical conversion
    
    
    
    %% STEP 2: Conversion from functional DICOM to NIfTI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do.funcConversion
        
        %% Read TSV file, get task names, runs and number of images
        if isfile(fullfile(rawSubDir,[rawSubNames{ss} '_scans.tsv']))
            movefile(fullfile(rawSubDir,[rawSubNames{ss} '_scans.tsv']),fullfile(rawSubDir,[rawSubNames{ss} '_scans.txt']),'f');
            [tmpTSV]= readtable(fullfile(rawSubDir,[rawSubNames{ss} '_scans.txt']));
            movefile(fullfile(rawSubDir,[rawSubNames{ss} '_scans.txt']),fullfile(rawSubDir,[rawSubNames{ss} '_scans.tsv']),'f');
        else
            error('No *_scans.tsv file found');
        end
        
        % Get task names:
        taskNames = tmpTSV.task;
        
        % Get runs
        runs = tmpTSV.run;
        
        folderContent   = dir(fullfile(sourceSubFuncDir,'run*'));
        
        if length(runs)~=length(folderContent) % sanity check
            error ('Number of runs found in *_scan.tsv file does not match the number of runs found in your subject folder')
        end
        
        %% Conversion Dicom 2 nifti
        % Check if the 'func' folder of your current subject is empty.
        if isempty(folderContent) % no dicoms in folder
            error('FOLDER CONTENT IS EMPTY - PROBABLY WRONG PATH.')
        else
            fprintf('FOLDER CONTENT FOUND \n')
            fprintf('========STARTING CONVERTION FROM DICOM TO NIFTI========\n\n');
            
            for i = 1:length(folderContent)            % loop through all folders found
                try
                    tmpTaskName = taskNames{i};
                    % get the directory of the current run in your 'sourcedata'
                    currDir = fullfile(sourceSubFuncDir, folderContent(i).name);
                    
                    % fill an array with all DICOM files
                    dirfiles = [];
                    for ext = 1:length(extensions)
                        dirfiles = [dirfiles; spm_select('FPList', currDir, extensions{ext})]; %#ok<*AGROW>
                    end
                    
                    if isempty(dirfiles)
                        error('FOLDER CONTENT IS EMPTY - PROBABLY WRONG PATH.')
                    end
                    
                    if strcmp(do.softwareFlag,'SPM')
                        niftinames = fullfile (rawSubFuncDir, [rawSubNames{ss} '_task-' tmpTaskName '_run-' num2str(i,formatSpec) '_bold.nii']);
                        warning("DUMMIES ARE DISCARDED WHEN USING SPM.")
                        % specify spm options
                        matlabbatch{1}.spm.util.import.dicom.data               = cellstr(dirfiles(nDummies+1:end,:));
                        matlabbatch{1}.spm.util.import.dicom.root               = 'flat';
                        matlabbatch{1}.spm.util.import.dicom.outdir             = cellstr(rawSubFuncDir);
                        matlabbatch{1}.spm.util.import.dicom.protfilter         = '.*';
                        matlabbatch{1}.spm.util.import.dicom.convopts.format    = 'nii';
                        matlabbatch{1}.spm.util.import.dicom.convopts.meta      = 0;
                        matlabbatch{1}.spm.util.import.dicom.convopts.icedims   = 0;
                        
                        fprintf('=> importing dicoms from %s\n', currDir);
                        fprintf('                      to %s\n', rawSubFuncDir);
                        
                        % start the actual job that converts DICOMs to NIfTIs
                        spm_jobman('run', matlabbatch);
                        clear matlabbatch
                        
                        % after conversion we select all created NIfTIs and
                        % convert them into one single 4D NIfTI
                        dicomFile                           = dirfiles(1,:);
                        hdr                                 = spm_dicom_headers(dicomFile);
                        currentNiftis                       = spm_select ('FPList', rawSubFuncDir, ['^f' '.*nii']);
                        
                        matlabbatch{1}.spm.util.cat.vols    = cellstr(currentNiftis);
                        matlabbatch{1}.spm.util.cat.name    = niftinames;
                        matlabbatch{1}.spm.util.cat.dtype   = 4;
                        
                        matlabbatch{1}.spm.util.cat.RT      = hdr{1}.RepetitionTime/1000;   % Repetition Time in seconds
                        spm_jobman('run', matlabbatch);
                        
                        clear matlabbatch
                        
                        % After conversion we delete the 3D images
                        delete(fullfile(rawSubFuncDir, 'f*.nii'))
                        
                    elseif strcmp(do.softwareFlag,'dcm2niix')
                        niftinames = [rawSubNames{ss} '_task-' tmpTaskName '_run-' num2str(i,formatSpec) '_bold'];
                        warning('NO DUMMY EXTRACTION PERFORMED WHEN USING DCM2NIIX');
                        str4dcm2niix = sprintf([fullfile(path2exe, 'dcm2niix'),  ' -b y -v 0 -z y -o', ' %s', ' -f', ' %s', ' -c', ' %s', ' %s'], rawSubFuncDir, niftinames, 'Dummy images included', currDir);
                        unix(str4dcm2niix);
                    end
                    
                catch ME
                    warning('Dicom in %s could not be converted.\n',folderContent(i).name)
                    fprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                        ME.stack(1).name, ME.stack(1).line, ME.message);
                end
            end
            
            if strcmp(do.softwareFlag,'SPM')
                % write JSON file if you are using SPM
                % give the function the directory where to store the
                % json file and the last dicom directory to read out
                % necessary information
                BIDS_bold_json (rawSubFuncDir,dirfiles(1,:),[rawSubNames{ss} '_task-' tmpTaskName '_bold.json']);
            end
        end
    end % do.funConversion
    
    
    %% STEP 3: Fieldmap converstion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Note: only (for now) dcm2niix implemented.
    % Note: it will rename in BIDS conform only for: magnitude + phase diff
    % OR EPI with inverted phase encoding direction.
    if do.fieldmapConversion         
        currDir = fullfile(sourceSubDir,fmapDir); % sourcedata dir for fmaps
        destDir = fullfile(rawSubDir, fmapDir);   % rawdata dir for fmaps
        
        for ext = 1:length(extensions)
            dirfiles = [dirfiles; spm_select('FPList', currDir, extensions{ext})];
        end
        
        if isempty(dirfiles)
            error('FOLDER CONTENT IS EMPTY - PROBABLY WRONG PATH.')
        else
            try
                fprintf('=> importing fieldmap to %s\n', destDir); % it should have only one folder (for both grep_fieldmap and different phase)
                
                if strcmp(do.softwareFlag,'dcm2niix') % only dcm2niix implemented
                    str4dcm2niix = sprintf([fullfile(path2exe, 'dcm2niix'), ' -b y -v 0 -z y -o', ' %s', '%s'], destDir, tmpDir);
                    unix(str4dcm2niix);
                elseif strcmp(do.softwareFlag,'SPM')
                    warning('SPM fieldmap convertion not implementet yet');
                else
                    warning('Not recognized software');
                end
                
                % Re-name the nifti and json files (to _magnitude1,
                % _magnitude2, _phasediff OR _epi)
                fmapImg     = spm_select('FPList',destDir,{'.json','nii.gz'});
                if length(fmapImg) > 1
                    for i = 1:length(fmapImg)
                        if contains(fmapImg(i,:),'_e1.nii.gz') % magnitude 1
                            movefile(fmapImg(i,:), fullfile(destDir,[rawSubNames{ss} '_magnitude1.nii.gz']));
                        elseif contains(fmapImg(i,:),'_e2.nii.gz') % magnitude 2
                            movefile(fmapImg(i,:), fullfile(destDir,[rawSubNames{ss} '_magnitude2.nii.gz']));
                        elseif contains(fmapImg(i,:),'_ph.nii.gz') % phase difference
                            movefile(fmapImg(i,:), fullfile(destDir,[rawSubNames{ss} '_phasediff.nii.gz']));
                        elseif contains(fmapImg(i,:),'_e1.json') % magnitude 1 json
                            movefile(fmapImg(i,:), fullfile(destDir,[rawSubNames{ss} '_magnitude1.json']));
                        elseif contains(fmapImg(i,:),'_e2.json') % magnitude 2 json
                            movefile(fmapImg(i,:), fullfile(destDir,[rawSubNames{ss} '_magnitude2.json']));
                        elseif contains(fmapImg(i,:),'_ph.json') % phase diff json
                            movefile(fmapImg(i,:), fullfile(destDir,[rawSubNames{ss} '_phasediff.json']));
                        elseif contains(fmapImg(i,:),'fmap_pepolar') % collected using inverted phase encoding
                            if contains(fmapImg(i,:),'.nii.gz') % nifti --> should be called '_epi.nii.gz'
                                movefile(fmapImg(i,:), fullfile(destDir,[rawSubNames{ss} 'dir_PA_epi.nii.gz']));
                            elseif contains(fmapImg(i,:),'.json') % sidecart json file defining what it is intended for.
                                movefile(fmapImg(i,:), fullfile(destDir,[rawSubNames{ss} 'dir_PA_epi.json']));
                            end
                        end
                    end
                end  
                
            catch ME
                warning('Fieldmap dicom in %s could not be converted.\n',currDir)
                fprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                    ME.stack(1).name, ME.stack(1).line, ME.message);
            end
        end
    end % Fieldmap Convertion
end % For all subjects
end % function