function DICOMconversionSPM12_BIDS()

%% DICOM to NifTI conversion using spm12
% VP 10/2020 (adapted from JBs preprocessingSPM12

%% Set root directory
%-------------------------------------------------------------------------%
% DEFINE path and get relevant filenames.
%-------------------------------------------------------------------------%
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "source_data")'])
root_dir    = uigetdir(homedir, 'Select Project Folder');
% stop the script when no folder was selected
if root_dir == 0
    error('No folder was selected --> I terminate the script')
end
source_dir  = fullfile(root_dir, 'source_data');
if ~isfolder(source_dir)
    fprintf(['It appears you do not have a "source_data" folder\n'...
        'Please select the folder that contains your dicoms'])
    source_dir  = uigetdir(root_dir, 'Select folder containing DICOMS');
    if source_dir == 0
        error('No folder was selected --> I terminate the script')
    end
end

prefix = input (['Please specify the prefix of your participant data.\n' ...
    '(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');

%-------------------------------------------------------------------------%
% DEFINE file extensions you want sort
extensions = {'**.IMA','**.ima'}; % extension you care about
% IMPORTANT: It seems like dir() (at least on MacOS) is not case sensitive,
% but spm_select('FPList',...) is
%-------------------------------------------------------------------------%

nDummies = input(['Please specify the number of dummy images.\n'...
    '!!!IMPORTANT!!!\n'...
    'Be super duper sure to give the right number of images, because those will be skipped in the conversion process!\n\n']);

% Safetynet, that ensures the input for dummy images is actually numeric
correctInput = isnumeric(nDummies);
while ~correctInput
    nDummies = input(['It seems your input for the number of dummy images is not numeric.\n'...
        'Please specify the number of dummy images again (do not forget to be super duper sure).\n']);
    correctInput = isnumeric(nDummies);
end

dic_struct_dir  = 'anat/';

raw_dir         = fullfile(root_dir, 'raw_data');

%% Decide what to do
%..............................WHAT TO DO.................................%
do.overwrite            = 0;
do.func_conversion      = 1;
do.struct_conversion    = 0;

% TODO: create a log file to save warning and error messages

%% OPEN SPM
spm fmri;

% create a BIDS conform directory structure for the NIFTIS
% first we need to create a cell containing the subject names
folders = dir(fullfile(source_dir,[prefix, '*']));
subNames = {folders(:).name};
% source data: DICOMS
% raw data: NIFTIS unprocessed
spm_mkdir(raw_dir, subNames, dic_struct_dir);
spm_mkdir(raw_dir, subNames, 'func'); 

% write the dataset description json file
BIDS_dataset_json(raw_dir);
%% start to perform the conversion
for ss = 1:length(subNames) % For all subjects do each ...
    
    nifti_dir = fullfile(raw_dir,subNames{ss});
    func_nifti_dir = fullfile(nifti_dir,'func');

    % where to find functional data
    dicom_dir       = fullfile(source_dir,subNames{ss});
    func_dicom_dir  = fullfile(dicom_dir,'func');
    
    %% Conversion from functional DICOM to NIfTI
    if do.func_conversion
        folderContent   = dir(fullfile(func_dicom_dir,'run*'));
        %% create a BIDS conform file structure for every subject
        % source data: DICOMS
        spm_mkdir (raw_dir, subNames{ss}, 'func', {folderContent(:).name});
        % ......DICOM to NIFTI Conversion...... %
        % Get folder Content
        
        if isempty(folderContent) % no dicoms in folder
            warning('FOLDER CONTENT IS EMPTY - PROBABLY WRONG PATH: process stopped. Press Enter to continue.')
            disp(['SUBJECTS PATH:  ' subj(ss).fol_dicom])
            pause;
        else
            fprintf('FOLDER CONTENT FOUND \n')
            fprintf('========STARTING CONVERTION FROM NIFTI TO DICOM========\n\n');
            
            for i = 1:length(folderContent)            % loop through all folders found
                try
                    curr_dir = fullfile(func_dicom_dir, folderContent(i).name);
                    dest_dir = fullfile(func_nifti_dir, folderContent(i).name);
                    
                    dirfiles = [];
                    for ext = 1:length(extensions)
                        dirfiles = [dirfiles; spm_select('FPList', curr_dir, extensions{ext})];
                    end
                    if isempty(dirfiles)
                        warning('NO FILES WITH YOUR SPECIFIED FILEEXTENTION SELECTED - PROBABLY WRONG PATH/FILENAME: process stopped. Press Enter to continue.')
                        disp(['CURRENT PATH:  ' curr_dir]);
                        pause;
                    end
                    
                    % specify spm options
                    matlabbatch{1}.spm.util.import.dicom.data               = cellstr(dirfiles(nDummies+1:end,:));
                    matlabbatch{1}.spm.util.import.dicom.root               = 'flat';
                    matlabbatch{1}.spm.util.import.dicom.outdir             = cellstr(dest_dir);
                    matlabbatch{1}.spm.util.import.dicom.protfilter         = '.*';
                    matlabbatch{1}.spm.util.import.dicom.convopts.format    = 'nii';
                    matlabbatch{1}.spm.util.import.dicom.convopts.meta      = 0;
                    matlabbatch{1}.spm.util.import.dicom.convopts.icedims   = 0;
                    
                    fprintf('=> importing dicoms from %s\n', curr_dir);
                    fprintf('                      to %s\n', dest_dir);
                    
                    spm_jobman('run', matlabbatch);
                    
                catch ME
                    warning('Dicom in %s could not be converted.\n',folderContent(i).name)
                    sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                        ME.stack(1).name, ME.stack(1).line, ME.message);
                end
                % write JSON file
                % give the function the directory where to store the
                % json file and the last dicom directory to read out
                % necessary information
                BIDS_bold_json (func_nifti_dir,dirfiles(1,:),[subNames{ss} '_bold.json']);
            end
        end
    end
    
    %% Conversion from DICOM to NIfTI of STRUCTURAL image
    % (just to have it in another folder and change name)
    if do.struct_conversion
        curr_dir = fullfile(dicom_dir,'anat');
        dest_dir = fullfile(nifti_dir, dic_struct_dir);
        
        % select files (either "ima" or "IMA")
        dirfiles = [];
        for ext = 1:length(extensions)
            dirfiles = [dirfiles; spm_select('FPList', curr_dir, extensions{ext})];
        end
        if isempty(dirfiles)
            warning('NO FILES WITH YOUR SPECIFIED FILEEXTENTION SELECTED - PROBABLY WRONG PATH/FILENAME: process stopped. Press Enter to continue.')
            disp(['CURRENT PATH:  ' curr_dir]);
            pause;
        end
        
        fprintf('=> importing structural to %s\n', dest_dir);
        % specify spm options
        matlabbatch{1}.spm.util.import.dicom.data               = cellstr(dirfiles);
        matlabbatch{1}.spm.util.import.dicom.root               = 'flat';
        matlabbatch{1}.spm.util.import.dicom.outdir             = cellstr(dest_dir);
        matlabbatch{1}.spm.util.import.dicom.protfilter         = '.*';
        matlabbatch{1}.spm.util.import.dicom.convopts.format    = 'nii';
        matlabbatch{1}.spm.util.import.dicom.convopts.meta      = 0;
        matlabbatch{1}.spm.util.import.dicom.convopts.icedims   = 0;
        
        spm_jobman('run', matlabbatch);
        
        % After creating the anatomical image, creat its corresponding json
        % file
        BIDS_anatT1w_json(dest_dir, dirfiles(1,:),[subNames{ss} '_T1w']);
    end
end