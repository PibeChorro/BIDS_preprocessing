function preprocessingSPM12_BIDS()

%% Preprocess fMRI data using spm12
% JB 08/2019 (adapted from PRG 05/2019)
% Further changed by VP 01/2021

%% Define important details of your file structure and location
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "raw_data")\n\n'])
rootDir    = '/Users/vpl/Documents/Master_Thesis/DATA/MRI_II'; % uigetdir(homedir, 'Select Project Folder');
if rootDir == 0
    error('No folder was selected --> I terminate the script')
end

% Set source_data directory. This is only needed for slice time correction
sourceDir  = fullfile(rootDir,'sourcedata');
if ~isfolder(sourceDir)
    fprintf(['It appears you do not have a "source_data" folder.\n'...
        'Please select the folder that contains your unprocessed niftis.\n\n'])
    sourceDir  = uigetdir(rootDir, 'Select DICOM folder');
    if sourceDir == 0
        error('No folder was selected --> I terminate the script')
    end
end

DICOMprefix = 'sMag'; % input (['Please specify the prefix of your participant data in your SOURCEDATA directory.\n' ...
%     '(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n\n'],'s');

DICOMfolders = dir(fullfile(sourceDir,[DICOMprefix, '*']));
DICOMsubNames = {DICOMfolders(:).name}; 

% Set raw_data directory.
rawDir     = fullfile(rootDir, 'rawdata');
if ~isfolder(rawDir)
    fprintf(['It appears you do not have a "raw_data" folder.\n'...
        'Please select the folder that contains your unprocessed niftis.'])
    rawDir  = uigetdir(rootDir, 'Select unprocessed nifti folder');
    if rawDir == 0
        error('No folder was selected --> I terminate the script')
    end
end

derivativesDir = fullfile (rootDir, 'derivatives');
if ~isfolder(derivativesDir)
    mkdir (derivativesDir);
end

sub_prefix      = 'sub-';   % subject prefix - this one is REQUIRED
session_prefix  = 'ses';    % session prefix - only needed if more than one session per subject
session_labels  = '';       % labels of sessions - must be alphanumeric
task_prefix     = 'task-';  % this one is REQUIRED (it does not make too much sense, when having only one task, but it gives more information in the name
task_labels     = {'magic'};  % labels of tasks
anat_dir        = 'anat';   % where the structural image is stored
anat_modality   = 'T1w';

%% Set name of data
% A first draft of how one could automatise the names of data
% sessionNames = {};
% if isempty(session_labels)
%     sessionNames = '';
% else
%     for ses = 1:length(session_labels)
%         sessionNames{end+1} = [session_prefix session_labels{ses}];
%     end
% end
% 
% taskNames = {};
% for task = 1:length(task_labels)
%     taskNames{end+1} = [task_prefix task_labels{task}];
% end
% 
% dataNames = {};
% for ses = 1:length(sessionNames)
%     for task = 1:length(taskNames)
%         dataNames{end+1} = [sessionNames{sees} '_' taskNames{task}];
%     end
% end

taskName = [task_prefix task_labels{1}];

%-------------------------------------------------------------------------%
% DEFINE file extensions of DICOMs you care about
extensions = {'**.IMA','**.ima'}; % extension you care about
% IMPORTANT: It seems like dir() (at least on MacOS) is not case sensitive,
% but spm_select('FPList',...) is
%-------------------------------------------------------------------------%


%% Decide what to do
%..............................WHAT TO DO.................................%
do.overwrite        = 1;
do.realignment      = 1; % 1 = realigning and unwarp;
do.slice_time_corr  = 1; % slice time correction (using slice TIMES); 
do.coregistration   = 'auto'; % 'manual' or 'auto';
do.segmentation     = 1;
do.normalisation    = 1; 
do.smoothing        = 1; %Smoothing Flag, set to 1 if you want to smooth. 
do.smoothNorm       = 'mni'; % Smooth normalize data = 'mni', native data = 'native' or 'both'
do.smoothingSize    = 6;

% already assign realignment parameter names
raParamNames = {'x-Axis', 'y-Axis', 'z-Axis',...
    'Pitch','Roll','Yaw'};

%% OPEN SPM
spm fmri;

% create a BIDS conform directory structure for the NIFTIS
% first we need to create a cell containing the subject names
softwareName = 'spm12';
pipelineName = 'spm12-preproc';
folders = dir(fullfile(rawDir,[sub_prefix, '*']));
subNames = {folders(:).name}; 

%% start to perform the preprocessing
for ss = 1:1% length(subNames) % For all subjects do each ...
    % get the unprocessed niftis
    raw_nifti_dir = fullfile(rawDir,subNames{ss});
    raw_func_nifti_dir = fullfile(raw_nifti_dir,'func');
    
    % check if structural volume exists if needed
    dir_raw_structImg   = dir(fullfile(raw_nifti_dir, anat_dir, ['*' anat_modality '.nii']));
    dir_raw_structImg   = fullfile(dir_raw_structImg.folder, dir_raw_structImg.name);
    
    %% create a BIDS conform file structure for every subject
    % !!!only the derivative_data folder is created here. The rest
    % (raw_data and source_data) already needs to be like this !!!
    % 
    % project/
    %   derivative_data/
    %       <software-name>         // spm12 in this case
    %           <pipeline-name>         // spm12-preproc in this case
    %               <processing-step1>
    %               <processing-step2>
    %                   sub<nr>/
    %                       anat/
    %                       func/
    %                           run<nr>/
    %               ...
    %   raw_data/
    %       sub<nr>/
    %           anat/
    %           func/
    %               run<nr>
    %   source_data
    %       sub<nr>/
    %           anat/
    %           func/
    %               run<nr>
    %       
    % create a directory name for all preprocessing steps
    realigned_dir               = fullfile (derivativesDir, softwareName, pipelineName, 'realigned/');
    slice_time_corrected_dir    = fullfile (derivativesDir, softwareName, pipelineName, 'slice_time_corrected/');
    coregistered_dir            = fullfile (derivativesDir, softwareName, pipelineName, 'coregistered/');
    normalized_dir              = fullfile (derivativesDir, softwareName, pipelineName, 'normalized/');
    smoothed_dir                = fullfile (derivativesDir, softwareName, pipelineName, [num2str(do.smoothingSize) 'mm-smoothed/']);
    segmented_dir               = fullfile (derivativesDir, softwareName, pipelineName, 'segmented/');
    % establish BIDS conform data structure for each step
    spm_mkdir (realigned_dir,               subNames{ss}, 'func');
    spm_mkdir (slice_time_corrected_dir,    subNames{ss}, 'func');
    spm_mkdir (coregistered_dir,            subNames{ss}, 'func');
    spm_mkdir (segmented_dir,               subNames{ss}, 'anat');
    spm_mkdir (normalized_dir,              subNames{ss}, 'func');
    spm_mkdir (normalized_dir,              subNames{ss}, 'anat');
    spm_mkdir (smoothed_dir,                subNames{ss}, 'func');
    
    %% STARTING PREPROCESSING
    %% Realignment: Estimate & unwarp
    if do.realignment
        
        % getting the raw functional NIfTIs out of the 'raw_data' directory
        folderContent = dir(fullfile(raw_func_nifti_dir,[subNames{ss} '_' taskName '_' 'run*'])); 
        nruns = length(folderContent);
        fprintf('STARTING REALIGNMENT AND UNWARPING \n\n')
        
        % iterate over runs
        for run = 1:nruns 
            dirfiles     = spm_select('ExtFPList', raw_func_nifti_dir, folderContent(run).name, Inf);
            
            if strcmp(dirfiles,'')
                warning('No files selected!');
                return;
            end
            
            % Get files
            matlabbatch{1}.spm.spatial{1}.realignunwarp.data(run).scans  = cellstr(dirfiles);
            
            % No field map is being used
            matlabbatch{1}.spm.spatial{1}.realignunwarp.data(run).pmscan = {''} ;
        end
        
        % Specify Estimation options
        matlabbatch{1}.spm.spatial{1}.realignunwarp.eoptions.quality    = 0.9;       % Quality vs speed trade-off (1 = highest quality)
        matlabbatch{1}.spm.spatial{1}.realignunwarp.eoptions.sep        = 4;         % Sampling of reference in mm (smaller is better but slower)
        matlabbatch{1}.spm.spatial{1}.realignunwarp.eoptions.fwhm       = 5;         % Smoothing kernel before realignment (5 mm typical for MRI)
        matlabbatch{1}.spm.spatial{1}.realignunwarp.eoptions.rtm        = 0;         % 1 = Register to mean (2-pass), 0 = Register to first only
        matlabbatch{1}.spm.spatial{1}.realignunwarp.eoptions.einterp    = 2;         % 2nd-degree B-spline interpolation
        matlabbatch{1}.spm.spatial{1}.realignunwarp.eoptions.ewrap      = [0 0 0];   % Y-Wrapping % PRG [0 0 0]
        matlabbatch{1}.spm.spatial{1}.realignunwarp.eoptions.weight     = '';        % No weighting of voxels in realignment process
        
        % Specify Unwarp Estimation options
        matlabbatch{1}.spm.spatial{1}.realignunwarp.uweoptions.basfcn   = [12 12];   % Basis function for each dimension (3rd dimension left open)
        matlabbatch{1}.spm.spatial{1}.realignunwarp.uweoptions.regorder = 1;         % Regularisation
        matlabbatch{1}.spm.spatial{1}.realignunwarp.uweoptions.lambda   = 100000;    % Regularisation factor = medium
        matlabbatch{1}.spm.spatial{1}.realignunwarp.uweoptions.jm       = 0;         % No distortion-based intensity correction (because not a good idea apparently)
        matlabbatch{1}.spm.spatial{1}.realignunwarp.uweoptions.fot      = [4 5];     % 1st order effects, model only Pitch & Roll (1:6 for all movements)
        matlabbatch{1}.spm.spatial{1}.realignunwarp.uweoptions.sot      = [];        % No 2nd order effects modelled
        matlabbatch{1}.spm.spatial{1}.realignunwarp.uweoptions.uwfwhm   = 4;         % Smoothing kernel for unwarp
        matlabbatch{1}.spm.spatial{1}.realignunwarp.uweoptions.rem      = 1;         % Movement parameters reestimated after every unwarp iteration
        matlabbatch{1}.spm.spatial{1}.realignunwarp.uweoptions.noi      = 5;         % Maximum number of iterations
        matlabbatch{1}.spm.spatial{1}.realignunwarp.uweoptions.expround = 'Average'; % Point around which to perform Taylor-expansion
        
        % Specify Reslice options
        matlabbatch{1}.spm.spatial{1,1}.realignunwarp.uwroptions.uwwhich  = [2 1];     % Create mean image & all images resliced
        matlabbatch{1}.spm.spatial{1,1}.realignunwarp.uwroptions.rinterp  = 4;         % 4th degree B-spline interpolation (for high-res 'Inf' = Fourier interpolation?)
        matlabbatch{1}.spm.spatial{1,1}.realignunwarp.uwroptions.wrap     = [0 0 0];   % Y-Wrapping % PRG [0 0 0]
        matlabbatch{1}.spm.spatial{1,1}.realignunwarp.uwroptions.mask     = 1;         % Use masking (search for voxels that cannot be sampled)
        matlabbatch{1}.spm.spatial{1,1}.realignunwarp.uwroptions.prefix   = 'u';
        
        fprintf('=> realigning, UNWARPING and reslicing\n');
        
        % run realign & unwarp
        spm_jobman('run', matlabbatch);
        clear matlabbatch;
        
        % move created realligned NIfTIs and other files (the mean EPI
        % NIfTI, the realignment text-files and the .mat file that is
        % created in the process from "raw_data" to "derivative_data"
        % move all the single functional niftis
        [success,message] = movefile(string(fullfile(raw_func_nifti_dir, 'u*')), ...
            fullfile(realigned_dir, subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
        % move the mean realigned nifti immage
        [success,message] = movefile(string(fullfile(raw_func_nifti_dir, 'meanu*')), ...
            fullfile(realigned_dir,subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
        % move the realignment parameter .txt files
        [success,message] = movefile(string(fullfile(raw_func_nifti_dir, 'rp*.txt')), ...
            fullfile(realigned_dir, subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
        
        % move the created .mat files
        [success,message] = movefile(string(fullfile(raw_func_nifti_dir, '*.mat')), ...
            fullfile(realigned_dir, subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
        
        ra_files = dir(fullfile(realigned_dir, subNames{ss}, 'func', '*.txt'));
        
        for run = 1:nruns 
            % read in the created realignment text file and plot them in a
            % firgure and save the figure 
            try
                [fid, mes] = fopen(fullfile(ra_files(run).folder,ra_files(run).name));
                realignmentMatrix = textscan(fid, '%f%f%f%f%f%f');
                % create figure to plot realigment parameters in
                fig = figure;
                % create a subplot for the tranlation
                subplot(2,1,1);
                hold on
                for param = 1:3
                    plot(realignmentMatrix{param},'DisplayName',raParamNames{param});
                end
                hold off
                legend();       % NOTICE - this is broken. The legend shows the first label correctly and the rest overlaps
                
                % create a subplot for the rotation
                subplot(2,1,2);
                hold on
                for param = 4:6
                    plot(realignmentMatrix{param},'DisplayName',raParamNames{param});
                end
                hold off
                legend();       % NOTICE - this is broken. The legend shows the first label correctly and the rest overlaps
                
                % save and close figure. Close realignment file
                savefig(fullfile(realigned_dir, subNames{ss}, 'func', ['realignmentPlot_run-' num2str(run)]));
                close(fig);
                fclose(fid);
            catch
                warning(mes)
            end
        end
        % TODO: create JSON file containing processing information and
        % store it BIDS conform
    end
    
    %% Slice Time Correction
    % slice time correction for every run individually because slice timing
    % differs slightly
    if do.slice_time_corr
        currentDir          = fullfile(realigned_dir, subNames{ss}, 'func');
        folderContent       = dir(fullfile(currentDir, ['u' subNames{ss} '_' taskName '_' 'run*'])); 
        DICOMfolderContent  = dir(fullfile(sourceDir, DICOMsubNames{ss}, 'func', 'run*'));
        nruns               = length(folderContent);
        fprintf('SLICE TIME CORRECTION\n\n')
        for run=1:nruns % for number of runs
            % load a dicom header that contains information needed for analysis
            dicom_dir           = fullfile(sourceDir,DICOMsubNames{ss},'func',DICOMfolderContent(run).name);
            % select all dicom files from the current run
            dicom_files = [];
            for ext = 1:length(extensions)
                dicom_files = [dicom_files; spm_select('FPList', dicom_dir, extensions{ext})];
            end
            if isempty(dicom_files)
                % TODO: check what is going on here!
                warning('NO *ima nor *IMA  FILES SELECTED FOR SLICE TIME CORRECTION - PROBABLY WRONG PATH/FILENAME: process stopped. Press Enter to continue.')
                disp(['CURRENT PATH:  ' dicom_dir]);
                pause;
            end

            % select only one from them
            dicom_file   = dicom_files(end,:); % Get one image, e.g. the last one
            fprintf('=> determining acquisition parameters from: \n %s \n', dicom_file);

            hdr                         = spm_dicom_headers(dicom_file);
            sequence.N_slices           = hdr{1}.Private_0019_100a;         % Number of slices
            sequence.TR                 = hdr{1}.RepetitionTime/1000;       % TR
            [~, sequence.slice_order]   = sort(hdr{1}.Private_0019_1029);   % slice order
            sequence.sliceTstamps       = hdr{1}.Private_0019_1029;         % slice times in milliseconds
            
            % select files to slice time correct
            nifti_files = spm_select('ExtFPList', currentDir, folderContent(run).name, Inf); % AFTER REALINGMENT!
            
            % find reference slice (the one in the middle) - if using
            % median one has two values --> select one of them using "min"
            tmp = abs(sequence.sliceTstamps - median(sequence.sliceTstamps));   % find the slice time in the middle! (in ms)
            [~, refInd] = min(tmp);
            
            % sanity check TR.
            if max(sequence.sliceTstamps)/1000 > sequence.TR 
                error('Found a slice timestamp that exceeds our TR! Make sure you did not select the first image in you run folder to assess slice timing!')
            end
            
            % sanity check Number of slices
            if length(sequence.sliceTstamps) ~= sequence.N_slices
                error('The slice time stamps found in the DICOM header do not correspond to Nslices!');
            end
            
            % sanity check Middle slice. if the slice you found is one of the middle slices.
            if logical(refInd ~= sequence.slice_order(sequence.N_slices/2)) && logical(refInd ~= sequence.slice_order(sequence.N_slices/2+1))
                error('Problem finding you middle slice! It doesnt correspond to the N_slices/2 nor (N_slices/2)+1')
            end
            
            % generate matlabbatch for slice time correction with everthing
            matlabbatch{1}.spm.temporal.st.scans      = {cellstr(nifti_files)}; % nifti files
            matlabbatch{1}.spm.temporal.st.nslices    = sequence.N_slices; % nr. slices
            matlabbatch{1}.spm.temporal.st.tr         = sequence.TR; % TR
            matlabbatch{1}.spm.temporal.st.ta         = 0; % will be ignored because we use slice_times (from SPM: "if the next two items (slice order & reference slice) are entered in milliseconds, this entry will not be used and can be set to 0")
            matlabbatch{1}.spm.temporal.st.so         = sequence.sliceTstamps;  % slice time stamps (in ms)
            matlabbatch{1}.spm.temporal.st.refslice   = sequence.sliceTstamps( refInd ); % reference time stamp (in ms), the one in the middle
            matlabbatch{1}.spm.temporal.st.prefix     = 'a';
            
            jobs = matlabbatch;
            spm('defaults', 'FMRI');
            spm_jobman('run', jobs);
            clearvars matlabbatch
            
        end
        [success,message] = movefile(string(fullfile(realigned_dir, subNames{ss}, 'func', 'a*')), ...
            fullfile(slice_time_corrected_dir, subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
    end
    
    %% Co-registration
    if do.coregistration
        % copy the slice time corrected images in the coregistration folder
        % then perform the coregistration on the copied images
        % TODO: add a prefix ('c') to files in the coregistation folder
        folderContent = dir(fullfile(slice_time_corrected_dir, subNames{ss}, 'func', ['au' subNames{ss} '_' taskName '_' 'run*'])); 
        nruns               = length(folderContent);
        for run = 1:length(folderContent)
            [success,message] = copyfile(fullfile(slice_time_corrected_dir, subNames{ss}, 'func', folderContent(run).name),...
                fullfile(coregistered_dir, subNames{ss}, 'func'));
            if ~success
                warning(message)
            end
        end
        
        % also copy the mean EPI image from the realignment folder into the
        % coregistation folder
        [success,message] = copyfile(fullfile(realigned_dir, subNames{ss}, 'func',['meanu' '*.nii']),...
            fullfile(coregistered_dir, subNames{ss}, 'func'));
        if ~success
                warning(message)
        end
        
        meanEpi = dir(fullfile(coregistered_dir, subNames{ss}, 'func' ,['meanu' '*.nii']));
        meanEpi = fullfile(coregistered_dir, subNames{ss}, 'func', meanEpi.name);
        if strcmpi(do.coregistration, 'manual')
            % TODO: this is not up to date and won't work
            fprintf('MANUAL COREGISTRATION\n')
            fprintf('EPI scans [meanEPI] -> Structural \n');
            
            mancoreg(cellstr(dir_raw_structImg),sourceimage)
            
        elseif strcmpi(do.coregistration, 'auto') || do.coregistration 
            
            fprintf('AUTOMATIC COREGISTRATION\n')
            fprintf('EPI scans [meanEPI] -> Structural \n');
            
            alltargets = {}; %JB It's necessary to initialize alltargets as a cell array, as this can prevent the vertcat error if paths from different runs have different character lengths.
            for run = 1:nruns
                dir_sessdata = fullfile(coregistered_dir, subNames{ss}, 'func');
                dirfiles     = spm_select('ExtFPList', dir_sessdata, folderContent(run).name, Inf);
                if strcmp(dirfiles,'')
                    warning('No files selected!');
                    return;
                end
                alltargets = [alltargets; cellstr(dirfiles)];
            end
            
            % Get files
            matlabbatch{1}.spm.spatial.coreg.estimate.ref               = cellstr(dir_raw_structImg);           % ANATOMICAL SCAN
            matlabbatch{1}.spm.spatial.coreg.estimate.source            = cellstr(meanEpi);                     % Mean EPI = source
            matlabbatch{1}.spm.spatial.coreg.estimate.other             = cellstr(alltargets);                  % Other files to be moved (all the realigned EPIs)
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';                                % Normalized mutual information
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];                                % Sampling in mm, coarse to fine
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 ... % Accuracy for each parameter
                0.01 0.01 0.01 0.001 0.001 0.001];   % Iterations stop when less than tolerance
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];                                % Gaussian smoothing to be applied

            jobs = matlabbatch;
            spm('defaults', 'FMRI');
            spm_jobman('run', jobs);
            clearvars matlabbatch
        end
    end
    
    %% Segmentation of anatomical image
    if do.segmentation
        % TODO: change directory of spm .nii images to a more "general"
        % directory
        
        fprintf('SEGMENTATION\n\n')
        matlabbatch{1}.spm.spatial.preproc.channel.vols     = cellstr(dir_raw_structImg);
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm    = {'/Users/jasper/spm12/tpm/TPM.nii,1'};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus  = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm    = {'/Users/jasper/spm12/tpm/TPM.nii,2'};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus  = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm    = {'/Users/jasper/spm12/tpm/TPM.nii,3'};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus  = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm    = {'/Users/jasper/spm12/tpm/TPM.nii,4'};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus  = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm    = {'/Users/jasper/spm12/tpm/TPM.nii,5'};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus  = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm    = {'/Users/jasper/spm12/tpm/TPM.nii,6'};
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus  = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf         = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup     = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg      = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm        = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp        = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write       = [1 1];
        
        jobs = matlabbatch;
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs);
        
        % move the segmentation images into the derivative folder
        [success,message] = movefile(fullfile(raw_nifti_dir, anat_dir, 'c*.nii'),...
            fullfile(segmented_dir, subNames{ss}, 'anat'));
        if ~success
                warning(message)
        end
        [success,message] = movefile(fullfile(raw_nifti_dir, anat_dir, 'i*.nii'),...
            fullfile(segmented_dir, subNames{ss}, 'anat'));
        if ~success
                warning(message)
        end
        [success,message] = movefile(fullfile(raw_nifti_dir, anat_dir, 'y*.nii'),...
            fullfile(segmented_dir, subNames{ss}, 'anat'));
        if ~success
                warning(message)
        end
        [success,message] = movefile(fullfile(raw_nifti_dir, anat_dir, '*.mat'),...
            fullfile(segmented_dir, subNames{ss}, 'anat'));
        if ~success
                warning(message)
        end
        
        clearvars matlabbatch
    end
    
    %% Normalisation: write only (uses deformations field from segmentation)
    if do.normalisation
        
        fprintf('NORMALISATION USING SEGMENTATION\n\n')
        
        struct_defForward   = dir(fullfile (segmented_dir, subNames{ss}, 'anat','y_*.nii'));
        folderContent       = dir(fullfile(coregistered_dir, subNames{ss}, 'func', ['au' subNames{ss} '_' taskName '_' 'run*.nii']));
        nruns               = length(folderContent);
        
        meanEpi = dir(fullfile(coregistered_dir, subNames{ss}, 'func', ['meanu' '*.nii']));
        meanEpi = fullfile(coregistered_dir, subNames{ss}, 'func', meanEpi.name);
        
        alltargets = {}; %JB It's necessary to initialize alltargets as a cell array, as this can prevent the vertcat error if paths from different runs have different character lengths.
        for r = 1:nruns
            dir_sessdata = fullfile(coregistered_dir, subNames{ss}, 'func');
            dirfiles     = spm_select('ExtFPList', dir_sessdata, folderContent(r).name, Inf);
            if strcmp(dirfiles,'')
                warning('No files selected!');
                return;
            end
            alltargets = [alltargets; cellstr(dirfiles)]; 
        end
        
        alltargets = cellstr(alltargets);
        alltargets{end+1} = meanEpi;
        alltargets{end+1} = dir_raw_structImg; % add anatomical image to normalise
        
        matlabbatch{1}.spm.spatial.normalise.write.subj.def         = {fullfile(struct_defForward.folder, struct_defForward.name)};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample    = alltargets;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb      = [-78 -112 -70
            78  76   85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox     = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp  = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix  = 'w';
        
        jobs = matlabbatch;
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs);
        clearvars matlabbatch
        
        % move created normalized niftis in a seperate folder
        [success,message] = movefile(fullfile(coregistered_dir, subNames{ss}, 'func', 'w*.nii'),...
            fullfile(normalized_dir, subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
        
        % move the created normalized anatomical image in a seperate folder
        [success,message] = movefile(fullfile(raw_nifti_dir, anat_dir, 'w*.nii'),...
            fullfile(normalized_dir, subNames{ss}, 'anat'));
        if ~success
                warning(message)
        end
    end
    
    %% Smoothing
    if do.smoothing
        fprintf('SMOOTHING\n\n')
        folderContent = dir(fullfile(normalized_dir, subNames{ss}, 'func', ['wau' subNames{ss} '_' taskName '_' 'run*']));
        nruns         = length(folderContent);
        
        alltargets = {};
                
        for r = 1:nruns
            
            if strcmp(do.smoothNorm,'mni') == 1
                dir_sessdata = fullfile(normalized_dir, subNames{ss}, 'func');
            elseif strcmp(do.smoothNorm, 'native') == 1
                dir_sessdata = fullfile(coregistered_dir, subNames{ss}, 'func');
            elseif strcmp(do.smoothNorm, 'both') == 1
                 % TO-DO
            end
            dirfiles     = spm_select('ExtFPList', dir_sessdata, folderContent(r).name, Inf);
                
            if strcmp(dirfiles,'')
                warning('No files selected!');
                return;
            end
            alltargets = [alltargets; cellstr(dirfiles)];
        end
        
        matlabbatch{1}.spm.spatial.smooth.data = cellstr(alltargets);
        matlabbatch{1}.spm.spatial.smooth.fwhm = repmat(do.smoothingSize,1,3); % Set at the beginning.
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = ['s' num2str(do.smoothingSize)];
        
        jobs = matlabbatch;
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs);
        clearvars matlabbatch
        for r = 1:nruns
            [success,message] = movefile(fullfile(normalized_dir, subNames{ss}, 'func', ['s' num2str(do.smoothingSize) '*.nii']),...
                fullfile(smoothed_dir, subNames{ss}, 'func'));
            if ~success
                warning(message)
            end
        end
    end
end
end