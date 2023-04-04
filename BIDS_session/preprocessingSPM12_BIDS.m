function preprocessingSPM12_BIDS()

%% Preprocess fMRI data using spm12
% JB 08/2019 (adapted from PRG 05/2019)
% Further changed by VP 01/2021

%% Define important details of your file structure and location
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "raw_data")\n\n'])
rootDir    = '/Users/vpl/Documents/Master_Thesis/DATA/MRI'; % uigetdir(homedir, 'Select Project Folder');
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

% Set rawdata directory.
rawDir     = fullfile(rootDir, 'rawdata');
if ~isfolder(rawDir)
    fprintf(['It appears you do not have a "rawdata" folder.\n'...
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

subPrefix      = 'sub-';   % subject prefix - this one is REQUIRED
sessionPrefix  = 'ses';    % session prefix - only needed if more than one session per subject
sessionLabels  = '';       % labels of sessions - must be alphanumeric
taskPrefix     = 'task-';  % this one is REQUIRED (it does not make too much sense, when having only one task, but it gives more information in the name
taskLabels     = {'magic'};  % labels of tasks
anatDir        = 'anat';   % where the structural image is stored
anatModality   = 'T1w';

%% Set name of data
% A first draft of how one could automatise the names of data
% sessionNames = {};
% if isempty(sessionLabels)
%     sessionNames = '';
% else
%     for ses = 1:length(sessionLabels)
%         sessionNames{end+1} = [sessionPrefix sessionLabels{ses}];
%     end
% end
% 
% taskNames = {};
% for task = 1:length(taskLabels)
%     taskNames{end+1} = [taskPrefix taskLabels{task}];
% end
% 
% dataNames = {};
% for ses = 1:length(sessionNames)
%     for task = 1:length(taskNames)
%         dataNames{end+1} = [sessionNames{sees} '_' taskNames{task}];
%     end
% end

taskName = [taskPrefix taskLabels{1}];

%-------------------------------------------------------------------------%
% DEFINE file extensions of DICOMs you care about
DICOMsExtensions = {'**.IMA','**.ima'}; % extension you care about
% IMPORTANT: It seems like dir() (at least on MacOS) is not case sensitive,
% but spm_select('FPList',...) is
%-------------------------------------------------------------------------%


%% Decide what to do
%..............................WHAT TO DO.................................%
do.overwrite        = 1;
do.realignment      = 0; % 1 = realigning and unwarp;
do.sliceTimeCorr    = 1; % 1 = slice time correction (using slice TIMES); 
do.coregistration   = ''; % 'manual' or 'auto';
do.segmentation     = 0;
do.normalisation    = 0; 
do.smoothing        = 1; %Smoothing Flag, set to 1 if you want to smooth. 
do.smoothNorm       = 'mni'; % Smooth normalize data = 'mni', native data = 'native' or 'both'
do.smoothingSize    = 9; % in mm 

% already assign realignment parameter names
raParamNames = {'x-Axis', 'y-Axis', 'z-Axis',...
    'Pitch','Roll','Yaw'};

%% OPEN SPM -- necessary ?
spm fmri;

% create a BIDS conform directory structure for the NIFTIS
% first we need to create a cell containing the subject names
softwareName    = 'spm12';
pipelineName    = 'spm12-preproc';
folders         = dir(fullfile(rawDir,[subPrefix, '*']));
subNames        = {folders(:).name}; 

%% start to perform the preprocessing
for ss = 1:length(subNames) % For all subjects do each ...
    % get the unprocessed niftis
    rawSubDir       = fullfile(rawDir,subNames{ss});
    rawSubFuncDir   = fullfile(rawSubDir,'func');
    
    % check if structural volume exists if needed
    rawSubAnatImg   = dir(fullfile(rawSubDir, anatDir, ['*' anatModality '.nii']));
    rawSubAnatImg   = fullfile(rawSubAnatImg.folder, rawSubAnatImg.name);
    
    %% create a BIDS conform file structure for every subject
    % !!!only the derivatives folder is created here. The rest
    % (rawdata and sourcedata) already needs to be like this !!!
    % 
    % project/
    %   derivatives/
    %       <software-name>         // spm12 in this case
    %           <pipeline-name>         // spm12-preproc in this case
    %               <processing-step1>
    %               <processing-step2>
    %                   sub<nr>/
    %                       anat/
    %                       func/
    %                           run<nr>/
    %               ...
    %   rawdata/
    %       sub<nr>/
    %           anat/
    %           func/
    %               run<nr>
    %   sourcedata
    %       sub<nr>/
    %           anat/
    %           func/
    %               run<nr>
    %       
    % create a directory name for all preprocessing steps
    realignedDir            = fullfile (derivativesDir, softwareName, pipelineName, 'realigned/');
    sliceTimeCorrectedDir   = fullfile (derivativesDir, softwareName, pipelineName, 'slice_time_corrected/');
    coregisteredDir         = fullfile (derivativesDir, softwareName, pipelineName, 'coregistered/');
    normalizedDir           = fullfile (derivativesDir, softwareName, pipelineName, 'normalized/');
    smoothedDir             = fullfile (derivativesDir, softwareName, pipelineName, [num2str(do.smoothingSize) 'mm-smoothed-' do.smoothNorm 'space/']);
    segmentedDir            = fullfile (derivativesDir, softwareName, pipelineName, 'segmented/');
    % establish BIDS conform data structure for each step
    spm_mkdir (realignedDir,            subNames{ss}, 'func');
    spm_mkdir (sliceTimeCorrectedDir,   subNames{ss}, 'func');
    spm_mkdir (coregisteredDir,         subNames{ss}, 'func');
    spm_mkdir (segmentedDir,            subNames{ss}, 'anat');
    spm_mkdir (normalizedDir,           subNames{ss}, 'func');
    spm_mkdir (normalizedDir,           subNames{ss}, 'anat');
    spm_mkdir (smoothedDir,             subNames{ss}, 'func');
    
    %% STARTING PREPROCESSING
    %% Realignment: Estimate & unwarp
    if do.realignment
        
        % getting the raw functional NIfTIs out of the 'rawdata' directory
        folderContent = dir(fullfile(rawSubFuncDir,[subNames{ss} '_' taskName '_' 'run*.nii'])); 
        nruns = length(folderContent);
        fprintf('STARTING REALIGNMENT AND UNWARPING \n\n')
        
        % iterate over runs
        for run = 1:nruns 
            dirfiles     = spm_select('ExtFPList', rawSubFuncDir, folderContent(run).name, Inf);
            
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
        % created in the process from "rawdata" to "derivatives"
        % move all the single functional niftis
        [success,message] = movefile(string(fullfile(rawSubFuncDir, 'u*')), ...
            fullfile(realignedDir, subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
        % move and then rename the mean realigned nifti immage
        [success,message] = movefile(string(fullfile(rawSubFuncDir, 'meanu*')), ...
            fullfile(realignedDir,subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
        
        [success,message] = movefile(string(fullfile(realignedDir,subNames{ss}, 'func', ['meanu' subNames{ss} '_task-' taskLabels{1} '_run-01_bold.nii'])), ...
            fullfile(realignedDir,subNames{ss}, 'func', ['meanu' subNames{ss} '_task-' taskLabels{1} '_bold.nii']),'f');
        if ~success
            warning(message)
        end
        
        % move the realignment parameter .txt files
        [success,message] = movefile(string(fullfile(rawSubFuncDir, 'rp*.txt')), ...
            fullfile(realignedDir, subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
        
        % move the created .mat files
        [success,message] = movefile(string(fullfile(rawSubFuncDir, '*.mat')), ...
            fullfile(realignedDir, subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
        
        raFiles = dir(fullfile(realignedDir, subNames{ss}, 'func', '*.txt'));
        
        for run = 1:nruns 
            % read in the created realignment text file and plot them in a
            % firgure and save the figure 
            try
                [fid, mes] = fopen(fullfile(raFiles(run).folder,raFiles(run).name));
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
                savefig(fullfile(realignedDir, subNames{ss}, 'func', ['realignmentPlot_run-' num2str(run)]));
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
    if do.sliceTimeCorr
        currentDir          = fullfile(realignedDir, subNames{ss}, 'func');
        folderContent       = dir(fullfile(currentDir, ['u' subNames{ss} '_' taskName '_' 'run*'])); 
        DICOMfolderContent  = dir(fullfile(sourceDir, DICOMsubNames{ss}, 'func', 'run*'));
        nruns               = length(folderContent);
        fprintf('SLICE TIME CORRECTION\n\n')
        for run=1:nruns % for number of runs
            % load a dicom header that contains information needed for analysis
            dicomDir           = fullfile(sourceDir,DICOMsubNames{ss},'func',DICOMfolderContent(run).name);
            % select all dicom files from the current run
            dicomFiles = [];
            for ext = 1:length(DICOMsExtensions)
                dicomFiles = [dicomFiles; spm_select('FPList', dicomDir, DICOMsExtensions{ext})];
            end
            if isempty(dicomFiles)
                % TODO: check what is going on here!
                warning('NO *ima nor *IMA  FILES SELECTED FOR SLICE TIME CORRECTION - PROBABLY WRONG PATH/FILENAME: process stopped. Press Enter to continue.')
                disp(['CURRENT PATH:  ' dicomDir]);
                pause;
            end

            % select only one from them
            dicomFile   = dicomFiles(end,:); % Get one image, e.g. the last one
            fprintf('=> determining acquisition parameters from: \n %s \n', dicomFile);

            hdr                         = spm_dicom_headers(dicomFile);
            sequence.NumSlices          = hdr{1}.Private_0019_100a;         % Number of slices
            sequence.TR                 = hdr{1}.RepetitionTime/1000;       % TR in SECONDS
            [~, sequence.sliceOrder]    = sort(hdr{1}.Private_0019_1029);   % slice order
            sequence.sliceTstamps       = hdr{1}.Private_0019_1029;         % slice times in milliseconds
            
            % select files to slice time correct
            niftiFiles = spm_select('ExtFPList', currentDir, folderContent(run).name, Inf); % AFTER REALINGMENT!
            
            % find reference slice (the one in the middle) - if using
            % median one has two values --> select one of them using "min"
            tmp = abs(sequence.sliceTstamps - median(sequence.sliceTstamps));   % find the slice time in the middle! (in ms)
            [~, refInd] = min(tmp);
            
            % sanity check TR.
            if max(sequence.sliceTstamps)/1000 > sequence.TR 
                error('Found a slice timestamp that exceeds our TR! Make sure you did not select the first image in you run folder to assess slice timing!')
            end
            
            % sanity check Number of slices
            if length(sequence.sliceTstamps) ~= sequence.NumSlices
                error('The slice time stamps found in the DICOM header do not correspond to Nslices!');
            end
            
            % sanity check Middle slice. if the slice you found is one of the middle slices.
            if logical(refInd ~= sequence.sliceOrder(sequence.NumSlices/2)) && logical(refInd ~= sequence.sliceOrder(sequence.NumSlices/2+1))
                error('Problem finding your middle slice! It doesnt correspond to the NumSlices/2 nor (NumSlices/2)+1')
            end
            
            % generate matlabbatch for slice time correction with everthing
            matlabbatch{1}.spm.temporal.st.scans      = {cellstr(niftiFiles)}; % nifti files
            matlabbatch{1}.spm.temporal.st.nslices    = sequence.NumSlices; % nr. slices
            matlabbatch{1}.spm.temporal.st.tr         = sequence.TR; % TR
            matlabbatch{1}.spm.temporal.st.ta         = 0; % will be ignored because we use slice_times (from SPM: "if the next two items (slice order & reference slice) are entered in milliseconds, this entry will not be used and can be set to 0")
            matlabbatch{1}.spm.temporal.st.so         = sequence.sliceTstamps;  % slice time stamps (in ms)
            matlabbatch{1}.spm.temporal.st.refslice   = sequence.sliceTstamps( refInd ); % reference time stamp (in ms), the one in the middle
            matlabbatch{1}.spm.temporal.st.prefix     = 'a';
            
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
            clearvars matlabbatch
            
        end
        [success,message] = movefile(string(fullfile(realignedDir, subNames{ss}, 'func', 'a*')), ...
            fullfile(sliceTimeCorrectedDir, subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
    end
    
    %% Co-registration
    if do.coregistration
        % copy the slice time corrected images in the coregistration folder
        % then perform the coregistration on the copied images
        % TODO: add a prefix ('c') to files in the coregistation folder
        folderContent = dir(fullfile(sliceTimeCorrectedDir, subNames{ss}, 'func', ['au' subNames{ss} '_' taskName '_' 'run*'])); 
        nruns               = length(folderContent);
        for run = 1:length(folderContent)
            [success,message] = copyfile(fullfile(sliceTimeCorrectedDir, subNames{ss}, 'func', folderContent(run).name),...
                fullfile(coregisteredDir, subNames{ss}, 'func'));
            if ~success
                warning(message)
            end
        end
        
        % also copy the mean EPI image from the realignment folder into the
        % coregistation folder
        [success,message] = copyfile(fullfile(realignedDir, subNames{ss}, 'func',['meanu' '*.nii']),...
            fullfile(coregisteredDir, subNames{ss}, 'func'));
        if ~success
                warning(message)
        end
        
        meanEpi = dir(fullfile(coregisteredDir, subNames{ss}, 'func' ,['meanu' '*.nii']));
        meanEpi = fullfile(coregisteredDir, subNames{ss}, 'func', meanEpi.name);
        if strcmpi(do.coregistration, 'manual')
            % TODO: this is not up to date and won't work
            fprintf('MANUAL COREGISTRATION\n')
            fprintf('EPI scans [meanEPI] -> Structural \n');
            
            mancoreg(cellstr(rawSubAnatImg),sourceimage)
            
        elseif strcmpi(do.coregistration, 'auto') || do.coregistration 
            
            fprintf('AUTOMATIC COREGISTRATION\n')
            fprintf('EPI scans [meanEPI] -> Structural \n');
            
            alltargets = {}; %JB It's necessary to initialize alltargets as a cell array, as this can prevent the vertcat error if paths from different runs have different character lengths.
            for run = 1:nruns
                sessionDir = fullfile(coregisteredDir, subNames{ss}, 'func');
                dirfiles     = spm_select('ExtFPList', sessionDir, folderContent(run).name, Inf);
                if strcmp(dirfiles,'')
                    warning('No files selected!');
                    return;
                end
                alltargets = [alltargets; cellstr(dirfiles)];
            end
            
            % Get files
            matlabbatch{1}.spm.spatial.coreg.estimate.ref               = cellstr(rawSubAnatImg);               % ANATOMICAL SCAN
            matlabbatch{1}.spm.spatial.coreg.estimate.source            = cellstr(meanEpi);                     % Mean EPI = source
            matlabbatch{1}.spm.spatial.coreg.estimate.other             = cellstr(alltargets);                  % Other files to be moved (all the realigned EPIs)
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';                                % Normalized mutual information
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];                                % Sampling in mm, coarse to fine
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 ... % Accuracy for each parameter
                0.01 0.01 0.01 0.001 0.001 0.001];                                                              % Iterations stop when less than tolerance
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];                                % Gaussian smoothing to be applied

            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
            clearvars matlabbatch
        end
    end
    
    %% Segmentation of anatomical image
    if do.segmentation
        % TODO: change directory of spm .nii images to a more "general"
        % directory
        
        fprintf('SEGMENTATION\n\n')
        matlabbatch{1}.spm.spatial.preproc.channel.vols     = cellstr(rawSubAnatImg);
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
        
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        
        % move the segmentation images into the derivative folder
        [success,message] = movefile(fullfile(rawSubDir, anatDir, 'c*.nii'),...
            fullfile(segmentedDir, subNames{ss}, 'anat'));
        if ~success
                warning(message);
        end
        [success,message] = movefile(fullfile(rawSubDir, anatDir, 'i*.nii'),...
            fullfile(segmentedDir, subNames{ss}, 'anat'));
        if ~success
                warning(message)
        end
        [success,message] = movefile(fullfile(rawSubDir, anatDir, 'y*.nii'),...
            fullfile(segmentedDir, subNames{ss}, 'anat'));
        if ~success
                warning(message)
        end
        [success,message] = movefile(fullfile(rawSubDir, anatDir, '*.mat'),...
            fullfile(segmentedDir, subNames{ss}, 'anat'));
        if ~success
                warning(message)
        end
        
        clearvars matlabbatch
    end
    
    %% Normalisation: write only (uses deformations field from segmentation)
    if do.normalisation
        
        fprintf('NORMALISATION USING SEGMENTATION\n\n')
        
        anatDeformationImg  = dir(fullfile (segmentedDir, subNames{ss}, 'anat','y_*.nii'));
        folderContent       = dir(fullfile(coregisteredDir, subNames{ss}, 'func', ['au' subNames{ss} '_' taskName '_' 'run*.nii']));
        nruns               = length(folderContent);
        
        meanEpi = dir(fullfile(coregisteredDir, subNames{ss}, 'func', ['meanu' '*.nii']));
        meanEpi = fullfile(coregisteredDir, subNames{ss}, 'func', meanEpi.name);
        
        alltargets = {}; %JB It's necessary to initialize alltargets as a cell array, as this can prevent the vertcat error if paths from different runs have different character lengths.
        for r = 1:nruns
            sessionDir = fullfile(coregisteredDir, subNames{ss}, 'func');
            dirfiles     = spm_select('ExtFPList', sessionDir, folderContent(r).name, Inf);
            if strcmp(dirfiles,'')
                warning('No files selected!');
                return;
            end
            alltargets = [alltargets; cellstr(dirfiles)]; 
        end
        
        alltargets = cellstr(alltargets);
        alltargets{end+1} = meanEpi;
        alltargets{end+1} = rawSubAnatImg; % add anatomical image to normalise
        
        matlabbatch{1}.spm.spatial.normalise.write.subj.def         = {fullfile(anatDeformationImg.folder, anatDeformationImg.name)};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample    = alltargets;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb      = [-78 -112 -70
            78  76   85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox     = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp  = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix  = 'w';
        
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clearvars matlabbatch
        
        % move created normalized niftis in a seperate folder
        [success,message] = movefile(fullfile(coregisteredDir, subNames{ss}, 'func', 'w*.nii'),...
            fullfile(normalizedDir, subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
        
        % move the created normalized anatomical image in a seperate folder
        [success,message] = movefile(fullfile(rawSubDir, anatDir, 'w*.nii'),...
            fullfile(normalizedDir, subNames{ss}, 'anat'));
        if ~success
                warning(message)
        end
    end
    
    %% Smoothing
    if do.smoothing
        fprintf('SMOOTHING\n\n')
        if strcmp(do.smoothNorm,'mni') == 1
                currDir = normalizedDir;
                folderContent = dir(fullfile(currDir, subNames{ss}, 'func', ['wau' subNames{ss} '_' taskName '_' 'run*.nii']));
            elseif strcmp(do.smoothNorm, 'native') == 1
                currDir = coregisteredDir;
                folderContent = dir(fullfile(currDir, subNames{ss}, 'func', ['au' subNames{ss} '_' taskName '_' 'run*.nii']));
            elseif strcmp(do.smoothNorm, 'both') == 1
                 % TO-DO
        end
        
        nruns         = length(folderContent);
        
        alltargets = {};
                
        for r = 1:nruns
            
            if strcmp(do.smoothNorm,'mni') == 1
                sessionDir = fullfile(normalizedDir, subNames{ss}, 'func');
            elseif strcmp(do.smoothNorm, 'native') == 1
                sessionDir = fullfile(coregisteredDir, subNames{ss}, 'func');
            elseif strcmp(do.smoothNorm, 'both') == 1
                 % TO-DO
            end
            dirfiles     = spm_select('ExtFPList', sessionDir, folderContent(r).name, Inf);
                
            if strcmp(dirfiles,'')
                warning('No files selected!');
                return;
            end
            alltargets = [alltargets; cellstr(dirfiles)];
        end
        
        matlabbatch{1}.spm.spatial.smooth.data      = cellstr(alltargets);
        matlabbatch{1}.spm.spatial.smooth.fwhm      = repmat(do.smoothingSize,1,3); % Set at the beginning.
        matlabbatch{1}.spm.spatial.smooth.dtype     = 0;
        matlabbatch{1}.spm.spatial.smooth.im        = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix    = ['s' num2str(do.smoothingSize)];
        
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clearvars matlabbatch
        [success,message] = movefile(fullfile(currDir, subNames{ss}, 'func', ['s' num2str(do.smoothingSize) '*.nii']),...
            fullfile(smoothedDir, subNames{ss}, 'func'));
        if ~success
            warning(message)
        end
    end
end
end