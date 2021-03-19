# BIDS_preprocessing
A matlab prepocessing pipeline to create a BIDS compatible data set <br>
_________________________________________________________________________________________________________________________________________________________
**IMPORTANT!!!** All the scripts here were not yet tested on a REAL BIDS-dataset and for now only works on a very simple dataset.<br>
It still has to be tested for datasets containing additional `tsv` and `json` files.<br>
_________________________________________________________________________________________________________________________________________________________

0: You need to have a project folder containing a "sourcedata" folder. This folder should contain your subject folders which themselves only contain the DICOMs<br>
- myProject
  - sourcedata
    - sub01
      - SUB001.....IMA
      - ...
      - ...
    - sub02
    - sub03

1: Run `sortDicomsIntoFolders.m` and select your "source_data" folder. The script does the following:<br>
1a: For every subject it sorts all DICOMs into folders named after the sequencenumber in which they were generated (i.e. 01 02 ...)<br>
1b: It then enters every created folder, reads out the sequence information from the first DICOM header and if this sequence is unknown (for the very first subject it will be unknown), it asks you to give it a proper name (functional data = 'func', anatomical = anat, diffusion wheight imaging = dwi, the rest you can call as you like).
You also need to deliver the number of images in this scan. ***!!!IMPORTANT!!!*** Make sure you deliver the right number of images. <br>
1c: Every folder is then moved to its corresponding sequence folder and named run-XX (if a folder contains more or less images than previously specified it will be saved as error-run-XX)<br>
1c-extra: If a sequence is called 'anat' or 'struct' you are asked to assign the modality of your anatomical scan (T1, T2 or T2star). The folder is then **RENAMED** to 'anat/<modality>' **NOT MOVED** into a 'anat/struct' folder and named 'run-XX'<br>
1d: It saves the delivered sequence information, names and coresponding number of images

- myProject
  - sourcedata
    - sequenceInfo.mat
    - sub01
      - anat
      - func
        - run-01
        - error-run-02
        - run-02
        - ...
      - loc
      - shim
    - sub02
    - sub03

2: Run `DICOMconversionSPM12_BIDS.m` (you may previously want to specify which data type you want to convert - for now only functional and anatomical data) The script does the following:<br>
2a: It askes you to select your project folder (the one CONTAINING your sourcedata folder) and to specify your subject prefix for the DICOMS.<br>
2b: The script creates a 'rawdata' folder in your project folder. This one is the folder you really want to be BIDS conform.<br>
2c: Your DICOMs in the 'sourcedata' folder will be converted into NIfTIs and stored in the 'rawdata' folder (only the 'non-error' data folders are included).<br>
2d: The created NIfTIs are then compressed into 4D NIfTIs and the existing 3D NIfTIs are deleted.<br>
2di: The 4D NIfTI files are named BIDS conform. To this point it only works for one session and one task.<br>
2e: For functional and anatomical scans a .json sidecar file is created. For now only one .json file is created per subject (we assume the same parameters for all runs).<br>

- myProject
  - rawdata
    - sub-01
      - func
        - sub-01_task-\<label\>.json
        - sub-01_task-\<label\>_run-\<index\>_bold.nii
        - sub-01_task-\<label\>_run-\<index\>_bold.nii
        - ...
      - anat
        - sub-01_\<modality\>.json
        - sub-01_\<modality\>.nii
    - sub-02
    - sub-03
  - sourcedata

**Before you can run `preprocessingSPM12_BIDS`you need to adjust the location of SPMs TMP.nii file**<br>
3: Run `preprocessingSPM12_BIDS.m` (you may previously want to specify preprocessing step you want to perform and which smoothing size you want to use) The script does the following:<br>
3a: It askes you to select your project folder (the one CONTAINING your sourcedata and rawdata folder).<br>
3b: It creates a 'derivatives' folder.<br>
3c: For every subject it creates an 'MRI' folder and within a 'func' and an 'anat' folder.
3d: For every preprocessing step it creates a corresponding folder and stores the resulting NIfTIs in their.

- myProject
  - derivative_data
    - \<other-pipeline\>
    - spm12
      - spm12-preproc
        - 6smoothed
          - sub-01
            - func
              - s6wausub-01_task\<label\>_run-\<index\>_bold.nii
              - ...
          - sub-02 
        - coregistered
        - normalized
        - realigned
        - segmented
        - slice_time_corrected
  - rawdata
  - sourcedata
