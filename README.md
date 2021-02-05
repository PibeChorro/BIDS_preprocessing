# BIDS_preprocessing
A matlab prepocessing pipeline to create a BIDS compatible data set <br>
0: You need to have a project folder containing a "source_data" folder. This folder should contain your subject folders which themselves only contain the DICOMs<br>
- myProject
  - source_data
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
1c-extra: If a sequence is called 'anat' or 'struct' the folder containing these images will be **RENAMED** to 'anat' **NOT MOVED** into a 'anat/struct' folder and named 'run-XX'<br>
1d: It saves the delivered sequence information, names and coresponding number of images

- myProject
  - source_data
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
2a: It askes you to select your project folder (the one CONTAINING your source_data folder)<br>
2b: The script creates a 'raw_data' folder in your project folder. This one is the folder you really want to be BIDS conform.<br>
2c: Your DICOMs in the 'source_data' folder will be converted into NIfTIs and stored in the 'raw_data' folder (only the 'non-error' data folders are included)

- myProject
  - raw_data
    - sub01
      - func
        - run-01
        - run-02
        - ...
      - anat
    - sub02
    - sub03
  - source_data

**Before you can run `preprocessingSPM12_BIDS`you need to adjust the location of SPMs TMP.nii file**<br>
3: Run `preprocessingSPM12_BIDS.m` (you may previously want to specify preprocessing step you want to perform and which smoothing size you want to use) The script does the following:<br>
3a: It askes you to select your project folder (the one CONTAINING your source_data and raw_data folder).<br>
3b: It creates a 'derivative_data' folder<br>
3c: For every subject it creates an 'MRI' folder and within a 'func' and an 'anat' folder.
3d: For every preprocessing step it creates a corresponding folder and stores the resulting NIfTIs in their.

- myProject
  - derivative_data
    - sub01
      - MRI
        - anat
          - normalized
          - segmented
        - func
          - 6smoothed
          - coregistered
          - normalized
          - realigned
          -slice_time_corrected
    - sub02
    - sub03
  - raw_data
  - source_data

4: Move all your additional data, such as behavioural-, eyetracking- or questionnaire data into the derivative folder

- myProject
  - derivative_data
    - sub01
      - eyetracking
      - MRI
      - psychoPhysic
    - sub02
    - sub03
  - raw_data
  - source_data
