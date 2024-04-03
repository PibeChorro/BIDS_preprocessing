# TODO-List 
## Minor changes

### All files

- [ ] Ask whether or not you are working on the cluster or not and adjust the file selection.<br>
Right now the selection of files is managed with a gui. Does not work on the cluster over the terminal
- [ ] Include an overwrite option that skips steps if they were performed already
- [ ] Ask user via command line which steps he/she would like to perform

### `sortDicomsIntoFolders.m`

- [ ] Change the method to get the scan number. 
  - Right now it is read out of the image name -> Seperating the image name and take the 4th entry
    - pro: it is fast
    - con: it is not really reliable
  - Alternative is to read in the header of the DICOM and read out the field `SeriesNumber`
    - pro: reliable
    - con: very slow
    
## Major changes

### `sortDicomsIntoFolders.`

- [ ] MAYBE merge the steps of moving files in folders corresponding to scan number and of moving those to sequence type
- [ ] add sorting for different conditions in functional data folder
- [ ] add subfolders for different anatomical scans

### `DICOMconversionSPM12_BIDS.m

- [ ] get the json-file done
  - [x] check for entries that are MUST HAVES and include them
  - [ ] add a `switch` `case` statement to ask for different scanner/institutes, hence differing DICOM headers
