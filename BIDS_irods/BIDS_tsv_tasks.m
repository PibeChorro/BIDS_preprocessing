%% Create .tsv file for your functional runs with task name and run for each subject
subjects = 6;

hdr  = {'Run' 'Task'}; % header
runs = {'1','2','3','4','5','6'}'; % which functional runs
tasks = {'TMSlow','TMShigh'}; % which task

ss = 1; % subject 1: 
vals(ss).hdr = hdr;
vals(ss).runs = runs;
vals(ss).tasks = {tasks{1},tasks{2},tasks{1},tasks{2},tasks{1},tasks{2}}; % Low-high-low-high-low-high


for i = 1:length(subjects)

    currDir = fullfile(sourceDir,subNames{ss}); % dicoms
    
    s = subjects(i);

    tmpfileName = ['sub-' num2str(s,formatSpec) '_scans'];

    tbl = cell2table(cat(2,vals(s).runs,vals(s).tasks));
    tbl.Properties.VariableNames = vals(s).hdr; 
    writetable(tbl,tmpfileName,'Delimiter','\t');
    movefile([tmpfileName '.txt'],[tmpfileName '.tsv'],'f');
    
end