function json_file = BIDS_select_sessions(rawdir)

        disp('===========================');
        disp('CREATING FMRI-PREP FILTER FILE');
        disp('===========================');
        BIDs_filter.bold.datatype = "func";
        BIDs_filter.bold.session  = "02";
        BIDs_filter.bold.suffix   = "bold";
        jsonwrite([params.rawDir '/fmriprep_BIDS_filter_ses2.json'], BIDs_filter, 'prettyPrint', true);
end
