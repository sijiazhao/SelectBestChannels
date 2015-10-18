# SelectBestChannels
Using localiser's data, this SPM script is used to select best channels in the butterfly plots.

1. Use SZ_SPMPreprocessing_loc.m to preprocess the EEG localiser data first.(Extra functions are required, you can find on our lab EEG analysis_PC.
2. then run SZ_run_find5BestChanns.m to use SZ_fSelectBestChannels.m to find best channels and plot the butterfly plots for each individual
subjects.
3. to show the channels (selected from localiser) in the butterfly of real experiment, run SZ_fSelectBestChannels again 
(notice: you have to have already done step 2 for the localiser, to get the file BestChanns_loc.mat)