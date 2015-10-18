% - This script uses function SZ_fSelectBestChannels.m
% - save BestChanns_loc.mat for SZ_fSelectBestChannels(subi,0,condi) and
% further analysis in real experiments.
%________________________________Written by Sijia Zhao (UCL), 17-10-2015
close all; clear; clc;

BestChanns_allsubs = []; N100Time_allsubs = []; 
sublist  = 1:20; % subjects to include
loc = 1; condi = 1; 
for subi = sublist
    %     SZ_SPMPreprocessing_loc(subi);
    [BestChanns, targetTime] = SZ_fSelectBestChannels(subi,loc,condi);
    BestChanns_allsubs = [BestChanns_allsubs; BestChanns];
    N100Time_allsubs = [N100Time_allsubs targetTime];
end

save('BestChanns_loc','BestChanns_allsubs','N100Time_allsubs'); % BestChanns_loc.mat 