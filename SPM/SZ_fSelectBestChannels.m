function [BestChanns, targetTime] = SZ_fSelectBestChannels(subi,loc,condi)
% This function finds the best channels in the buttterfly plots of
% localiser.
% - use this after preprocessing (separate code)
% * Best channels: channs will largest (absolute) amplitude at N100
% (auditory reaction) 
% %%%% INPUTS: (Please modify these to adapt your experiments)
% subi = subject index
% loc = whether is localiser or not? (always should run localiser first)
% condi = which condition?
% %%%% OUTPUTS:
% BestChanns: list of selected best channels
% targetime: time when M100 happened(determined in RMS, then used to find best channels)
%________________________________Written by Sijia Zhao (UCL), 17-10-2015
numChannsToSelect = 5;
saveplots=1; % 1 to save butterfly plots

if loc == 1
    directory  = './';
    data = [directory 'm' 'subject' num2str(subi) '_loc']; %output file from SZ_SPMPreprocessing_loc.m
    condlist = {'click'}; % list of condition
    filename = ['./Plots/butterfly_channs_sub' num2str(subi) '_loc'];
    stimLength  = .15; bleedin = .1; bleedout = .1; %xtime = [-.1 .25]
elseif loc == 0
    directory = './';
    data = [directory 'subject' num2str(subi)];
    condlist = {'H','M','E','N'};
    filename = ['./Plots/butterfly_channs_sub' num2str(subi)];
    stimLength  = 4; bleedin = .5; bleedout = 1; %xtime = [-.5 5];
end
startPlot = -bleedin; endPlot = stimLength + bleedout; %xlim for plot, in s

if nargin < 2; disp('Too few arguments'); return; end
if nargin < 3; condi = 1; end % if did not specify which condition to choose, choose the first in condlist
if nargin == 2; if (loc ~=0) || (loc ~= 1); disp('Please specify the data. localiser or real expeirments?'); return; end; end

spm('defaults', 'eeg');

%% Manually remove outlier channels (if applicable)
badchanns =[];
if subi == 3; badchanns =[1,2,33,34]; badchanns_twindow_in = -.1; badchanns_twindow_out = -.08; end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load EEG data
D      = spm_eeg_load(data);
subset = D.chanlabels(indchantype(D,'EEG')); % select all channels

g = figure(subi);

%% In localiser, find when N100 exactly happens, in RMS, range 100±30 ms
temp = selectdata(D, subset, [startPlot endPlot], condlist{1,condi});

if ~isempty(badchanns); temp(badchanns,:) = 0; end % if bad chann, the data wont be considered in selection.

% Compute RMS and plot RMS
rms_individual = sqrt(nanmean(temp.^2,1));
plot(linspace(startPlot,endPlot,length(rms_individual)),rms_individual,'-k','LineWidth',1); hold on

if loc == 1;
    N100_time = 100; N100_bleed = 30; % preset a range to find, in ms
    N100_startTime = N100_time-N100_bleed; N100_endTime = N100_time+N100_bleed;
    % In spm, the raw data saved in a matrix where each column = 3.75ms
    N100_startCol = round((N100_startTime-startPlot*1000)/3.75); N100_endCol = round((N100_endTime-startPlot*1000)/3.75);
    % N100_timeWindow = [N100_startTime, N100_endTime]/1000;
    N100_rms_window = rms_individual(N100_startCol:N100_endCol);
    
    [maxNum, maxIndex] = max(N100_rms_window);
    [row, col] = ind2sub(size(rms_individual), maxIndex);
    targetTime = (N100_startCol+ col)*3.75 + startPlot*1000; %convert the time in matrix column number to the time in ms
    % N100channs = row;

    BestChanns = [];
end
if loc == 0; load('BestChanns_loc_maria.mat'); targetTime = N100Time_allsubs(subi); BestChanns = BestChanns_allsubs(subi,:); end

searchwindow_bleed = 15;

% Update time window to the new target time (N100's peak time point)
searchwindow_start = (targetTime-searchwindow_bleed)/1000;  searchwindow_end = (targetTime+searchwindow_bleed)/1000;
searchwindow = [searchwindow_start searchwindow_end];

% ------------- Find the channels with max values in the time window given --------
temp        = selectdata(D, subset, searchwindow , condlist{1,condi});
% if ~isempty(badchanns); temp(badchanns,:) = 0; end

% Main part: Find the channels with max values in the time window given
% BestChanns = [];
if isempty(BestChanns) %loc == 1    
    for iChannel = 1:numChannsToSelect
        [maxNum, maxIndex] = max(abs(temp(:)));
        [row, col] = ind2sub(size(temp), maxIndex);
        BestChanns = [BestChanns row];
        temp(row,:) = 0;
    end
end

%% Butterfly plot: all channs (black), badchanns (blue), selected good channs (red)

% Shaded area: searchwindow (pink)
ylim_value = 5;%please specify this by yourself
v1= [searchwindow_start searchwindow_end searchwindow_end searchwindow_start];
v2= [-ylim_value -ylim_value ylim_value ylim_value];
patch(v1,v2,[255 182 193]/255,'FaceAlpha', .3, 'EdgeColor','none'); hold on %pink: [255 182 193]

% Shaded area: badchannels search window (blue) (pre-set)
if ~isempty(badchanns)
    v1=[badchanns_twindow_in badchanns_twindow_out badchanns_twindow_out badchanns_twindow_in];
    patch(v1,v2,[49 140 231]/255,'FaceAlpha', .2, 'EdgeColor','none'); hold on % bleu de france: [49 140 231]
end

% Plot best channs(red) and bad channs(blue)
if loc == 1
    plot(D.time, squeeze(D(indchantype(D,'EEG'),:,1)),'-','Color',[.7 .7 .7]); hold on
    if ~isempty(badchanns)
        plot(D.time, squeeze(D(badchanns,:,1)),'-b'); hold on
    end
    
    plot(D.time, squeeze(D(BestChanns,:,1)),'-r'); hold on % plot best channels (red)
    % plot(D.time, squeeze(D(fixedchanns,:,1)),':g'); hold off
    plot(linspace(startPlot,endPlot,length(rms_individual)),rms_individual,'-k','LineWidth',3); hold off % plot 
end

if loc == 0
    subset = D.chanlabels(indchantype(D,'EEG')); % select all channels
    temp = selectdata(D, subset, [startPlot endPlot], condlist{1,condi});
    plot(linspace(startPlot,endPlot,length(temp)), temp,'-','Color',[.7 .7 .7]); hold on
    
    if ~isempty(badchanns)
        BadChannsSet = D.chanlabels(badchanns);
        BadChannsSelection = selectdata(D,BadChannsSet,[-.5 5], condlist{1,condi});
        plot(linspace(startPlot,endPlot,size(BadChannsSelection,2)), BadChannsSelection,'-b'); hold on
    end
    BestChannsSet = D.chanlabels(BestChanns);
    BestChannsSelection = selectdata(D,BestChannsSet,[-.5 5], condlist{1,condi});
    plot(linspace(startPlot,endPlot,size(BestChannsSelection,2)), BestChannsSelection,'-r'); hold on
    
    rms_individual = sqrt(nanmean(temp.^2,1));
    plot(linspace(startPlot,endPlot,length(rms_individual)),rms_individual,'-k','LineWidth',3); hold off
end

xlim([startPlot endPlot]);

if ~isempty(badchanns)
    h = legend('(Black) RMS of ALL channels', ['(Red) Best channs: ' num2str(BestChanns)], ['(Blue) Bad channs: ' num2str(badchanns)]);
else
    h = legend('(Black) RMS of ALL channels',['(Red) Best channs: ' num2str(BestChanns)]);
end
h.Location = 'southwest';

if loc == 0
    title(['\rm Butterfly and RMS plot for subject ' num2str(subi) ', condition ' condlist{condi}]);
else
    title(['\rm Butterfly and RMS plot for subject ' num2str(subi)]);
end

% Save plots
if saveplots == 1
    if loc == 0; filename = [filename '_' condlist{condi}]; end
    filenamePNG = [filename '.png'];
    print(g,filenamePNG,'-dpng');
end
end