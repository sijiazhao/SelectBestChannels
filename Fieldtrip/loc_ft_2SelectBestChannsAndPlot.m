%% The 3rd function in FT localiser analysis
% Run this script after loc_ft_1Prepro_ICA (or loc_ft_2PlotERPandTopo)
% Requires
% * pre-saved data: [root 'cleaneddata_loc_sub' num2str(subi) '.mat']
% * (Optional) targettime: the peak time roughly measured by ERP RMS. You
% can specify for each subject (in loc_ft_2PlotERPandTopo), or just use default N100 time, 0.1s
% For further experimental analysis, please use the new output 'BestChanns_10_posneg.mat'
% This is the end of localiser analysis!
% ______________ Sijia Zhao (R) 2016-03-14; Edited on 2016-04-04 for demostration
close all; clear;clc;
sublist = [1]; % Please adjust for yours
root = '160404_demonstration/';

numChannsToBeSelected = 10;
selectionstyle = 'posneg';

%% Find Best channels and plot butterfly
% %%%%%%%%%% targettime 
% In loc_ft_2PlotERPandTopo, you can manually measure the peak time by looking at ERP, and then type them in here.
% Or, you can use the default value targettime = 0.1(s), by 'targettime_rough = []'
% At the end of this step, a more accurate target time will be saved targettime_accurate.
targettime_rough = [];
targettime_accurate = []; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeinterval_searching =0.01; % in seconds
BestChanns_allsubs =[];

figure(1);
for subi = sublist % Loop over all subjects
    load([root 'cleaneddata_loc_sub' num2str(subi) '.mat']); %%Set at the end of loc_ft_1Prepro_ICA
    % %%%%%%%%% Notes: Data used in this script (after loc_ft_1Prepro_ICA)
    %*gradavg: gradavg = ft_timelockgrandaverage(cfg,timelock); <--- data preprocessed and grand averaged; it has to include gradavg.time and gradavg.avg
    %*rms_individual: rms_individual = sqrt(nanmean((gradavg.avg).^2,1));
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if numel(sublist)>1; subplot(round(length(sublist)/5+1),5,subi); end
    
    % Step 1: Find the target time for the peak (N100 or ?)
    if isempty(targettime_rough); targettime = 0.1; else targettime=targettime_rough(subi); end
    % Covert time unit s into the column number
    time_xaxis = gradavg.time; % N*1 array, time in seconds %<------ From loaded file
    % To find the target time in the N*1 array (i.e. find the closest value
    tmp = abs(time_xaxis-targettime); [idx idx] = min(tmp); closest = gradavg.time(idx); %idx is the index of the closest value of target time %<------ From loaded file
    timeinterval_numIdx = round(timeinterval_searching/(time_xaxis(2)-time_xaxis(1)));
    N100_RMS = rms_individual(idx-timeinterval_numIdx : idx+timeinterval_numIdx);
    [maxNum, maxIndex] = max(rms_individual(idx-timeinterval_numIdx:idx+timeinterval_numIdx));
    [idx_max1,idx_max2] =find(rms_individual==maxNum);
    targettime=time_xaxis(idx_max1,idx_max2);
    targettime_accurate= [targettime_accurate targettime]; %Extra$: more accurate targettime. Save newly found target time for topography
    
    % Step 2: Find the best channels at the peak
    BestChanns=[];
    temp = gradavg.avg(:,idx_max2); %all 64 channels' values
    for iChannel = 1:numChannsToBeSelected
        if strcmp(selectionstyle,'posneg')
            if iChannel > numChannsToBeSelected/2
                [maxNum, maxIndex] = max(temp);  %         if absolute ==1; % else [maxNum, maxIndex] = max(PosOrNeg*(temp(:))); end
            else
                [maxNum, maxIndex] = min(temp);
            end
        elseif strcmp(selectionstyle,'abs')
            [maxNum, maxIndex] = max(abs(temp(:))); % all absolute
        elseif strcmp(selectionstyle,'neg')
            [maxNum, maxIndex] = max(-(temp(:))); % all negative
        end
        BestChanns = [BestChanns maxIndex];
        temp(maxIndex,:) = 0;
    end
    
    %% Step 3: Plot butterfly for all, (bad) and BEST channels
    % Plot the searching-best-channel time interval as pink
    ylim_value = 5;
    pink =[255 182 193]/255;
    v1=[targettime-timeinterval_searching targettime+timeinterval_searching targettime+timeinterval_searching targettime-timeinterval_searching];
    v2 = [-ylim_value -ylim_value ylim_value ylim_value];
    patch(v1,v2,pink,'FaceAlpha', .3, 'EdgeColor','none'); hold on
    
    % plot butterfly
    allchanns = gradavg.avg;
    for i=1:size(allchanns,1); plot(time_xaxis,allchanns(i,:),'-','Color',[.5 .5 .5]);hold on; end
    % Badchannels
    % Best channels
    for i=BestChanns; plot(time_xaxis,allchanns(i,:),'-r');hold on;end
    
    % plot RMS ERP
    plot(time_xaxis,rms_individual,'-k','LineWidth',3); hold off;
    
    xlim([-0.1 0.25]);
    ylim([-ylim_value ylim_value]);
    title(num2str(subi));
    
    % Step4: save the bestchannels
    BestChanns_allsubs = [BestChanns_allsubs; BestChanns];
end
loc_suffix=['_' num2str(numChannsToBeSelected) '_' selectionstyle];
save(['BestChanns' loc_suffix],'BestChanns_allsubs','targettime_accurate');
save([root 'BestChanns' loc_suffix],'BestChanns_allsubs','targettime_accurate');

%% Step 4: Plot topography (with updated target time)
if ~isempty(targettime_accurate)
    figure(2);
    load('Inputs\Electrodes_labels_Biosemi64.mat');
    for subi = sublist % Loop over all subjects
        load([root 'cleaneddata_loc_sub' num2str(subi) '.mat']);
        if numel(sublist)>1; subplot(round(length(sublist)/5+1),5,subi); end
        cfg = [];
        cfg.layout = 'biosemi64.lay';
        layout = ft_prepare_layout(cfg,gradavg);
        gradavg.label = layout.label(1:size(gradavg.avg,1));
      
        cfg = [];
        cfg.parameter = 'avg';
        cfg.channel = 'all';
        cfg.layout =layout;
        cfg.xlim = [targettime_accurate(subi)-0.01 targettime_accurate(subi)+0.01]; %     cfg.xlim = [0.08 0.12];
        cfg.zlim            = 'maxabs'; %[-3.5 2.5]*10^-13;
        cfg.marker          = 'on';
        cfg.interactive     = 'no';
        cfg.comment         = 'no';
        cfg.colorbar        = 'no';
        cfg.highlight       = 'on';
        cfg.highlightchannel = gradavg.label(BestChanns_allsubs(subi,:));
        cfg.highlightcolor  = 'k';
        cfg.highlightsymbol = '.';
        cfg.highlightsize   = 20;
        cfg.markersymbol    = '.';
        cfg.markersize      = 3;
       
        ft_topoplotER(cfg,gradavg);
        title(num2str(subi));       
    end
end