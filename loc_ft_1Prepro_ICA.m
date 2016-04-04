%% A complete fieldtrip script to run preprocessing on localiser
% The 1st step in FT localiser analysis.
% This script will run the following steps:
% 1. Preprocessing
% 2. Timelock analysis (ERP) -- for later comparison with ICA cleaned data
% 3. ICA (reject component 1)
% 4. Timelock analysis (ERP) again, plot topography for N1.
% The next step is loc_ft_2PlotERPandTopo.
% _____ 2016-03-13 Sijia Zhao; Edited on 2016-04-04 for demostration

sublist = [1]; % Please adjust for yours
root = '160404_demonstration/'; if exist(root)==0; mkdir(root); end % Please adjust for yours
filename_tail='_loc';

for subi = sublist % Loop over all subjects
    
    % Your data filename
    filename = ['Datasets/sub' num2str(subi) '/eeg2_sub' num2str(subi) filename_tail '.bdf'];
    
    %% Step 1: Define trial (has to be done before ft_preprocessing
    hdr = ft_read_header(filename);
    cfg = [];
    cfg.dataset =  filename; % your filename with file extension;
    cfg.trialdef.eventtype  = 'STATUS'; % Status notation maybe Biosemi's  tigger name
    cfg.trialdef.eventvalue = 61; % your event value
    cfg.trialdef.prestim    = 0.1;  % before stimulation (sec), only use positive value
    cfg.trialdef.poststim   = 0.25; % after stimulation (sec) , only use positive value
    cfg = ft_definetrial(cfg);
    
    %% Step 2: Jump Removal -- before all filters; required before ICA
    cfg.continuous = 'yes';
    cfg.artfctdef.jump.channel='EEG';
    [cfg, artifact] = ft_artifact_jump(cfg);
    
    %% Step 3: Preprocessing
    cfg.hpfilter ='yes';
    cfg.hpfreq =0.1;
    cfg.hpfiltdir = 'twopass';
    cfg.hpfiltord = 4;
    
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 30;
    cfg.lpfiltdir = 'twopass';
    cfg.lpfiltord = 4;
    
    cfg.detrend ='yes';
    cfg.channel = {'EEG' 'EOG'};
    cfg.continuous = 'yes';
    
    cfg.reref = 'yes';
    cfg.refchannel = 'EEG';
    raw_data  = ft_preprocessing(cfg);
    
    %% Step 4: Downsample to 256Hz
    cfg = [];
    cfg.resamplefs = 256;
    data = ft_resampledata(cfg,raw_data);
    
    %% Step 5: ERP Analysis
    cfg = [];
    cfg.channel = 'EEG';
    cfg.trials = 'all';
    cfg.covariance         = 'no';
    cfg.covariancewindow   = 'all';
    cfg.keeptrials         = 'yes';
    cfg.removemean         = 'no';
    cfg.vartrllength       = 0;
    timelock = ft_timelockanalysis(cfg,data);
    
    % Baseline correction
    cfg.baseline = [-0.1 0];
    cfg.channel = 1:64;
    cfg.parameter ='avg';
    timelock = ft_timelockbaseline(cfg,timelock);
    
    % Grand average
    gradavg = ft_timelockgrandaverage(cfg,timelock);
    
    % Plot topography for N1
    load('Inputs\Electrodes_labels_Biosemi64.mat');
    gradavg.label = Electrodes_labels_Biosemi64';
    cfg =[];
    cfg.xlim = [0.08 0.12];
    cfg.layout = 'biosemi64.lay';
    figure; ft_topoplotER(cfg,gradavg); colorbar;
    
    % Calculate RMS (and plot)
    rms_individual_beforeICA = sqrt(nanmean((gradavg.avg).^2,1));
    
    % % Step: ft_timelockstatistics
    %     stat = ft_timelockstatistics(cfg, gradavg);
    
    
    %% Run ICA to remove eye blinks
    cfg = [];
    cfg.channel = 'EEG';
    cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
    comp = ft_componentanalysis(cfg, data); % cdConclpc is the is the structure containing your data (probably called S in your case?)
    
    % view ICA across channels
    load('labelname');
    comp.topolabel = labelname;
    
    cfg = [];
    cfg.continuous = 'yes';
    %     cfg.blocksize = 30;
    cfg.component = 1:10;
    cfg.viewmode = 'component';
    cfg.layout = 'biosemi64.lay'; % specify the layout file that should be used for plotting
    ft_databrowser(cfg, comp);
    %
    
    % % Plot topographies for the components
    %     cfg =[];
    %     cfg.component =1:10;
    %     cfg.layout = 'biosemi64.lay';
    %     cfg.comment = 'no';
    %     ft_topoplotIC(cfg,comp);
    
    % Reject ICA components
    cfg = [];
    cfg.component = [1];
    data_icacleaned = ft_rejectcomponent(cfg, comp); %Always error: Reference to non-existent field 'balance'. Error in ft_rejectcomponent (line 217) if ~isfield(sens.balance, 'invcomp')
    
    %% ERP analysis again (after ICA)
    cfg = [];
    cfg.channel = 1:64;
    cfg.trials = 'all';
    cfg.covariance         = 'no';
    cfg.covariancewindow   = 'all';
    cfg.keeptrials         = 'yes';
    cfg.removemean         = 'no';
    cfg.vartrllength       = 0;
    timelock = ft_timelockanalysis(cfg,data_icacleaned);
    
    cfg.baseline = [-0.1 0];
    cfg.channel = 1:64;
    cfg.parameter ='avg';
    timelock = ft_timelockbaseline(cfg,timelock);
    
    gradavg = ft_timelockgrandaverage(cfg,timelock);
    rms_individual = sqrt(nanmean((gradavg.avg).^2,1));
    
    %% Save and plot
    filename_save = [root 'cleaneddata_loc_sub' num2str(subi) '.mat'];
    save(filename_save,...
        'data_icacleaned','gradavg','rms_individual','rms_individual_beforeICA');
    
    toplot =1;
    if toplot
        load('Inputs\Electrodes_labels_Biosemi64.mat');
        gradavg.label = Electrodes_labels_Biosemi64';
        cfg =[];
        cfg.xlim = [0.08 0.12];
        cfg.layout = 'biosemi64.lay';
        figure(subi); ft_topoplotER(cfg,gradavg); colorbar;
        
        figure(subi*10); plot(rms_individual_beforeICA,'-k');hold on; plot(rms_individual,'-r'); hold off
        legend('beforeICA','after ICA');
    end
    %     figure(1); for i=1:185; plot(raw_data.time{i},mean(raw_data.trial{i})); hold on;end
    %     figure(2); for i=1:185; plot(data_icacleaned.time{i},mean(data_icacleaned.trial{i})); hold on;end
end