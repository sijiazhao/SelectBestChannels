function SZ_SPMPreprocessing_loc(subi)
%% EEG preprocessing code in SPM
%________________________________Written by Sijia Zhao (UCL), 17-10-2015
spm('defaults', 'eeg');

root  = './';
outdirraw = [root 'Datasets/'];
outdirdata  = [root 'SPMdata/'];

exptName = 'RAND20TO10'; %for filenaming
sublist  = 1:20;

condlist = {'click'};
triglist = [220]; triglist = round(triglist*2048/44100);
stimulusLength  = 150; %in ms

timewindow_td       = [-100 stimulusLength+100]; %trial definition
timewindow_epoch    = [-100 stimulusLength+100];
timewindow_bc       = [-100 0];

jumps           = 0; if jumps == 1; prefix_jumps = 'j';else prefix_jumps = ''; end
eyeblink        = 0; if eyeblink == 1; prefix_eyelink = 'T'; else prefix_eyelink = ''; end
baselinec       = 0; if baselinec == 1; prefix_baselinec = 'b'; else prefix_baselinec = ''; end % Dont do BC twice (first in epoch)

outfile_folder =['sub' num2str(subi) '/'];
outfile = [exptName '_sub' num2str(subi) '_loc'];

%% Convert .bdf file into a .mat file and a .dat file
S                 = [];
S.dataset         = [outdirraw outfile_folder outfile '.bdf'];
S.mode            = 'continuous';
S.outfile         = [outdirdata outfile];
S.eventpadding    = 0;
S.blocksize       = 3276800;
S.checkboundary   = 1;
S.saveorigheader  = 0;
S.timewin         = [];
S.conditionlabels = {'Undefined'};
S.inputformat     = [];
D = spm_eeg_convert(S);
disp(['Converted subject ' num2str(subi) '/' num2str(max(sublist))]);
%% Prepare: define channel types
% EOG: EXG1 & EXG2
S               = [];
D               = spm_eeg_load([outdirdata outfile]);
S.D             = D;
S.task          = 'settype';
S.type          = 'EOG';
S.ind           = indchannel(D,{'EXG1', 'EXG2'}); % get indices of chans EXG1 and EXG2
S.save          = 1;
D = spm_eeg_prep(S);
disp(['Defined channels for subject ' num2str(subi) '/' num2str(max(sublist)) ]);

% (Optional) Mastoids EEG:EXG3 & EXG4 - Needed for mastoids reference or to
%use mastoids as extra temporal electrodes
%         S      = [];
%         S.D    = [outdirdata outfile '.mat']; % SS
%         S.task = 'settype';
%         S.type = 'EEG';
%         S.ind  = [67 68];
%         S.save = 1;
%         D = spm_eeg_prep(S);
%
%% High-pass filter (0.1 Hz) (Must before everything!)
S        = [];
S.D      = [outdirdata outfile];
S.type   = 'butterworth';
S.band   = 'high';
S.freq   = .1;
S.dir    = 'twopass';
S.order  = 5;
S.prefix = 'f';
D = spm_eeg_filter(S);
disp(['High-pass filtered subject ' num2str(subi) '/' num2str(max(sublist)) ]);
outfile = ['f' outfile];

%% Low-pass filter (30 Hz = for ERP) (if using manual artefacts, put this after prepro2)
S        = [];
S.D      = D;
S.type   = 'butterworth';
S.band   = 'low';
S.freq   = 30;
S.dir    = 'twopass';
S.order  = 5;
S.prefix = 'f';
D = spm_eeg_filter(S);
disp(['Low-pass filtered (1/2) subject ' num2str(subi) '/' num2str(max(sublist))]);
outfile = ['f' outfile];

%% Filter artefactual jumps in channels amplitudes
if jumps == 1
    S                   = [];
    S.D                 = [outdirdata outfile];
    S.channels          = 'EEG';
    D = spm_eeg_remove_jumps(S);
    disp(['Removed artefacts (jumps) for subject ' num2str(subi) '/' num2str(max(sublist)) ]);
    outfile = [prefix_jumps outfile];
end

%% Downsample (256 Hz)
S             = [];
S.D           = [outdirdata outfile];
S.fsample_new = 256;
S.prefix      = 'd';
D = spm_eeg_downsample(S);
disp(['Downsampled subject ' num2str(subi) '/' num2str(max(sublist)) ]);
outfile = ['d' outfile];

%% Montage
% Mark badchannels A3 &/or A4 -Extra step due to issues with these
% channels on the EEG machine- + Create average reference & average
% EOG montage structure

D = spm_eeg_load([outdirdata outfile]);

% if ~isempty(badchanind{subi})
%     D = badchannels(D,badchanind{subi},1);
% end

S = [];
S.D                  = D;

% re-reference (average)
% ---------------------------------------------------
%
if ~isempty(D.badchannels)
    G = setdiff(indchantype(D,'EEG'),D.badchannels);
    R = G;
    
    T = eye(numel(indchantype(D,'EEG')));
    T(G,R) = T(G,R) - 1/length(R);
    
    S.montage.tra = T; % mxn transformation matrix
    %         S.montage.labelorg = D.chanlabels(D.meegchannels);
    %         S.montage.labelnew = D.chanlabels(D.meegchannels);
    
else
    neeg                 = numel(indchantype(D,'EEG')); % get number of EEG channels (should be 64)
    S.montage.tra        = eye(neeg) - ones(neeg)/neeg; % diagonal matrix minus 1/64-th of all other channels
end
S.montage.tra(65,66) = -1; % set EXG1 as the difference between EXG1 and EXG2 (EOG channel)
S.montage.tra(65,65) = 1;
%     S.montage.labelorg = D.chanlabels(D.meegchannels);
%     S.montage.labelnew = D.chanlabels(D.meegchannels);
S.montage.labelorg   = D.chanlabels([indchantype(D,'EEG') indchantype(D,'EOG')]);
S.montage.labelnew   = D.chanlabels([indchantype(D,'EEG') 65]);
S.montage.name       = 'AverageRef';
% Now write montage
S.mode             = 'write';
S.blocksize        = 655360;
S.prefix           = 'M';
S.keepothers = 1;
S.updatehistory = 1;
spm_eeg_montage(S);

disp(['Defined montage for subject ' num2str(subi) '/' num2str(max(sublist)) ]);

outfile = ['M' outfile];

%     filename = 'dffRAND20TO10_sub20_run1.dat';

%% Trial definition
%Read the varying length of the triggers to define the trials
S               = [];
S.D             = [outdirdata outfile '.mat'];
S.timewin       = timewindow_td; % time-window: 500 ms before trigger, 1000 ms post-offset.
for ii=1:numel(condlist)
    S.trialdef(ii).conditionlabel = condlist(ii);
    S.trialdef(ii).eventtype      = 'trigger_up';
    S.trialdef(ii).eventvalue     = triglist(ii);
    S.trialdef(ii).trlshift       = 0;
    S.filename                    = [outdirraw outfile_folder exptName '_sub' num2str(subi) '_loc.bdf'];
end
[trl, labels] = spm_eeg_definetrial_chaitlab(S);
disp(['Defined trials for subject ' num2str(subi) '/' num2str(max(sublist)) ]);

%% Eyeblink removal (Must before epoch)
if eyeblink == 1
    %Fiducial coordinates
    TemplateFids = [1 85 -40; -87 -12 -50; 87 -12 -50];
    
    spm_jobman('initcfg');
    batch{1}.spm.meeg.source.headmodel.D = {[outdirdata outfile '.mat']};
    batch{1}.spm.meeg.source.headmodel.val = 1;
    batch{1}.spm.meeg.source.headmodel.comment = 'FidsTempl';
    batch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
    batch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
    batch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
    batch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = TemplateFids(1,:);
    % batch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
    batch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    batch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = TemplateFids(2,:);
    % batch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'FIL_CTF_L';
    batch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    batch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = TemplateFids(3,:);
    % batch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'FIL_CTF_R';
    batch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    batch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    batch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    %spm_jobman('run',batch);
    clear batch
    
    %             eyechan = 'MRF14'; % subject 9
    %eyechan = 'MRT31'; % subject 10
    %eyechan = 'MLT31'; % subject 11
    
    % Find blink artefacts in frontal channels
    batch{1}.spm.meeg.preproc.artefact.D = {[outdirdata outfile '.mat']};
    batch{1}.spm.meeg.preproc.artefact.mode = 'mark';
    batch{1}.spm.meeg.preproc.artefact.badchanthresh = 0.2;
    batch{1}.spm.meeg.preproc.artefact.append = false;
    batch{1}.spm.meeg.preproc.artefact.methods.channels{1}.chan = 'EXG1'; %'VEOG'; % eyechan --> 'VEOG'
    batch{1}.spm.meeg.preproc.artefact.methods.fun.eyeblink.threshold = 2;
    batch{1}.spm.meeg.preproc.artefact.methods.fun.eyeblink.excwin = 0;
    batch{1}.spm.meeg.preproc.artefact.prefix = 'a';
    spm_jobman('run',batch);
    clear batch
    
    eyefile = ['a' outfile];
    
    % Create new epoched dataset containing only eye blink data
    batch{1}.spm.meeg.preproc.epoch.D = {[outdirdata eyefile '.mat']};
    batch{1}.spm.meeg.preproc.epoch.trialchoice.define.timewin = [-500 500];
    batch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.conditionlabel = 'Eyeblink';
    batch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.eventtype = 'artefact_eyeblink';
    batch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.eventvalue = 'EXG1'; %'VEOG'; % eyechan --> 'VEOG'
    batch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef.trlshift = 0;
    batch{1}.spm.meeg.preproc.epoch.bc = 0;
    batch{1}.spm.meeg.preproc.epoch.eventpadding = 0;
    batch{1}.spm.meeg.preproc.epoch.prefix = 'eyeblink';
    spm_jobman('run',batch);
    clear batch
    
    eyefile = ['eyeblink' eyefile];
    
    % Compute spatial confounds in eyeblink dataset
    batch{1}.spm.meeg.preproc.sconfounds.D(1) = {[outdirdata eyefile '.mat']};
    batch{1}.spm.meeg.preproc.sconfounds.mode{1}.svd.timewin = [-Inf Inf];
    batch{1}.spm.meeg.preproc.sconfounds.mode{1}.svd.threshold = NaN;
    batch{1}.spm.meeg.preproc.sconfounds.mode{1}.svd.ncomp = 4;
    spm_jobman('run',batch);
    clear batch
    
    % Find spatial confounds in dataset of interest
    batch{1}.spm.meeg.preproc.sconfounds.D = {[outdirdata outfile '.mat']};
    batch{1}.spm.meeg.preproc.sconfounds.mode{1}.spmeeg.conffile(1) = {[outdirdata eyefile '.mat']};
    spm_jobman('run',batch);
    clear batch
    
    % Remove eye blink from dataset of interest
    batch{1}.spm.meeg.preproc.correct.D = {[outdirdata outfile '.mat']};
    batch{1}.spm.meeg.preproc.correct.mode = 'berg';
    batch{1}.spm.meeg.preproc.correct.prefix = 'T';
    spm_jobman('run',batch);
    clear batch
    
    outfile = [prefix_eyelink outfile];
    display(['removed blinks for subject ' num2str(subi) '/' num2str(length(sublist))]);
end

%% Epoch + baseline correct
S                   = [];
S.D                 = [outdirdata outfile];
S.trl               = trl;
S.conditionlabels   = labels;
S.bc                = 1;
S.timewin           = timewindow_epoch; % 1000 ms before onset to 1000 ms after offset.
spm_eeg_epochs(S);
disp(['Epoched subject ' num2str(subi) '/' num2str(max(sublist))]);
outfile = ['e' outfile];

%% Baseline correction
if baselinec == 1
    S                   = [];
    S.D                 = [outdirdata outfile];
    S.timewin           = timewindow_bc; % 1000 ms before onset to 1000 ms after offset.
    spm_eeg_bc(S);
    disp(['Baseline correction ' num2str(subi) '/' num2str(max(sublist)) ]);
    outfile = ['b' outfile];
end

%% Artefact removal: zscore-based automatic outliers removal method
S       = [];
S.D     = [outdirdata outfile];
visual  = 0; % value of 1 for manual outlier removal (better but more time-consuming)

switch visual
    case 0
        thresh                          = 7;
        S.D                             = [outdirdata outfile];
        S.methods(1).fun                = 'zscore';
        S.methods(1).channels           = 'EEG';
        S.methods(1).settings.threshold = thresh;
        S.methods(1).settings.excwin    = 100;
        S.methods(2).fun                = 'zscorediff';
        S.methods(2).channels           = 'EEG';
        S.methods(2).settings.threshold = thresh;
        S.methods(2).settings.excwin    = 100;
        S.mode                          = 'reject';
        S.append                        = 0; % overwrite previous markings
        S.badchanthresh                 = .2;
        S.prefix                        = 'a';
        D = spm_eeg_artefact(S);
        
        % Case too many bad trials
        while numel(D.badtrials) > 8
            disp([char(10) '*******'])
            disp(['Removing outliers from subject ' num2str(subi) ...
                '/' num2str(max(sublist)) ]);
            disp('Rejected too many trials, retrying with more liberal threshold...')
            disp(['*******' char(10)])
            thresh                          = thresh+.5;
            S.D                             = [outdirdata outfile];
            S.methods(1).settings.threshold = thresh;
            S.methods(2).settings.threshold = thresh;
            S.methods(2).settings.excwin    = 100;
            S.badchanthresh                 = .2;
            D = spm_eeg_artefact(S);
        end
        
    case 1
        S.method  = 'summary';
        S.prefix  = 'a';
        S.latency = [-Inf Inf];
        S.channel = 'EEG';
        D = fVisualOutliers(S);
    otherwise
        disp('Problem at artefact removal step. Aborting ! ')
        return;
end

%% Move file
temp            = dir([root 'c' 'a' prefix_baselinec 'e' prefix_eyelink 'Md' prefix_jumps 'ff' 'RAND20TO10_sub' num2str(subi) '_loc']); % caeMfdf --> aeMfdfSS_sub
D               = spm_eeg_load([outdirdata root 'a' prefix_baselinec 'e' prefix_eyelink 'Md' prefix_jumps 'ff' 'RAND20TO10_sub' num2str(subi) '_loc.mat']);
outfile         = ['subject' num2str(subi) '_loc'];
D = move(D,[root outfile]); % For time frequency
%% Rename channels to 10-20
load([root 'Inputs/Electrodes_labels_Biosemi64.mat']);
D = spm_eeg_load([outdirdata outfile]);
D = chanlabels(D,1:64,Electrodes_labels_Biosemi64);
save(D);

%% Average
S=[]; S.D = [root outfile];
S.robust = 0;
D = spm_eeg_average(S);
disp(['Average ' num2str(subi) '/' num2str(max(sublist)) ]);
%filename of the output will be msubject1
end