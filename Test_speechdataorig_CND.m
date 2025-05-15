% Test mTRF updates using original speech data from Natural Speech
% This computes a forward model using the envelope as input
% Nate Zuk (2020)

addpath(genpath('../mTRF-Toolbox/'));

sbj = 10; % subject 10 has good data, d-prime performance > 4
    % Subject8 also ok, d-prime ~2
    % Subject5 good, d-prime ~3, but interesting TRF shape
eegpth = 'C:\Users\zsxbo\Projects\cnsp-code\datasets\LalorNatSpeech\dataCND\';
stimpth = '..\cnsp-code\datasets\LalorNatSpeech\dataCND\';
Fs = 128; % sampling rate of envelope and EEG
maxtr = 20; % maximum trial to include in the analysis
chan = 85; % channel to use for prediction (Fz)
freq_range = [1 15]; % frequency range of the bandpass filter for the EEG
nperm = 100;

% Modeling parameters
lambdas = [0 10.^(-2:8)]; % regularization parameters to test
% lambdas = 10^18;
tmin = -100;
tmax = 350;
map = 1; % compute forward direction
t = floor(tmin/1000*Fs):ceil(tmax/1000*Fs);

% Load the stimulus envlope
d = load([stimpth 'dataStim']);
stims = d.stim.data(1,:); % contains the stimulus envelope on each trial
clear d

% Load the EEG for this participant
eeg_d = load([eegpth sprintf('dataSub%d',sbj)]);
EEGs = eeg_d.eeg.data;
refs = eeg_d.eeg.extChan{1}.data; % reference electrodes (Mastoids)
clear d

% Do some very basic preprocessing on the data
for n = 1:length(EEGs)
    % reference to the average of the mastoids
    ref_avg = mean(refs{n},2);
    EEGs{n} = EEGs{n}-ref_avg;
    % detrend linear slope
    EEGs{n} = detrend(EEGs{n}(:,chan));
    % bandpass filter
    EEGs{n} = eeg_bandpass(EEGs{n},Fs,'highpass_cf',freq_range(1),'lowpass_cf',freq_range(2));
    fprintf('Processed trial %d\n',n);
end

clear refs

% Truncate so the stimulus and EEG is the same length
for n = 1:maxtr
    min_len = min([length(stims{n}) length(EEGs{n})]);
    stims{n} = stims{n}(1:min_len);
    EEGs{n} = EEGs{n}(1:min_len);
    % normalize
%     stims{n} = zscore([0; diff(stims{n})]); % by root mean square, so it stays positive
    stims{n} = zscore(stims{n});
    EEGs{n} = zscore(EEGs{n}); % z-score, so it is centered and std 1
end

% Use leave-one-trial-out to train and test the prediction model on 1
% channel
mdl = cell(maxtr,1);
pred = cell(maxtr,1);
% true_cv_stats = cell(maxtr,1);
mean_cv_r = NaN(length(lambdas),maxtr);
opt_lmb = NaN(maxtr,1);
test_r = NaN(maxtr,1);
all_mdl = NaN(length(t),maxtr);
mdltm = tic;
fprintf('Testing and training on %d trials',maxtr);
for n = 1:maxtr
    fprintf('.');
    test_tr = n;
    train_trs = setxor(1:maxtr,n);
    true_cv_stats = mTRFcrossval(stims(train_trs),EEGs(train_trs),Fs,1,tmin,tmax,lambdas,...
        'verbose',0);
    mean_cv_r(:,n) = mean(true_cv_stats.r);
    % identify the optimal lambda based on the average correlation coefficient
%             opt_idx = find(mean(true_cv_stats.err)==min(mean(true_cv_stats.err)),1,'first');
        % based on minimum MSE instead of correlation
    opt_idx = find(mean_cv_r(:,n)==max(mean_cv_r(:,n)),1,'first');
%             opt_idx = 1;
%             fprintf('-- Optimal lambda = %d\n',lambdas(opt_idx));
    mdl{n} = mTRFtrain(stims(train_trs),EEGs(train_trs),Fs,1,tmin,tmax,lambdas(opt_idx),...
        'verbose',0);
    % save the prediction on each iteration (this is equivalent to
    % re-running the prediction on the same model
    [pred{n},test_stats] = mTRFpredict(stims(test_tr),EEGs(test_tr),mdl{n},'verbose',0);
    opt_lmb(n) = lambdas(opt_idx);
    test_r(n) = test_stats.r;
    all_mdl(:,n) = mdl{n}.w;
end
fprintf('Completed @ %.3f s\n',toc(mdltm));

disp('Creating the null distribution...');
circ_test_r = NaN(nperm,1);
for n = 1:nperm
    circ_pair_tr = randi(maxtr);
    rnd_shft = randi(length(EEGs{circ_pair_tr}));
    shft_resp = circshift(EEGs{circ_pair_tr},rnd_shft);
    circ_test_r(n) = corr(pred{circ_pair_tr},shft_resp);
end

% Display d-prime prediction accuracy
mndiff = mean(test_r)-mean(circ_test_r);
stdiff = sqrt(0.5*(var(test_r)+var(circ_test_r)));
dp = mndiff/stdiff;

figure
plot(t/Fs*1000,all_mdl)
set(gca,'FontSize',14);
xlabel('Delay (ms)')
ylabel('TRF weight');
title(sprintf('%s, d-prime = %.3f',sbj,dp));

if length(lambdas)>1 % if there is more than one lambda value, plot the CV tuning curves
    figure
    hold on
    plot(lambdas,mean_cv_r,'k');
    set(gca,'FontSize',14,'XScale','log');
    xlabel('\lambda');
    ylabel('Average CV r across folds');
    title(sprintf('%s, optimal lambda = %d',sbj,mode(opt_lmb)));
end

% Calculate SNR
var_s = NaN(maxtr,1);
var_n = NaN(maxtr,1);
for n = 1:maxtr
    var_s = var(pred{n});
    var_n = var(EEGs{n});
end
SNR = 10*log10(mean(var_s)/mean(var_n));
fprintf('SNR = %.1f\n',SNR);

% Plot an example prediction and EEG
figure
set(gcf,'Position',[100 100 900 350]);
hold on
tstim = (0:length(pred{1})-1)/Fs;
plot(tstim,EEGs{1},'k');
hold on
plot(tstim,pred{1},'r','LineWidth',1.5);
set(gca,'XLim',[5,20]);
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original EEG','Prediction');
exportgraphics(gca,sprintf('fig/ExamplePred1_sbj%d.pdf',sbj),'ContentType','vector');

% Save the results
svpth = 'ns_env_res/';
svfl = sprintf('Sbj%d_frwdtrf',sbj);
save([svpth svfl],'all_mdl','opt_lmb','test_r','circ_test_r','mean_cv_r',...
    'lambdas','chan','freq_range','maxtr','Fs','tmin','tmax','SNR');