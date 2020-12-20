% Test mTRF updates using original speech data from Natural Speech
% This computes a forward model using the envelope as input
% Nate Zuk (2020)

addpath('../mtrf/');

sbj = 'Subject19'; % subject 10 has good data, d-prime performance > 4
    % Subject8 also ok, d-prime ~2
    % Subject5 good, d-prime ~3, but interesting TRF shape
eegpth = ['/Volumes/Untitled/Natural Speech/EEG/' sbj '/'];
stimpth = '/Volumes/Untitled/Natural Speech/Stimuli/Envelopes/';
Fs = 128; % sampling rate of envelope and EEG
maxtr = 20; % maximum trial to include in the analysis
chan = 85; % channel to use for prediction
freq_range = [1 15]; % frequency range of the bandpass filter for the EEG
nperm = 100;

% Modeling parameters
lambdas = 10.^(-2:8); % regularization parameters to test
% lambdas = 10^18;
tmin = -100;
tmax = 350;
map = 1; % compute forward direction
t = floor(tmin/1000*Fs):ceil(tmax/1000*Fs);

stims = cell(maxtr,1); % one for each trial
EEGs = cell(maxtr,1);
for n = 1:maxtr
    fprintf('Run %d\n',n);
    % load stimulus
    sfl = sprintf('audio%d_128Hz',n);
    Sd = load([stimpth sfl]);
    stims{n} = Sd.env;
    % load eeg data
    efl = sprintf('%s_Run%d',sbj,n);
    Ed = load([eegpth efl]);
    % detrend and refrence data
    ref = mean(Ed.mastoids,2);
    eeg = Ed.eegData-ref*ones(1,128); % remove reference from mastoids
    % Filter the EEG
    % detrend
    eeg = detrend(eeg(:,chan)); % remove linear trend
    % bandpass filter
    eeg = eeg_bandpass(eeg,Fs,'highpass_cf',freq_range(1),'lowpass_cf',freq_range(2));
    % Truncate stim and eeg to the same length
    min_len = min([length(stims{n}) size(eeg,1)]);
    stims{n} = stims{n}(1:min_len);
    eeg = eeg(1:min_len,:); % set to the same length as the stimulus
    % Transform into the principal components of the scalp topography
    EEGs{n} = eeg;
    % zscore the principal components
    EEGs{n} = zscore(EEGs{n});
    % zscore the stimulus -- to make the orig and new mTRFcrossval
    % comparable
    stims{n} = zscore(stims{n});
end

% Remove variables to save RAM
clear eeg Sd Ed ref

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

% Save the results
svpth = 'ns_env_res/';
svfl = sprintf('%s_frwdtrf',sbj);
save([svpth svfl],'all_mdl','opt_lmb','test_r','circ_test_r','mean_cv_r',...
    'lambdas','chan','freq_range','maxtr','Fs','tmin','tmax');