% Examine the effect of using the optimal lambda when computing null
% prediction accuracies vs recomputing the lambda with cross-validation for each
% permutation of trials.
% Nate Zuk (2020)

addpath('../mtrf/');

dur = 30; % duration of the stimulus and response
Fs = 100; % sampling frequency of the signals
ntr = 30; % number of trials
freq_range = [1 10]; % frequency range of the signals
snr = 1; % signal to noise ratio in the response
lambdas = [0 10.^(-2:14)];
nperm = 50; % number of times to shuffle the data and get null testing values

%%% Generate a TRF model
resp_dur = 250; % duration of the response (in ms)
resp_t = (0:ceil(resp_dur/1000*Fs))/Fs;
resp_frq = 8; % frequency of the response
% the trf is an exponentially decaying sinusoid
true_trf = sin(2*pi*resp_frq*resp_t).*exp(-resp_t/(resp_dur/1000/5));

%%% Simulate low-frequency data, which requires regularization
%%% There should be an embedded signal that is produced by convolution with
%%% some TRF-like filter, and added noise with the same frequency
%%% distribution.
% stimulus is random bandpass noise
stim = cell(ntr,1);
resp = cell(ntr,1);
lags = round(resp_t*Fs);
for n = 1:ntr
    % make sure the variance of the stimulus and response are both 1
    stim{n} = zscore(bandpass_noise(freq_range,dur,Fs));
%     trfout = zscore(conv(stim{n},true_trf,'same')); %%% Didn't work the
%     same way? Using 'same' shifts the convolution backwards, so the
%     fitted model won't match up
    X = lagGen(stim{n},lags);
    trfout = zscore(X*true_trf');
    ns = zscore(bandpass_noise(freq_range,dur,Fs));
    resp{n} = zscore(trfout*snr+ns);
end

%%% Train, with cross-validation to get the optimum model, and test on
%%% left-out data
disp('Using leave-one-out to train and test model...')
test_r = NaN(ntr,1);
opt_lmb = NaN(ntr,1);
for n = 1:ntr
    disp('Testing trial %d...');
    true_tm = tic;
    train_trs = setxor(1:ntr,n);
    test_trs = n;
    true_cv_stats = mTRFcrossval(stim(train_trs),resp(train_trs),Fs,1,0,resp_dur,lambdas,...
        'verbose',0);
    % identify the optimal lambda based on the average correlation coefficient
    opt_idx = find(mean(true_cv_stats.r)==max(mean(true_cv_stats.r)),1,'last');
    fprintf('-- Optimal lambda = %d\n',lambdas(opt_idx));
    opt_lmb(n) = lambdas(opt_idx);
    mdl = mTRFtrain(stim(train_trs),resp(train_trs),Fs,1,0,resp_dur,opt_lmb(n),'verbose',0);
    [~,test_stats] = mTRFpredict(stim(test_trs),resp(test_trs),mdl,'verbose',0);
    test_r(n) = test_stats.r;
    fprintf('* Completed training and testing on true data in @ %.3f s\n',toc(true_tm));
end

% plot both models to see how similar they are
figure
hold on
plot(resp_t*1000,true_trf/rms(true_trf),'k');
plot(mdl.t,mdl.w/rms(mdl.w),'b','LineWidth',2);
set(gca,'FontSize',14);
xlabel('Delay (ms)');
ylabel('RMS normalized model weights');

%%% Permute 1: Permute trials, recompute model on all trials with one left
%%% out using lambda randomly selected from optimal lambdas, and test on
%%% left out trial
disp('Permute 1: Use optimal lambdas of true models...');
nulltm = tic;
shuff_test_r = NaN(nperm,1);
for n = 1:nperm
    null_stim_test = randi(ntr);
    null_stim_train = setxor(1:ntr,null_stim_test);
    % shuffle stimulus trials
    null_stim_train = null_stim_train(randperm(ntr-1));
    null_resp_test = randi(ntr);
    null_resp_train = setxor(1:ntr,null_resp_test);
    lmb = opt_lmb(randi(ntr));
    null_mdl = mTRFtrain(stim(null_stim_train),resp(null_resp_train),Fs,1,0,resp_dur,lmb,'verbose',0);
    [~,null_stats] = mTRFpredict(stim(null_stim_test),resp(null_resp_test),null_mdl,'verbose',0);
    shuff_test_r(n) = null_stats.r;
end
fprintf('Completed @ %.3f s\n',toc(nulltm));

%%% Permute 2: Permute trials, recompute optimal lambda on all trials with
%%% one left out, and test on left out trial
disp('Permute 2: Recompute lambda for each permutation...');
nulltm = tic;
shuff_cv_r = NaN(nperm,1);
shuff_cv_lmb = NaN(nperm,1);
for n = 1:nperm
    null_stim_test = randi(ntr);
    null_stim_train = setxor(1:ntr,null_stim_test);
    % shuffle stimulus trials
    null_stim_train = null_stim_train(randperm(ntr-1));
    null_resp_test = randi(ntr);
    null_resp_train = setxor(1:ntr,null_resp_test);
    null_cv_stats = mTRFcrossval(stim(null_stim_train),resp(null_resp_train),Fs,1,0,resp_dur,lambdas,...
        'verbose',0);
    % identify the optimal lambda based on the average correlation coefficient
    null_idx = find(mean(null_cv_stats.r)==max(mean(null_cv_stats.r)),1,'last');
    shuff_cv_lmb(n) = lambdas(null_idx);
    null_mdl = mTRFtrain(stim(null_stim_train),resp(null_resp_train),Fs,1,0,resp_dur,shuff_cv_lmb(n),...
        'verbose',0);
    [~,null_stats] = mTRFpredict(stim(null_stim_test),resp(null_resp_test),null_mdl,'verbose',0);
    shuff_cv_r(n) = null_stats.r;
end
fprintf('Completed @ %.3f s\n',toc(nulltm));

% Plot reconstruction accuracies
figure
hold on
plot(zeros(ntr,1),test_r,'b.','MarkerSize',16);
plot(ones(nperm,1),shuff_test_r,'k.','MarkerSize',16);
plot(2*ones(nperm,1),shuff_cv_r,'k.','MarkerSize',16);
set(gca,'FontSize',14,'XTick',[0 1 2],'XTickLabel',{'True','Pick \lambda','Recomp \lambda'});
ylabel('Reconstruction accuracy');