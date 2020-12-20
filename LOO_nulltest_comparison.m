% Simulate the effect of different types of permutations by simulated a bunch
% of trials with an embedded signal created by a convolution with added
% noise.  Use leave-one-out testing, which is closer to what we normally do.
% After testing, get the null distribution by either:
% 1) shuffling the testing data and using the same model to test
% 2) circularly shifting the testing data, but keeping stimulus-response
% pairs ths same and using the original model
% 3) shuffling both training and testing data to recompute a model (with the
% same amount of regularization) and test
% 4) circularly shifting both the training and testing data to recompute a
% model and test
% Nate Zuk (2020)

addpath('../mtrf/');

dur = 180; % duration of the stimulus and response
Fs = 128; % sampling frequency of the signals
ntr = 20; % number of trials
freq_range = [0.1 10]; % frequency range of the signals
snr = 0.1; % signal to noise ratio in the response
lambdas = [0 10.^(0:14)];
nperm = 100; % number of times to shuffle the data and get null testing values
quantiles_to_plot = [0.05 0.95];
% drift_mag = 10;

%%% Generate a TRF model
resp_dur = 250; % duration of the response (in ms)
resp_t = (0:ceil(resp_dur/1000*Fs))/Fs;
resp_frq = 6; % frequency of the response
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
%     stim{n} = zscore(bandpass_noise(freq_range,dur,Fs));
    s = bandpass_noise(freq_range,2*dur,Fs);
    stim{n} = zscore(s(1:dur*Fs));
    % add a linear drift
%     stim{n} = stim{n}+(0:dur*Fs-1)'/(dur*Fs)*drift_mag;
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
disp('Training on true data...')
true_tm = tic;
pred = cell(ntr,1);
mdl = cell(ntr,1);
opt_idx = NaN(ntr,1);
test_r = NaN(ntr,1);
for n = 1:ntr
    test_tr = n;
    train_trs = setxor(1:ntr,n);
    true_cv_stats = mTRFcrossval(stim(train_trs),resp(train_trs),Fs,1,0,resp_dur,lambdas,'verbose',0);
    % identify the optimal lambda based on the average correlation coefficient
    opt_idx(n) = find(mean(true_cv_stats.r)==max(mean(true_cv_stats.r)),1,'first');
%     fprintf('-- Optimal lambda = %d\n',lambdas(opt_idx));
    mdl{n} = mTRFtrain(stim(train_trs),resp(train_trs),Fs,1,0,resp_dur,lambdas(opt_idx(n)),'verbose',0);
    % save the prediction on each iteration (this is equivalent to
    % re-running the prediction on the same model
    [pred{n},test_stats] = mTRFpredict(stim(test_tr),resp(test_tr),mdl{n},'verbose',0);
    test_r(n) = test_stats.r;
end
fprintf('* Completed training and testing on true data in @ %.3f s\n',toc(true_tm));

%% Plot both models to see how similar they are
figure
hold on
plot(resp_t*1000,true_trf/rms(true_trf),'k');
for n = 1:ntr
    plot(mdl{n}.t,mdl{n}.w/rms(mdl{n}.w),'b');
end
set(gca,'FontSize',14);
xlabel('Delay (ms)');
ylabel('RMS normalized model weights');

%% Permutations
%%% Permute 1: Shuffle the testing data to get a null distribution of test
%%% values
disp('** Permute 1 **');
perm_i = tic;
shuff_test_r = NaN(nperm,1);
for n = 1:nperm
    % randomly pick two trials to pair for testing
    shuff_test_tr = randi(ntr);
    shuff_mdl_tr = randi(ntr); 
    shuff_test_r(n) = corr(pred{shuff_mdl_tr},resp{test_tr});
end
fprintf('* Completed @ %.3f s\n',toc(perm_i));

%%% Permute 2: Randomly select a pairing of trials (this time the correct
%%% pairing) but circularly shift the response
disp('** Permute 2 **');
perm_ii = tic;
circ_test_r = NaN(nperm,1);
for n = 1:nperm
    circ_pair_tr = randi(ntr);
    rnd_shft = randi(length(resp{circ_pair_tr}));
    shft_resp = circshift(resp{circ_pair_tr},rnd_shft);
    circ_test_r(n) = corr(pred{circ_pair_tr},shft_resp);
end
fprintf('* Completed @ %.3f s\n',toc(perm_ii));

%%% Permute 3: Shuffle both the training and testing data, recompute the
%%% model, and test on the left out (shuffled) data
disp('** Permute 3 **');
perm_iii = tic;
use_lambda = lambdas(mode(opt_idx));
shuff_all_r = NaN(nperm,1);
for n = 1:nperm
    shuff_trials = randperm(ntr);
    % randomly select a pair for testing
    test_pair = randi(ntr);
    train_pairs = setxor(1:ntr,test_pair);
    shuff_mdl = mTRFtrain(stim(train_pairs),resp(shuff_trials(train_pairs)),...
        Fs,1,0,resp_dur,use_lambda,'verbose',0);
    [~,shuff_stats] = mTRFpredict(stim(test_pair),resp(shuff_trials(test_pair)),...
        shuff_mdl,'verbose',0);
    shuff_all_r(n) = shuff_stats.r;
end
fprintf('* Completed @ %.3f s\n',toc(perm_iii));

%%% Permute 4: Randomly circularly shift the training and testing data (use
%%% the same shift for all trials), recompute the model, and test on the
%%% left out pair of data
disp('** Permute 4 **');
perm_iv = tic;
circ_all_r = NaN(nperm,1);
for n = 1:nperm
    rnd_shft = randi(dur*Fs); % use the length of all trials
    % circularly shift all of the responses
    shft_resp = cell(ntr,1);
    for m = 1:ntr
        shft_resp{m} = circshift(resp{m},rnd_shft);
    end
    % randomly select a pair for testing
    test_pair = randi(ntr);
    train_pairs = setxor(1:ntr,test_pair);
    circ_mdl = mTRFtrain(stim(train_pairs),shft_resp(train_pairs),...
        Fs,1,0,resp_dur,use_lambda,'verbose',0);
    [~,circ_stats] = mTRFpredict(stim(test_pair),shft_resp(test_pair),...
        circ_mdl,'verbose',0);
    circ_all_r(n) = circ_stats.r;
end
fprintf('* Completed @ %.3f s\n',toc(perm_iv));

% Show all of the results
fprintf('\n-- Median [%dth %dth]\n',quantiles_to_plot(1)*100,quantiles_to_plot(2)*100);
fprintf('True r: %.3f [%.3f %.3f]\n',median(test_r),quantile(test_r,...
    quantiles_to_plot(1)),quantile(test_r,quantiles_to_plot(2)));
fprintf('Permute 1 (shuffle testing): %.3f [%.3f %.3f]\n',...
    median(shuff_test_r),quantile(shuff_test_r,quantiles_to_plot(1)),...
    quantile(shuff_test_r,quantiles_to_plot(2)));
fprintf('Permute 2 (circshift testing): %.3f [%.3f %.3f]\n',...
    median(circ_test_r),quantile(circ_test_r,quantiles_to_plot(1)),...
    quantile(circ_test_r,quantiles_to_plot(2)));
fprintf('Permute 3 (shuffle all): %.3f [%.3f %.3f]\n',...
    median(shuff_all_r),quantile(shuff_all_r,quantiles_to_plot(1)),...
    quantile(shuff_all_r,quantiles_to_plot(2)));
fprintf('Permute 4 (circshift all): %.3f [%.3f %.3f]\n',...
    median(circ_all_r),quantile(circ_all_r,quantiles_to_plot(1)),...
    quantile(circ_all_r,quantiles_to_plot(2)));

% Plot the results
figure
hold on
plot([1 5],[0 0],'k--');
% True
md = median(test_r);
uq = quantile(test_r,quantiles_to_plot(2));
lq = quantile(test_r,quantiles_to_plot(1));
errorbar(1,md,md-lq,uq-md,'k.','MarkerSize',20,'LineWidth',2);
% Permute 1
md = median(shuff_test_r);
uq = quantile(shuff_test_r,quantiles_to_plot(2));
lq = quantile(shuff_test_r,quantiles_to_plot(1));
errorbar(2,md,md-lq,uq-md,'k.','MarkerSize',20,'LineWidth',2);
% Permute 2
md = median(circ_test_r);
uq = quantile(circ_test_r,quantiles_to_plot(2));
lq = quantile(circ_test_r,quantiles_to_plot(1));
errorbar(3,md,md-lq,uq-md,'k.','MarkerSize',20,'LineWidth',2);
% Permute 3
md = median(shuff_all_r);
uq = quantile(shuff_all_r,quantiles_to_plot(2));
lq = quantile(shuff_all_r,quantiles_to_plot(1));
errorbar(4,md,md-lq,uq-md,'k.','MarkerSize',20,'LineWidth',2);
% Permute 4
md = median(circ_all_r);
uq = quantile(circ_all_r,quantiles_to_plot(2));
lq = quantile(circ_all_r,quantiles_to_plot(1));
errorbar(5,md,md-lq,uq-md,'k.','MarkerSize',20,'LineWidth',2);

set(gca,'FontSize',14,'XTick',1:5,'XTickLabelRotation',45,...
    'XTickLabel',{'true','shuffle test','circshift test','shuffle all','circshift all'});
ylabel('Pearsons r');