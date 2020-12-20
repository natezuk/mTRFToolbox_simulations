% Simulate the effect of different types of permutations by simulated a bunch
% of trials with an embedded signal created by a convolution with added
% noise.  Use leave-one-out testing, which is closer to what we normally do.
% Try this for various SNRs and numbers of trials (amount of data)
% and plot how the correlation between the predicted and true model varies
% with d-prime prediction accuracy.
% Nate Zuk (2020)

addpath('../mtrf/');
%addpath('~/Documents/Matlab/shadedErrorBar/');

dur = 60; % duration of the stimulus and response in each trial
Fs = 128; % sampling frequency of the signals
freq_range = [0.1 4]; % frequency range of the signals
%%% highpass cutoff of signal has higher frequency
resp_dur = 650; % duration of the response (in ms) (response starts at 0 delay)
snr = logspace(-1,-2.5,7); % signal to noise ratio in the response
ntr = 2.^(2:6); % number of trials (must have at least two trials)
trf_range = [-100 750]; % range of delays to include in the TRF estimate
lambdas = 10.^(-2:8);
niter = 100; % number of times to repeat the analysis with new simulated data
nperm = 100; % number of times to shuffle the data and get null testing values

% Load example EEG (taken from NeuroMusic)
s = load('example_eeg');
eeg_sample = s.sample_eeg(:,1); % use channel 1
eFs = s.eFs; % sampling rate of example EEG
clear s
% detrend
eeg_sample = detrend(eeg_sample); % remove linear trend
% bandpass filter
eeg_sample = eeg_bandpass(eeg_sample,eFs,'highpass_cf',freq_range(1),'lowpass_cf',freq_range(2));
%%% (23-10-2020) EEG noise and stimulus have same spectrum
eeg_sample = resample(eeg_sample,Fs,eFs); % resample the EEG to the correct Fs
eeg_sample = eeg_sample(1:dur*Fs); % truncate to signal duration

%%% Generate a TRF model
resp_t = (0:ceil(resp_dur/1000*Fs))/Fs;
resp_frq = 2; % frequency of the response
% the trf is an exponentially decaying sinusoid
true_trf = sin(2*pi*resp_frq*resp_t).*exp(-resp_t/(resp_dur/1000/5));

% compute the range of trf delays, in indexes
lags = floor(trf_range(1)/1000*Fs):ceil(trf_range(2)/1000*Fs);

total_mdl = cell(length(ntr),1);
test_r = cell(length(ntr),1);
opt_lmb = cell(length(ntr),1);
circ_test_r = cell(length(ntr),1);
for kk = 1:length(ntr)
    % get the number of trials
    N = ntr(kk);
    fprintf('-- %d trials --\n',N);
    % preallocate arrays for this iteration of number of trials
    total_mdl{kk} = NaN(length(lags),niter,length(snr));
    test_r{kk} = NaN(N,niter,length(snr));
    opt_lmb{kk} = NaN(N,niter,length(snr));
    circ_test_r{kk} = NaN(nperm,niter,length(snr));
    for jj = 1:length(snr)
        fprintf('** SNR = %.3g **\n',snr(jj));
        snr_tm = tic;
        fprintf('Iterating %d simulations',niter);
        for ii = 1:niter
            fprintf('.');

            %%% Simulate low-frequency data, which requires regularization
            %%% There should be an embedded signal that is produced by convolution with
            %%% some TRF-like filter, and added noise with the same frequency
            %%% distribution.
            % stimulus is random bandpass noise
            stim = cell(N,1);
            resp = cell(N,1);
            for n = 1:N
                % make sure the variance of the stimulus and response are both 1
                stim{n} = zscore(bandpass_noise(freq_range,dur,Fs));
                X = lagGen(stim{n},resp_t*Fs);
                trfout = zscore(X*true_trf');
                %%% Generate noise by randomizing the phases of the
                %%% EEG sample, so the noise has a EEG-like spectrum
                ns = zscore(randomize_phases(eeg_sample)); 
                resp{n} = zscore(trfout*snr(jj)+ns);
            end

            %%% Train, with cross-validation to get the optimum model, and test on
            %%% left-out data
            true_tm = tic;
            pred = cell(N,1);
            mdl = cell(N,1);
            cv_iter_r = NaN(length(lambdas),N);
            for n = 1:N
                test_tr = n;
                train_trs = setxor(1:N,n);
                true_cv_stats = mTRFcrossval(stim(train_trs),resp(train_trs),Fs,1,...
                    trf_range(1),trf_range(2),lambdas,'verbose',0);
                % identify the optimal lambda based on the average correlation coefficient
                opt_idx = find(mean(true_cv_stats.r)==max(mean(true_cv_stats.r)),1,'last');
                    %%% changed to get the last index, which increases the
                    %%% amount of regularization
                % Save the cross-validation tuning curve, averaged across folds
                cv_iter_r(:,n) = mean(true_cv_stats.r);
                % Create the optimal model
                mdl{n} = mTRFtrain(stim(train_trs),resp(train_trs),Fs,1,trf_range(1),trf_range(2),...
                    lambdas(opt_idx),'verbose',0);
                % save the prediction on each iteration (this is equivalent to
                % re-running the prediction on the same model
                [pred{n},test_stats] = mTRFpredict(stim(test_tr),resp(test_tr),mdl{n},'verbose',0);
                opt_lmb{kk}(n,ii,jj) = lambdas(opt_idx);
                test_r{kk}(n,ii,jj) = test_stats.r;
            end

            %%% Refit a model using the lambda values that was optimal most
            %%% often across trials
            total_opt_lmb = mode(opt_lmb{kk}(:,ii,jj)); % find the 
            opt_mdl = mTRFtrain(stim,resp,Fs,1,trf_range(1),trf_range(2),total_opt_lmb,'verbose',0);
            total_mdl{kk}(:,ii,jj) = opt_mdl.w;

            %%% Circularly shift data to get a null distribution
            for n = 1:nperm
                circ_pair_tr = randi(N);
                rnd_shft = randi(length(resp{circ_pair_tr}));
                shft_resp = circshift(resp{circ_pair_tr},rnd_shft);
                circ_test_r{kk}(n,ii,jj) = corr(pred{circ_pair_tr},shft_resp);
            end
        end
        fprintf('Completed @ %.3f s\n',toc(snr_tm));
    end
end

save('TRFestimate_w_SNRnTr_delta','total_mdl','test_r','opt_lmb','circ_test_r',...
    'snr','lambdas','nperm','dur','Fs','freq_range','resp_dur','resp_frq',...
    'trf_range','ntr');
