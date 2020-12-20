function [mdl_corr,dpr,mn_r,all_mdls,lags] = compute_pred_acc_measures(S)
% For a data structure output from Test_TRFestimate_SNR, compute the
% average r, d-prime accuracies, and the true-to-predicted TRF correlation

% Get analysis variables
snr = S.snr;
ntr = S.ntr;
niter = size(S.test_r{1},2); % number of times the simulation was run
test_r = S.test_r;
circ_test_r = S.circ_test_r;
total_mdl = S.total_mdl;

% get the true response
resp_t = (0:ceil(S.resp_dur/1000*S.Fs))/S.Fs;
% the trf is an exponentially decaying sinusoid
true_trf = sin(2*pi*S.resp_frq*resp_t).*exp(-resp_t/(S.resp_dur/1000/5));

% get the model lags
lags = floor(S.trf_range(1)/1000*S.Fs):ceil(S.trf_range(2)/1000*S.Fs);

dpr = NaN(niter,length(snr),length(ntr));
mdl_corr = NaN(niter,length(snr),length(ntr));
mn_r = NaN(niter,length(snr),length(ntr));
all_mdls = NaN(niter,length(snr),length(ntr));
for kk = 1:length(ntr)
    for jj = 1:length(snr)
        for ii = 1:niter
            % compute dprim
            mdiff = mean(test_r{kk}(:,ii,jj))-mean(circ_test_r{kk}(:,ii,jj));
            stdcmb = sqrt(0.5*(var(test_r{kk}(:,ii,jj))+var(circ_test_r{kk}(:,ii,jj))));
            dpr(ii,jj,kk) = mdiff/stdcmb;
            % get the indexes that correspond to the start and end of the
            % true trf
            trf_start_idx = find(resp_t(1)*S.Fs==lags);
            trf_end_idx = find(resp_t(end)*S.Fs==lags);
            mdl_corr(ii,jj,kk) = corr(total_mdl{kk}(trf_start_idx:trf_end_idx,ii,jj),true_trf');
            % save the original model as well
            all_mdls{ii,jj,kk} = total_mdl{kk}(:,ii,jj);
            % compute the average r value across trials
            mn_r(ii,jj,kk) = mean(test_r{kk}(:,ii,jj));
        end
    end
end