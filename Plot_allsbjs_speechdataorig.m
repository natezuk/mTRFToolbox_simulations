% Load the forward modeling envelope analysis from Natural Speech, and plot
% all of the subjects along with their d-prime prediction accuracies
% Nate Zuk (2020)

sbjs = 1:19;

respth = 'ns_env_res/';

tmin = -100;
tmax = 350;
Fs = 128;
lags = floor(tmin/1000*Fs):ceil(tmax/1000*Fs);
maxtr = 20;

mdls = NaN(length(lags),length(sbjs));
opt_lmb = NaN(maxtr,length(sbjs));
dpr = NaN(length(sbjs),1);
SNR = NaN(length(sbjs),1);
test_r = NaN(length(sbjs),maxtr);
circ_test_r = NaN(length(sbjs),100);
for s = 1:length(sbjs)
    % load data
    resfl = sprintf('Sbj%d_frwdtrf',sbjs(s));
    d = load([respth resfl]);
    % get the resulting model
    mdls(:,s) = mean(d.all_mdl,2); % average across leave-one-trial-out iterations
    opt_lmb(:,s) = d.opt_lmb; % optimal lambdas for each iteration
    % compute d-prime prediction accuracy
    mndiff = mean(d.test_r)-mean(d.circ_test_r);
    stdiff = sqrt(0.5*(var(d.test_r)+var(d.circ_test_r)));
    dpr(s) = mndiff/stdiff;
    SNR(s) = d.SNR;
    test_r(s,:) = d.test_r;
    circ_test_r(s,:) = d.circ_test_r;
end

% Sort d-prime descending
[srtdpr,dpr_idx] = sort(dpr,'ascend');

% Plot d-prime and prediction accuracies
figure
set(gcf,'Position',[100 100 900 400]);
hold on
yyaxis left
plot(dpr(dpr_idx),'bo','MarkerSize',14,'LineWidth',1.5);
plot(sbjs([1 end]),[2 2],'c--','LineWidth',1.5);
ylabel('d-prime prediction accuracy');
yyaxis right
md_r = median(test_r,2);
uq_r = quantile(test_r,0.75,2);
lq_r = quantile(test_r,0.25,2);
md_nr = median(circ_test_r,2);
uq_nr = quantile(circ_test_r,0.75,2);
lq_nr = quantile(circ_test_r,0.25,2);
errorbar(1:19,md_r(dpr_idx),md_r(dpr_idx)-lq_r(dpr_idx),uq_r(dpr_idx)-md_r(dpr_idx),'r.','MarkerSize',16);
errorbar(1:19,md_nr(dpr_idx),md_nr(dpr_idx)-lq_nr(dpr_idx),uq_nr(dpr_idx)-md_nr(dpr_idx),'kx','MarkerSize',16);
ylabel('Prediction accuracy');
set(gca,'FontSize',14);
xlabel('Subject');
exportgraphics(gca,'fig/AllSbjPerf_dpr_r.pdf','ContentType','vector')

% Plot each model
% EEGplot(mdls(:,dpr_idx(10:end))*3,Fs,1,lags/Fs*1000);
[yt,h] = plot_stacked_traces(lags/Fs*1000,mdls(:,dpr_idx),dpr(dpr_idx));
set(gcf,'Position',[360,25,400,675]);
set(gca,'YTickLabel',sbjs,'XLim',[tmin tmax]);
xlabel('Delay (ms)');
ylabel('Subject');
exportgraphics(gca,'fig/NatSpeech_trfbydpr.pdf','ContentType','vector');