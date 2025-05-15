% Plot the results from Test_TRFestimate_SNRnTr -- show effect of SNR and
% amount of data on d-prime and TRF estimate.
% Nate Zuk (2020)

addpath('../mTRF-toolbox/');
% U = userpath;
% addpath([U '/shadedErrorBar/']);

% Load the data
disp('Loading data...');
hf_d = load('TRFestimate_w_SNRnTr_highfreq.mat');
lf_d = load('TRFestimate_w_SNRnTr_delta.mat');

% It's assumed both datasets were calculated using the same condition
% parameters:
snr = hf_d.snr;
ntr = hf_d.ntr;
niter = size(hf_d.test_r{1},2);

% Set colors for both conditions
frqcnd_clrs = {'b','r'};

% Create the legend for each frequency region
frqrange_leg{1} = sprintf('%.2g - %.2g Hz',hf_d.freq_range(1),hf_d.freq_range(2));
frqrange_leg{2} = sprintf('%.2g - %.2g Hz',lf_d.freq_range(1),lf_d.freq_range(2));

% Get d-prime, r performance, and correlation with true model
[mdl_corr(:,:,:,1),dpr(:,:,:,1),mn_r(:,:,:,1)] = compute_pred_acc_measures(hf_d);
[mdl_corr(:,:,:,2),dpr(:,:,:,2),mn_r(:,:,:,2)] = compute_pred_acc_measures(lf_d);

% Reshape so all datapoints (irrespective of SNR/nTr) are along one column
ndatapts = size(mn_r,1)*size(mn_r,2)*size(mn_r,3); % get the total number of data points in each results file
R = reshape(mn_r,[ndatapts,2]);
DPR = reshape(dpr,[ndatapts,2]);
MDL_CORR = reshape(mdl_corr,[ndatapts,2]);

% Plot model correlation vs r values
figure
hold on
plot(R(:,1),MDL_CORR(:,1),'.','Color',frqcnd_clrs{1},'MarkerSize',12);
plot(R(:,2),MDL_CORR(:,2),'.','Color',frqcnd_clrs{2},'MarkerSize',12);
set(gca,'FontSize',14);
xlabel('Prediction accuracy (Pearson''s r)');
ylabel('Correlation with true TRF');
legend(frqrange_leg);
exportgraphics(gca,'fig/RvsTRFcorr_alldata.pdf','ContentType','vector');

% Collapse all conditions onto one plot, and correlation with true TRF for
% various d-prime bins
plot_quantiles = [0.1 0.9];
nbins = 20;
mdlcorr_md = NaN(nbins,2);
mdlcorr_qnt = NaN(nbins,2,2);
[mdlcorr_md(:,1),mdlcorr_qnt(:,:,1),hf_dpr_bins,hf_nvals] = distr_by_xbins(DPR(:,1),MDL_CORR(:,1),...
    plot_quantiles,nbins);
[mdlcorr_md(:,2),mdlcorr_qnt(:,:,2),lf_dpr_bins,lf_nvals] = distr_by_xbins(DPR(:,2),MDL_CORR(:,2),...
    plot_quantiles,nbins);
figure
set(gcf,'Position',[100 100 650 500]);
pl_mds = NaN(2,1);
hold on
cntr_hf = hf_dpr_bins(1:end-1)+diff(hf_dpr_bins)/2;
pl_mds(1) = plot(cntr_hf,mdlcorr_md(:,1),'Color',frqcnd_clrs{1},'LineWidth',2);
plot(cntr_hf,mdlcorr_qnt(:,1,1),'--','Color',frqcnd_clrs{1},'LineWidth',1.5);
plot(cntr_hf,mdlcorr_qnt(:,2,1),'--','Color',frqcnd_clrs{1},'LineWidth',1.5);
% errorbar(cntr_hf,mdlcorr_md(:,1),mdlcorr_md(:,1)-mdlcorr_qnt(:,1,1),mdlcorr_qnt(:,2,1)-mdlcorr_md(:,1),...
%     '.','Color',frqcnd_clrs{1},'MarkerSize',18,'LineWidth',1.5);
cntr_lf = lf_dpr_bins(1:end-1)+diff(lf_dpr_bins)/2;
pl_mds(2) = plot(cntr_lf,mdlcorr_md(:,2),'Color',frqcnd_clrs{2},'LineWidth',2);
plot(cntr_lf,mdlcorr_qnt(:,1,2),'--','Color',frqcnd_clrs{2},'LineWidth',1.5);
plot(cntr_lf,mdlcorr_qnt(:,2,2),'--','Color',frqcnd_clrs{2},'LineWidth',1.5);
% errorbar(cntr_lf,mdlcorr_md(:,2),mdlcorr_md(:,2)-mdlcorr_qnt(:,1,2),mdlcorr_qnt(:,2,2)-mdlcorr_md(:,2),...
%     '.','Color',frqcnd_clrs{2},'MarkerSize',18,'LineWidth',1.5);
set(gca,'FontSize',14);
xlabel('d-prime');
ylabel(sprintf('Correlation with true TRF (median [%d%% %d%%])',plot_quantiles(1)*100,plot_quantiles(2)*100));
legend(pl_mds,frqrange_leg,'Location','southeast');
exportgraphics(gca,'fig/DprvsTRFcorr.pdf','ContentType','vector');

% Collapse all conditions onto one plot, and correlation with true TRF for
% various Pearson's r bins
mdlcorr_r_md = NaN(nbins,2);
mdlcorr_r_qnt = NaN(nbins,2,2);
[mdlcorr_r_md(:,1),mdlcorr_r_qnt(:,:,1),hf_r_bins] = distr_by_xbins(R(:,1),MDL_CORR(:,1),...
    plot_quantiles,nbins);
[mdlcorr_r_md(:,2),mdlcorr_r_qnt(:,:,2),lf_r_bins] = distr_by_xbins(R(:,2),MDL_CORR(:,2),...
    plot_quantiles,nbins);
figure
set(gcf,'Position',[100 100 650 500]);
pl_mds = NaN(2,1);
hold on
cntr_hf = hf_r_bins(1:end-1)+diff(hf_r_bins)/2;
pl_mds(1) = plot(cntr_hf,mdlcorr_r_md(:,1),'Color',frqcnd_clrs{1},'LineWidth',2);
plot(cntr_hf,mdlcorr_r_qnt(:,1,1),'--','Color',frqcnd_clrs{1},'LineWidth',1.5);
plot(cntr_hf,mdlcorr_r_qnt(:,2,1),'--','Color',frqcnd_clrs{1},'LineWidth',1.5);
% errorbar(cntr_hf,mdlcorr_md(:,1),mdlcorr_md(:,1)-mdlcorr_qnt(:,1,1),mdlcorr_qnt(:,2,1)-mdlcorr_md(:,1),...
%     '.','Color',frqcnd_clrs{1},'MarkerSize',18,'LineWidth',1.5);
cntr_lf = lf_r_bins(1:end-1)+diff(lf_r_bins)/2;
pl_mds(2) = plot(cntr_lf,mdlcorr_r_md(:,2),'Color',frqcnd_clrs{2},'LineWidth',2);
plot(cntr_lf,mdlcorr_r_qnt(:,1,2),'--','Color',frqcnd_clrs{2},'LineWidth',1.5);
plot(cntr_lf,mdlcorr_r_qnt(:,2,2),'--','Color',frqcnd_clrs{2},'LineWidth',1.5);
% errorbar(cntr_lf,mdlcorr_md(:,2),mdlcorr_md(:,2)-mdlcorr_qnt(:,1,2),mdlcorr_qnt(:,2,2)-mdlcorr_md(:,2),...
%     '.','Color',frqcnd_clrs{2},'MarkerSize',18,'LineWidth',1.5);
set(gca,'FontSize',14);
xlabel('Pearson''s r');
ylabel(sprintf('Correlation with true TRF (median [%d%% %d%%])',plot_quantiles(1)*100,plot_quantiles(2)*100));
legend(pl_mds,frqrange_leg,'Location','southeast');
exportgraphics(gca,'fig/RvsTRFcorr.pdf','ContentType','vector');