% Plot the results from Test_TRFestimate_SNRnTr -- show effect of SNR and
% amount of data on d-prime and TRF estimate.
% Nate Zuk (2020)

%addpath('../mTRF-toolbox/');
U = userpath;
addpath([U '/shadedErrorBar/']);

% Load the data
disp('Loading data...');
sim_cond = 'delta'; % suffix defining the frequency range for the simulation
load(sprintf('TRFestimate_w_SNRnTr_%s.mat',sim_cond));

% Compute the true trf
resp_t = (0:ceil(resp_dur/1000*Fs))/Fs;
% resp_frq = 6; % frequency of the response
% the trf is an exponentially decaying sinusoid
true_trf = sin(2*pi*resp_frq*resp_t).*exp(-resp_t/(resp_dur/1000/5));

% Compute the predicted trf lags
lags = floor(trf_range(1)/1000*Fs):ceil(trf_range(2)/1000*Fs);

% Number of iterations for each SNR/ntr pair
niter = size(test_r{1},2);

% Plot each SNR as a separate color
cmap = colormap('jet');
% create a legend for each snr color plotted
snr_leg = cell(length(snr),1);
snr_lbl = cell(length(snr),1); % label for axis ticks of image plots
for jj = 1:length(snr)
    snr_leg{jj} = sprintf('SNR = %d dB',20*log10(snr(jj))); 
    snr_lbl{jj} = sprintf('%d',20*log10(snr(jj)));
end

%% First, just plot the true model:
figure
set(gcf,'Position',[100 100,500,250]);
plot(resp_t*1000,true_trf,'LineWidth',2,'Color','k');
set(gca,'FontSize',14);
xlabel('Delay (ms)');
ylabel('True TRF');
exportgraphics(gca,sprintf('fig/TrueModel_%s.pdf',sim_cond),'ContentType','vector');

%% Plot the true model and the model estimates for each number of trials
figure
set(gcf,'Position',[45,315,1200,385]);
for kk = 1:length(ntr)
    subplot(ceil(length(ntr)/3),3,kk)
    hold on
    mdl_plt = NaN(length(snr)+1,1);
    mdl_plt(1) = plot(resp_t*1000,true_trf/rms(true_trf),'k--','LineWidth',2);
    for jj = 1:length(snr)
        clr_idx = get_color_idx(jj,length(snr),cmap);
        mdl_rms = repmat(rms(total_mdl{kk}(:,:,jj),1),[length(lags),1]);
            % normalize each estimated model by its rms, to account for scaling
            % issues
        mn_m = mean(total_mdl{kk}(:,:,jj)./mdl_rms,2);
        std_m = std(total_mdl{kk}(:,:,jj)./mdl_rms,[],2);
%         mn_m = mean(total_mdl{kk}(:,:,jj),2);
%         std_m = std(total_mdl{kk}(:,:,jj),[],2);
        pl = shadedErrorBar(lags/Fs*1000,mn_m,std_m/sqrt(niter),'lineProps',{'Color',cmap(clr_idx,:),'LineWidth',2});
        mdl_plt(jj+1) = pl.mainLine;
    end
    set(gca,'FontSize',14);
    xlabel('Delay (ms)');
    ylabel('TRF weight (RMS=1)');
    title(sprintf('%d trials',ntr(kk)));
end
legend(mdl_plt,[{'True model'}; snr_leg]);

%% Plot two images
% one of d-prime as a function of SNR and number of
% trials, and one of the correlation between true and prediction TRF as a
% function of SNR and number of trials
dpr = NaN(niter,length(snr),length(ntr));
mdl_corr = NaN(niter,length(snr),length(ntr));
mn_r = NaN(niter,length(snr),length(ntr));
for kk = 1:length(ntr)
    for jj = 1:length(snr)
        for ii = 1:niter
            % compute dprim
            mdiff = mean(test_r{kk}(:,ii,jj))-mean(circ_test_r{kk}(:,ii,jj));
            stdcmb = sqrt(0.5*(var(test_r{kk}(:,ii,jj))+var(circ_test_r{kk}(:,ii,jj))));
            dpr(ii,jj,kk) = mdiff/stdcmb;
            % get the indexes that correspond to the start and end of the
            % true trf
            trf_start_idx = find(resp_t(1)*Fs==lags);
            trf_end_idx = find(resp_t(end)*Fs==lags);
            mdl_corr(ii,jj,kk) = corr(total_mdl{kk}(trf_start_idx:trf_end_idx,ii,jj),true_trf');
            % compute the average r value across trials
            mn_r(ii,jj,kk) = mean(test_r{kk}(:,ii,jj));
        end
    end
end

figure
set(gcf,'Position',[45,315,1200,300]);

subplot(1,3,1); % ** trf correlation
imagesc(1:length(snr),1:length(ntr),squeeze(median(mdl_corr,1))');
colorbar;
caxis([0 1]);
set(gca,'FontSize',14,'XTick',1:length(snr),'XTickLabel',snr_lbl,'XTickLabelRotation',45,...
    'YTick',1:length(ntr),'YTickLabel',ntr);
xlabel('SNR (dB)');
ylabel('Number of trials');
title('Median correlation with true TRF');

subplot(1,3,2); % average r value
imagesc(1:length(snr),1:length(ntr),squeeze(median(mn_r,1))');
colorbar;
caxis([0 0.1]);
set(gca,'FontSize',14,'XTick',1:length(snr),'XTickLabel',snr_lbl,'XTickLabelRotation',45,...
    'YTick',1:length(ntr),'YTickLabel',ntr);
xlabel('SNR (dB)');
ylabel('Number of trials');
title('Median prediction accuracy');

subplot(1,3,3); % ** d-prime plot
imagesc(1:length(snr),1:length(ntr),squeeze(median(dpr,1))');
colorbar;
caxis([0 5]);
set(gca,'FontSize',14,'XTick',1:length(snr),'XTickLabel',snr_lbl,'XTickLabelRotation',45,...
    'YTick',1:length(ntr),'YTickLabel',ntr);
xlabel('SNR (dB)');
ylabel('Number of trials');
title('Median d-prime');

exportgraphics(gca,sprintf('fig/TRFcor_R_Dpr_%s.pdf',sim_cond),'ContentType','vector');

% Collapse all conditions onto one plot, and correlation with true TRF for
% various d-prime bins
DPR = reshape(dpr,[numel(dpr),1]);
MDL_CORR = reshape(mdl_corr,[numel(mdl_corr),1]);
figure
plot(DPR,MDL_CORR,'k.','MarkerSize',12);
set(gca,'FontSize',14);
xlabel('d-prime');
ylabel('Correlation with true TRF');

figure % plot each condition as a different color/size
set(gcf,'Position',[100 100 800 600]);
hold on
ntr_sizes = [10 14 18 22 26]; % marker sizes for each # trials
for kk = 1:length(ntr)
    for jj = 1:length(snr)
        clr_idx = get_color_idx(jj,length(snr),cmap);
        plot(dpr(:,jj,kk),mdl_corr(:,jj,kk),'.','Color',cmap(clr_idx,:),...
            'MarkerSize',ntr_sizes(kk)*1.5);
    end
end
set(gca,'FontSize',14,'XLim',[-2 6]);
xlabel('d-prime');
ylabel('Correlation with true TRF');
exportgraphics(gca,sprintf('fig/DprvsTRFcorr_snrcolor_ntrsize_%s.pdf',sim_cond),'ContentType','vector');

figure % plot each condition as a different color/size
set(gcf,'Position',[100 100 800 600]);
hold on
ntr_sizes = [10 14 18 22 26]; % marker sizes for each # trials
for kk = 1:length(ntr)
    for jj = 1:length(snr)
        clr_idx = get_color_idx(jj,length(snr),cmap);
        plot(dpr(:,jj,kk),atanh(mdl_corr(:,jj,kk)),'.','Color',cmap(clr_idx,:),...
            'MarkerSize',ntr_sizes(kk)*1.5);
    end
end
set(gca,'FontSize',14,'XLim',[-2 6]);
xlabel('d-prime');
ylabel('Correlation with true TRF (Fisher transform)');
exportgraphics(gca,sprintf('fig/DprvsTRFcorr_snrcolor_ntrsize-fisherz_%s.pdf',sim_cond),'ContentType','vector');

% Plot model correlation vs r values
R = reshape(mn_r,[numel(mn_r),1]);
figure
plot(R,MDL_CORR,'k.','MarkerSize',12);
set(gca,'FontSize',14);
xlabel('Prediction accuracy (Pearson''s r)');
ylabel('Correlation with true TRF');

% Plot a few examples of models
corr_examples = [0.1 0.6 0.999];
allr_idx = NaN(length(corr_examples),1);
figure
set(gcf,'Position',[360 1 500 700]);
for n = 1:length(corr_examples)
    subplot(length(corr_examples),1,n);
    % get the index for the closest DPR value to the 
    diffR = abs(MDL_CORR-corr_examples(n));
    allr_idx(n) = find(diffR==min(diffR),1,'last');
    % get the corresponding iter, snr, and ntr indexes for this simulation
    iter_idx = mod(allr_idx(n)-1,niter)+1;
    snr_idx = mod(floor((allr_idx(n)-1)/niter),length(snr))+1;
    ntr_idx = floor((allr_idx(n)-1)/(niter*length(snr)))+1;
    % get the model calculated for that simulation
    trf_start_idx = find(resp_t(1)*Fs==lags);
    trf_end_idx = find(resp_t(end)*Fs==lags);
    mdl = total_mdl{ntr_idx}(:,iter_idx,snr_idx);
    % plot the true TRF and the fit
    hold on
    plot(resp_t*1000,true_trf,'k--','LineWidth',3);
    plot(lags/Fs*1000,mdl/rms(mdl)*rms(true_trf),'b','LineWidth',2);
    set(gca,'FontSize',14,'YLim',[-0.5 0.7],'XLim',trf_range);
    xlabel('Delay (ms)');
    ylabel('TRF weight');
    title(sprintf('Iter = %d, SNR = %d, #trials = %d; dpr = %.3f, r = %.3f',...
        iter_idx,20*log10(snr(snr_idx)),ntr(ntr_idx),...
        DPR(allr_idx(n)),MDL_CORR(allr_idx(n))));
end

% Plot median and quantiles within each DPR bin, and plot dots for the
% example TRFs
plot_quantiles = [0.1 0.9];
[mdl_corr_md,mdl_corr_quantiles,dpr_bins,ndpr] = distr_by_xbins(DPR,MDL_CORR,plot_quantiles);
figure
cntr = dpr_bins(1:end-1)+diff(dpr_bins)/2;
hold on
plot(cntr,mdl_corr_md,'k','LineWidth',2);
plot(cntr,mdl_corr_quantiles,'k--','LineWidth',1.5);
plot(DPR(allr_idx),MDL_CORR(allr_idx),'r.','MarkerSize',32);
% errorbar(cntr,mdl_corr_md,mdl_corr_md-mdl_corr_quantiles(:,1),mdl_corr_quantiles(:,2)-mdl_corr_md,...
%     'k.','MarkerSize',18,'LineWidth',1.5);
set(gca,'FontSize',14);
xlabel('d-prime');
ylabel(sprintf('Correlation with true TRF (median [%d%% %d%%])',plot_quantiles(1)*100,plot_quantiles(2)*100));
title(sprintf('%.3g datapoints per bin',mean(ndpr)));