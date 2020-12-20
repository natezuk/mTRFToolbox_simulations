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
for s = 1:length(sbjs)
    % load data
    resfl = sprintf('Subject%d_frwdtrf',sbjs(s));
    d = load([respth resfl]);
    % get the resulting model
    mdls(:,s) = mean(d.all_mdl,2); % average across leave-one-trial-out iterations
    opt_lmb(:,s) = d.opt_lmb; % optimal lambdas for each iteration
    % compute d-prime prediction accuracy
    mndiff = mean(d.test_r)-mean(d.circ_test_r);
    stdiff = sqrt(0.5*(var(d.test_r)+var(d.circ_test_r)));
    dpr(s) = mndiff/stdiff;
end

% Sort d-prime descending
[srtdpr,dpr_idx] = sort(dpr,'ascend');

% Plot d-prime
figure
hold on
plot(dpr(dpr_idx),'bo','MarkerSize',14,'LineWidth',1.5);
plot(sbjs([1 end]),[2 2],'k--','LineWidth',1.5);
set(gca,'FontSize',14);
xlabel('Subject');
ylabel('d-prime prediction accuracy');

% Plot each model
% EEGplot(mdls(:,dpr_idx(10:end))*3,Fs,1,lags/Fs*1000);
[yt,h] = plot_stacked_traces(lags/Fs*1000,mdls(:,dpr_idx),dpr(dpr_idx));
set(gcf,'Position',[360,25,400,675]);
set(gca,'YTickLabel',sbjs,'XLim',[tmin tmax]);
xlabel('Delay (ms)');
ylabel('Subject');