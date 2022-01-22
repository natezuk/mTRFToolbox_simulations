% Analyze effect of fitting SVM and testing / training on shuffled data

nclass = 4;
nchan = 57;
ntr = 100;
snr = 0.1;

% Assign weightings corresponding to the 4 different classes
class_weights = randn(nchan,nclass);

% Generate trial data for each class, introducing noise
data = NaN(ntr,nchan,nclass);
for c = 1:nclass
    % Note: I can use 1/SNR to get the noise magnitude because the signal
    % std is, on average, 1 because it was generated with randn
    data(:,:,c) = ones(ntr,1)*class_weights(:,c)' + (1/snr)*randn(ntr,nchan);
end
% Rearrange the data so the first dimension is ntr x nclass
% data = permute(data,[1 3 2]);
% data = reshape(data,[ntr*nclass nchan]);
% lbls = repelem([1 2 3 4]',ntr,1);

% Fit the model
nfolds = 3;
% niter = 10;
truelbls = repelem([1 2 3 4]',ntr,1);
% Identify the starting and ending trials for each fold (start = folds(n)+1, end
% = folds(n+1));
folds = round(linspace(0,ntr,nfolds+1));
prdL = NaN(ntr*nclass,1);
%%% Cut repeating with 10 iterations, because the folds are always the
%%% same...
disp('Fitting true model...')
for n = 1:nfolds
    test_trs = folds(n)+1:folds(n+1);
    train_trs = setxor(1:ntr,test_trs);
    % Rearrange the training and testing data so it's ntr x nclass
    [trainD,trainL] = rearrange_data(data(train_trs,:,:));
    [testD,testL] = rearrange_data(data(test_trs,:,:));
    svm_mdl = fitcecoc(trainD,trainL);
    % Get the correct indexes in prdL for the testing data
    test_fullarr_idx = test_trs'*ones(1,nclass) + ones(length(test_trs),1)*((0:nclass-1)*ntr);
    test_fullarr_idx = reshape(test_fullarr_idx,[length(test_trs)*nclass, 1]);
    prdL(test_fullarr_idx) = predict(svm_mdl,testD);
end
acc = sum(prdL==truelbls)/(ntr*nclass);

%% Estimating null distribution
nnull = 100;
% Option 1) compute accuracy using predicted labels and shuffled true
% labels
nullacc_shufflbl = NaN(nnull,1);
disp('Calculating null distribution from testing predictions rel shuffled labels...');
for jj = 1:nnull
    shuffidx = randperm(ntr*nclass);
    nullacc_shufflbl(jj) = sum(prdL==truelbls(shuffidx))/(ntr*nclass);
end

%% Option 2) refit the model to shuffled labels, then calculate performance
nullacc_refit = NaN(nnull,1);
nullprdL = NaN(ntr*nclass,nnull);
disp('Calculating null distribution by refitting to shuffled labels...');
for jj = 1:nnull
    fprintf('*Null %d/%d*\n',jj,nnull);
    shuffidx = randperm(ntr*nclass);
    shuffL = truelbls(shuffidx);
    for n = 1:nfolds
        test_trs = folds(n)+1:folds(n+1);
        train_trs = setxor(1:ntr,test_trs);
        % Rearrange the training and testing data so it's ntr x nclass
        trainD = rearrange_data(data(train_trs,:,:));
        testD = rearrange_data(data(test_trs,:,:));
        % Get the correct indexes in prdL for the testing data
        test_fullarr_idx = test_trs'*ones(1,nclass) + ones(length(test_trs),1)*((0:nclass-1)*ntr);
        test_fullarr_idx = reshape(test_fullarr_idx,[length(test_trs)*nclass, 1]);
        train_fullarr_idx = train_trs'*ones(1,nclass) + ones(length(train_trs),1)*((0:nclass-1)*ntr);
        train_fullarr_idx = reshape(train_fullarr_idx,[length(train_trs)*nclass, 1]);
        % Reassign training and testing labels based on the shuffling
        testL = shuffL(test_fullarr_idx);
        trainL = shuffL(train_fullarr_idx);
        % Fit the model
        svm_mdl = fitcecoc(trainD,trainL);
        nullprdL(test_fullarr_idx,jj) = predict(svm_mdl,testD);
    end
    nullacc_refit(jj) = sum(nullprdL(:,jj)==shuffL)/(ntr*nclass);
end

%% Plotting
figure
hold on
% plot chance level (25%)
plot([0 3],[0.25 0.25],'k--');
% plot true accuracy
plot([0 3],[acc acc],'r--','LineWidth',2);
% Plot the null distribution for testing on shuffled labels
md_shufflbl = median(nullacc_shufflbl);
uq_shufflbl = quantile(nullacc_shufflbl,0.95);
lq_shufflbl = quantile(nullacc_shufflbl,0.05);
errorbar(1,md_shufflbl,md_shufflbl-lq_shufflbl,uq_shufflbl-md_shufflbl,'ko','MarkerSize',12,'LineWidth',2);
% Plot the null distribution for training on shuffled labels
md_refit = median(nullacc_refit);
uq_refit = quantile(nullacc_refit,0.95);
lq_refit = quantile(nullacc_refit,0.05);
errorbar(2,md_refit,md_refit-lq_refit,uq_refit-md_refit,'ko','MarkerSize',12,'LineWidth',2);
set(gca,'FontSize',14','XTick',1:2,'XTickLabel',{'Permute after test','Permute before test'},...
    'XTickLabelRotation',45);
ylabel('Proportion correct');

%% Functions %%
function [rd,l] = rearrange_data(d)
% rearrange data and assign labels
nclass = size(d,3);
nchan = size(d,2);
ntr = size(d,1);
rd = permute(d,[1 3 2]);
rd = reshape(rd,[ntr*nclass nchan]);
l = repelem([1 2 3 4]',ntr);
end