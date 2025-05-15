% Plot the null distribution and true r values for a subject
sbj = 10;
svpth = 'ns_env_res/';

% load the data
load([svpth sprintf('Sbj%d_frwdtrf',sbj)]);

edges = -0.05:0.005:0.05;
centers = edges(1:end-1)+diff(edges)/2;
h_null = histcounts(circ_test_r,edges);

figure
set(gcf,'Position',[100 100 800 400]);
hold on
bar(centers,h_null/sum(h_null),1,'k');
stem(test_r,ones(length(test_r),1)*0.3,'r','LineWidth',1.5);
set(gca,'FontSize',14);
xlabel('Prediction accuracy (Pearson''s r)');
ylabel(sprintf('Proportion\nnull iterations'))
exportgraphics(gca,sprintf('fig/ExmpNullHist_wTrue_Sbj%d.pdf',sbj),'ContentType','vector');