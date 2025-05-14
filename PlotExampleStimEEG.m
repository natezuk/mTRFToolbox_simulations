% Plot example stimulus and EEG data
% Expected that the data has been loaded with Test_speechdataorig_CND.m

t_range = [4 6]; % in s
t = (0:6-1/Fs)/Fs;
idx = find(t>=t_range(1)&t<t_range(2));

figure
set(gcf,'Position',[100 100 1000 600])
% stim1
subplot(2,2,1)
plot(t(idx),stims{1},'k');
xlabel('Time (s)')
ylabel('Stim 1')
% eeg1
subplot(2,2,2)
plot(t(idx),EEGs{1}(:,chan),'b');
xlabel('Time (s)')
ylabel('EEG 1')
% stim1
subplot(2,2,3)
plot(t(idx),stims{end},'k');
xlabel('Time (s)')
ylabel('Stim end')
% eeg1
subplot(2,2,4)
plot(t(idx),EEGs{end}(:,chan),'b');
xlabel('Time (s)')
ylabel('EEG end')

exportgraphics(gca,sprintf('ExmpStimEEG_sbj%d.pdf',sbj),'ContentType','vector');