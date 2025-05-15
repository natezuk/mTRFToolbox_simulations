% Plot example stimulus and EEG data
% Expected that the data has been loaded with Test_speechdataorig_CND.m

t_range = [4 8]; % in s
t = (0:6*Fs-1)/Fs;
idx = find(t>=t_range(1)&t<t_range(2));

figure
set(gcf,'Position',[100 100 800 400])
% stim1
subplot(2,2,1)
plot(t(idx),stims{1}(idx),'k');
xlabel('Time (s)')
ylabel('Stim 1')
% eeg1
subplot(2,2,2)
plot(t(idx),EEGs{1}(idx),'b');
xlabel('Time (s)')
ylabel('EEG 1')
% stim1
subplot(2,2,3)
plot(t(idx),stims{end}(idx),'k');
xlabel('Time (s)')
ylabel('Stim end')
% eeg1
subplot(2,2,4)
plot(t(idx),EEGs{end}(idx),'b');
xlabel('Time (s)')
ylabel('EEG end')

exportgraphics(gcf,sprintf('fig/ExmpStimEEG_sbj%d.pdf',sbj),'ContentType','vector');