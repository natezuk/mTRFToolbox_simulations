% Evoked response model and phase
erp_dur = 500; % duration of erp, in ms
fs = 512;
ferp = 6; % in Hz
tauerp = 50;
terp = 0:1/fs:erp_dur/1000;
erp = sin(2*pi*terp*ferp).*exp(-terp/(tauerp/1000));

% Plot the original ERP
figure
plot(terp,erp,'k','LineWidth',2);
xlabel('Time (ms)');
ylabel('ERP');

% Plot the signal at each event rate, overlaid by a matched sinusoid from
% the fourier transform
evt_rate = 1:8;
stim_dur = 2;
t = (0:stim_dur*fs-1)/fs;
rate = NaN(length(evt_rate),1);
phase = NaN(length(evt_rate),1);
figure;
for ii = 1:length(evt_rate)
    evt = zeros(fs*stim_dur,1);
    evt_idx = round((0:1/evt_rate(ii):stim_dur-1/fs)*fs)+1;
    evt(evt_idx) = 1;
    % create the response with regular erps
    y = conv(evt,erp);
    y = y(1:stim_dur*fs);
    % get the phase at the event rate
    Y = fft(y);
    frq = (0:length(Y)-1)/length(Y)*fs;
    % in the spectrum, identify the closest index to the event rate
    df = abs(frq-evt_rate(ii));
    rate_idx = find(df==min(df),1,'first');
    rate(ii) = frq(rate_idx);
    phase(ii) = angle(Y(rate_idx)); % phase of sinusoid
    % Plotting
    subplot(length(evt_rate),1,ii);
    hold on
    sin_fit = real(exp(2*pi*1i*t*rate(ii) + phase(ii)*1i));
        % the sinusoid to plot is the real component of the complex
        % sinusoid at the event rate
    plot(t,sin_fit,'r');
    plot(t,y/max(y),'k'); % adjust so max magnitude is 1, comparable to the sinusoid
    plot(t(evt_idx),sin_fit(evt_idx),'r.','MarkerSize',16);
        % place dots on the sinusoid at the event onset times
    xlabel('Time (s)');
    title(sprintf('Rate = %d Hz',evt_rate(ii)));
end

% Plot phase rel rate
figure;
plot(rate,rad2deg(phase),'b');
set(gca,'FontSize',14);
xlabel('Event rate');
ylabel('Phase (degrees)');