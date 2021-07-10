function y = bandpass_noise(modfreq,dur,Fs,ph)
% Create bandpass filtered noise (using a brick frequency filter)
% Nate Zuk (2020) (based on my code from the Subcortical-beats project)
% Inputs:
% - modfreq = array of two values for the lower and upper cutoff of the
% bandpass filter, or one value for a single modulation frequency
% - dur = duration of the signal (in s)
% - Fs = sampling frequency of the signal (in Hz)
% - (Optional) ph = if one modulation frequency is specified, set the phase
% here (in periods, default = 0)

if nargin<4, ph = 0; end

% Create the modulation
if length(modfreq)==1,
    y = 0.5*(1+cos(2*pi*t*modfreq+2*pi*ph))'; % modulation
elseif length(modfreq)==2, % filter between those two frequencies for the envelope
    rndenv = randn(dur*Fs,1);
    RNDENV = fft(rndenv);
    f = (0:length(RNDENV)-1)/length(RNDENV)*Fs;
    idx = (f>=modfreq(1)&f<=modfreq(2))|(f<=Fs-modfreq(1)&f>=Fs-modfreq(2));
        % need to filter evenly on both sides of Fs/2 for real signals
    RNDENV(~idx) = 0; % remove those frequencies
    y = real(ifft(RNDENV));
%     rndenv = (rndenv-mean(rndenv))/range(rndenv)*2; % remove the mean and set to a range = 2
%     y = 0.5*(-min(rndenv)+rndenv); % set to a range of 0 to 1
end