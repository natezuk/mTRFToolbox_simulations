function flt_eeg = eeg_bandpass(eeg,Fs,varargin)
% Bandas filter the EEG using a 1st order Butterworth with cutoff
% frequencies at 1 and 15 Hz
% Use SOS-based zero-phase filtering
% Nate Zuk (2020)

N = 2; % order of filter
lowpass_cf = 15; % lowpass cutoff frequency
highpass_cf = 0.1; % highpass cutoff frequency

if ~isempty(varargin)
    for n = 2:2:length(varargin)
        eval([varargin{n-1} '=varargin{n};']);
    end
end

% Make the filter
D = fdesign.bandpass('N,F3dB1,F3dB2',N,highpass_cf,lowpass_cf,Fs);
flt = design(D,'butter');
flt_eeg = filtfilthd(flt,eeg);