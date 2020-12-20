function rnd_x = randomize_phases(x)
% Take a signal x and create a new signal rnd_x with the same magnitude
% spectrum but randomized phases. But keep the phases relatively correlated
% across columns by adding the random phases to the original phases.
% Nate Zuk (2020)

X = fft(x);
% identify the number of points to synthesize
% npts_to_synth = ceil(length(x)/2)-1;
npts_to_synth = ceil((length(x)+1)/2)-1;
ncols = size(x,2); % the number of columns in x
% generate random phases
rnd_ph = (exp(2*pi*1i*rand(npts_to_synth,1)))*ones(1,ncols);
% recreate the frequency spectrum with randomized phases
% NZ (7-6-2020): If X has an even number of points, remove the last phase
% value when the phase values are inverted. Otherwise, keep it in.
invrt_ph = flipud(1./rnd_ph(1:end-mod(length(x)-1,2),:));
% RND_X = [X(1,:); ...
%     X(2:npts_to_synth+1,:).*rnd_ph;...
%     X(npts_to_synth+2:end,:).*flipud(1./invrt_ph)];
RND_X = [X(1,:); ...
    X(2:npts_to_synth+1,:).*rnd_ph;...
    X(npts_to_synth+2:end,:).*invrt_ph];
rnd_x = real(ifft(RND_X)); % get the real component of the signal