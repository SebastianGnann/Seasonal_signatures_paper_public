function [FFT_max_F, FFT_max_Q] = calcFourierSpectrum(F,Q,w)
%PLOTFOURIERSPECTRUM 
%
% INPUT
% F: forcing, typically approximated by P - PET [mm]
% Q: streamflow [mm]
% w: angular frequency [1/days]
%
% OUTPUT 
% fft_max_F: maximum Fourier mode of forcing
% fft_max_Q: maximum Fourier mode of streamflow
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk


%% Fourier transform
N = length(Q(:,2)); % length of time series
fs = 1; % 1d sample frequency
% f = (0:N-1)*fs/N; % frequency range 1/d

% forcing
fft_F = fft(F(:,2)); % compute Fourier transform
P2_F = abs(fft_F/N); % compute Fourier transform
P1_F = P2_F(2:floor(N/2)+1); % compute one-sided spectrum
P1_F(2:end-1) = 2*P1_F(2:end-1);

% streamflow
fft_Q = fft(Q(:,2)); % compute Fourier transform
P2_Q = abs(fft_Q/N); % compute Fourier transform
P1_Q = P2_Q(2:floor(N/2)+1); % compute one-sided spectrum
P1_Q(2:end-1) = 2*P1_Q(2:end-1);

f = fs*(1:floor(N/2))/N; % define frequency domain
t = 1./f'/365; % and asociated time domain in years

[~,ind] = max(P1_F);
FFT_max_F = t(ind);
[~,ind] = max(P1_Q);
FFT_max_Q = t(ind);

end

