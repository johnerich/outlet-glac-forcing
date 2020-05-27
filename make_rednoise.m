%%% Script to make forcing timeseries with varying degrees of persistence.
%%% Used for outlet glacier forcing in Christian et al. (Cryosphere, 2020).
%%% Persistence can be an AR-1 type (one "memory" timescale) 
%%% or power-law type. Method to generate noise follows that of Roe and
%%% Baker (J. Glac., 2016) and Percival (J. Clim., 2001). The same set of
%%% random phases will be used, so that anomalies can be compared in the
%%% time domain. JEC May 2020.

% Assumed that dt = 1 year always.

clear all

% time domain:
N = 1e5; % number of time steps of final timeseries
dt = 1;  % script assumes this is 1, no guarantees it is robust if changed!
% frequency domain:
df = 1/(dt*N);
f0 = 0.5/dt;
f1 = [0:df:f0];

% persistence parameters
beta = 0.5; % exponent for power law noise (slope of spectrum in log-log space)
r = 0.95; % lag-1 autocorrelation for AR1 noise
tau = dt/(1-r); % "memory" of AR1 process
p0 = 1;

% power spectra (analytical)
pb = sqrt(p0*(f0./f1(2:end)).^beta); 
pr = sqrt(p0./(1 + r^2 - 2*r*cos(2*pi*dt*f1(2:end))));


% random phases - must be conjugate symmetric
phase = zeros(1,N);
pos_freq = pb(2:ceil(N/2)); % does NOT include nyquist or DC
ranphase = rand(size(pos_freq)); % note this is uniformly, not gaussian distributed from 0 to 1
phase_half = 1i*2*pi*ranphase;
phase(2:length(pos_freq)+1)=phase_half;
phase(N-length(pos_freq)+1:end) = conj(flip(phase_half));

pb2 = [pb,flip(pb(1:floor(N/2)))];
pr2 = [pr,flip(pr(1:floor(N/2)))];

% combine phase and power spectra:
pbrand = pb2.*exp(phase);
prrand = pr2.*exp(phase);
pwhite = exp(phase);

% back out timeseries from power spectrum, and normalize to unit variance
% for whole timeseries.
pl = ifft(pbrand,'symmetric');
pl = pl-mean(pl);
pl = pl/std(pl);

ar1 = ifft(prrand,'symmetric');
ar1 = ar1-mean(ar1);
ar1 = ar1/std(ar1);

white = ifft(pwhite,'symmetric');
white = white - mean(white);
white = white/std(white);

%% Analyze synthetic timeseries

% Autocorrelation:
maxlag = 1e2;
[acf_white,lags_dt] = xcorr(detrend(white),maxlag/dt,'coeff');
[acf_pl,lags_dt] = xcorr(detrend(pl),maxlag/dt,'coeff');
[acf_ar1,lags_dt] = xcorr(detrend(ar1),maxlag/dt,'coeff');

% Power spectrum via Welch's method:
trans_pl = fft(pl);
trans_ar1 = fft(ar1);
trans_white = fft(white);

window = 8;
power_white = abs(trans_white(1:length(f1))).^2;
[P_white,f_white] = pwelch(white,floor(N/window),[],N,1/dt,'onesided');
power_pl = abs(trans_pl(1:length(f1))).^2;
[P_pl,f_pl] = pwelch(pl,floor(N/window),[],N,1/dt,'onesided');
power_ar1 = abs(trans_ar1(1:length(f1))).^2;
[P_ar1,f_ar1] = pwelch(ar1,floor(N/window),[],N,1/dt,'onesided');

%% Plotting
figure(1);
% plot in time domain:
subplot 211; hold on
offset = 3;
plot(white,'k','linewidth',1)
plot(pl + offset,'b','linewidth',1)
plot(ar1 - offset,'r','linewidth',1)
xlim([0 500])
legend('white noise','power-law noise','AR-1 noise')
title('Anomalies in time domain')
xlabel('Time (yrs)')
ylabel('Anomalies (offset, arbitrary units)')

% Autocorrelation function:
subplot 223; hold on; grid on
plot(lags_dt,acf_white,'k','linewidth',1)
plot(lags_dt,acf_pl,'b','linewidth',1)
plot(lags_dt,acf_ar1,'r','linewidth',1)
xlim([0 30])
ylim([-0.1 1])
title('Autocorrelation function')
ylabel('Autocorr. coeff.')
xlabel('Lag (yr)')

% Power spectrum:
subplot 224; 
loglog(f1,P_white,'k','linewidth',1); hold on; grid on
loglog(f1,P_pl,'b','linewidth',1); 

loglog(f1,P_ar1,'r','linewidth',1); 

title('Power spectral density')
xlabel('Frequency (1/yrs)')
ylabel('Power (arbitrary units)')