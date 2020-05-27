%%% This script creates a large ensemble of synthetic ~century-long
%%% "observations" to build up a null distribution for trend testing. It
%%% uses the 2-D vector version of the linearized two-stage glacier model
%%% (Robel et al., 2018 JGR) as implemented in Christian et al., 2020 TC.
%%% Here, it is set up to force the glacier with stochastic variability, which
%%% can be white or AR1 red noise. Most relevant outputs are "synth_out",
%%% which is the entire collection of synthetic "observations", and
%%% "deltaL", which is the linear trends in all those observations. These
%%% are the basis for Fig 6c in Christian et al., 2020 TC. 

%%% JEC May 2020

clear all; %close all;


dt = 1*3.15e7;  % typically 1 year (in seconds)
nts = 1e4; % full length of each synthetic; should be long enough to decorrelate, ~ 5*Ts
numits = 1e4; % number of synthetics to generate
winsize = 100; % numnber of years to save as observation

% initialize arrays
synth_out = zeros(numits,winsize);

%%%% Glacier Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% steady-state values from nonlinear model: "glacier 1" in Christian et al 2020, Tf = 76 yrs, Ts ~ 2000 yrs
Hbar = 1413.2;              % (m)
Lbar = 184.75e3;            % (m)
hg = 526.3;                 % thickness @ grounding line (m)

% other parameters for this glacier:
Sbar = 0.5/3.15e7;          % interior accumulation rate (m/s)
bx = -2e-3;                 % bed slope
n = 3;                      % creep exponent
m = 1/n;                    % sliding exponent
beta = (m+n+3)/(m+1);       % grounding line flux exponent (Schoof 2007)
alpha = 2*n + 1;
gamma = n;
theta=0.7;                  % buttressing parameter (between 0 and 1; see Schoof 2007)
g = 9.81;
rho_i = 917; rho_w = 1028;  % ice, water densities
lambda = rho_w/rho_i;       
A_glen = 4.22e-25;          % Nye-Glen coeff (Pa^-3 s^-1)
C = 7.624e6;                % Weertman coeff (Pa m^-1/3 s^1/3)

omega = (A_glen*(rho_i*g)^(n+1)*(theta*(1-lambda^-1))^n*(4^n*C)^-1)^(1/(m+1)); % grounding line flux coeff. (Schoof 2007)
Qg = omega*hg^beta; % grounding line flux

% linearized couplings between length and thickness 
ah = -Qg*alpha/(hg*Lbar);
al = Qg/Lbar^2*(1 + gamma*Hbar/hg + beta*lambda*bx*Lbar/hg*(1 - Hbar/hg));
bh = Qg*alpha/(Hbar*hg);
bl = Qg/hg*(beta*lambda*bx/hg - gamma/Lbar);

% re-define some combos to save space...
aa = 1/(1 - ah*dt);
ee = 1/(1 - bl*dt);
kk = 1/(1-al*bh*aa*ee*dt^2);

%% LOOP OVER SYNTHETICS
disp('generating synthetics')
for nn = 1:numits
    H_out = zeros(1,nts);
    L_out = zeros(1,nts);
    %%% new forcing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = 0.95; % lag-1 autocorrelation for AR1 noise; set to 0 for white noise
    tau_noise = dt/(1-r)/dt; % memory timescale of noise (yrs)

    if r == 0
        noise = randn(1,nts);
        noise = noise-mean(noise);
        noise = noise/std(noise);
    else % generates red noise with the specified memory
        p0 = 1;
        df = 1/(dt*nts);
        phase = zeros(1,nts);
        f0 = 0.5/dt;
        f1 = [0:df:f0];
        pr = sqrt(p0./(1 + r^2 - 2*r*cos(2*pi*dt*f1(2:end))));
        pos_freq = pr(2:ceil(nts/2)); % does NOT include nyquist or DC
        phase_half = 1i*2*pi*rand(size(pos_freq));
        phase(2:length(pos_freq)+1)=phase_half;
        phase(nts-length(pos_freq)+1:end) = conj(flip(phase_half));
        pr2 = [pr,flip(pr(1:floor(nts/2)))];
        prrand = pr2.*exp(phase);
        noise = ifft(prrand,'symmetric');
        noise = noise-mean(noise);
        noise = noise/std(noise);
    end

    % forcing GL flux:
    sigO = 0.2;             % fractional noise
    omegap = sigO*omega*noise;
    Qgp = omegap*hg^beta;

    % fractional noise in SMB.:
    sigS = 0.0;             % fractional noise
    Sp = sigS*Sbar*noise; 


    % solve vector AR equations:
    for ii = 2:nts
        % length:
        L_out(ii) = kk * ee * [L_out(ii-1) + aa*bh*dt*H_out(ii-1) + dt*(bh*dt*(Hbar/hg-1)/Lbar - 1/hg)*Qgp(ii) + bh*aa*dt^2*Sp(ii)];
        % thickness:
        H_out(ii) = kk * aa * [H_out(ii-1) + ee*al*dt*L_out(ii-1) + dt*((Hbar/hg-1)/Lbar - al*ee*dt/hg)*Qgp(ii) + dt*Sp(ii)];
    end

    synth_out(nn,:) = L_out(end-winsize+1:end);
    if nn/1000 == floor(nn/1000)
        disp(strcat('iteration__',num2str(nn),'/',num2str(numits)))
    end
end
%% 

[a b] = size(synth_out); % a should be # of synthetic runs, b the observational window

b = 50; % reduce observational period, if desired 
% (if comparing results for different obs intervals, do this instead of running entire loop again)
t_win = (1:b);

deltaL = zeros(a,1);

% Determine trends in synthetic observations
disp('finding trends')
for ii = 1:a
    if floor(ii/1000) == ii/1000;
        disp(strcat('iteration__',num2str(ii),'/',num2str(a)))
    end
    X = synth_out(ii,1:b);
    trend = polyfit(t_win,X,1); % best-fit linear slope
    trend = trend(1);
    trend = trend*b; % convert to length change
    deltaL(ii) = trend; % save delta L
end  

%% schematic figure that illustrates process:
figure(1);
plotlim = nts; % in years

subplot 221; hold on
plot([1:plotlim]/1e3,L_out(1:plotlim),'b')
plot([plotlim-b+1:plotlim]/1e3,L_out(plotlim-b+1:end),'r','linewidth',2)
legend('Full synthetic segment','retained for trend testing')
xlabel('Time (kyr)')
ylabel('L'' (m)')
title('One member of ensemble')

subplot 223; hold on
histogram(deltaL,'normalization','pdf')
xlabel(strcat('\DeltaL in observational period (',num2str(b),'-yr obs.)'))
title('Distribution of linear trends in ensemble')

subplot 122; hold on
plot(synth_out(end-20:end,end-b+1:end)')
xlabel('Time (yr)')
ylabel('L'' (m)')
title('Sample of ensemble used for trend testing')
