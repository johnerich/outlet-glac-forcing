%%% 2-D vector version of Alex Robel's linearized glacier model (see Robel, 
% Roe, Haseloff JGR 2018) as implemented in Christian et al. 2020, TC. 
% Linearization here includes perturbations in grounding-line flux. 
% Several canned geometries and initial conditions are provided. 

% Idealized forcings are easy to play with, and are set up mainly to run trends, steps, or white 
% noise in Precipitation or Grounding line flux coefficient. These are defined as fractional 
% anomalies. Model outputs L and H anomalies. JEC Dec. 2019.

clear all; %close all;

dt = 1*3.15e7;  % timestep, typically 1 year (in seconds)
nts = 1e4;      % number of time steps

% initialize output arrays:
H_out = zeros(1,nts);
L_out = zeros(1,nts);

%%%% Glacier Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% steady-state values obtained from nonlinear model: 
% defualt: "glacier 1" from Christian et al. 2020, TC
% Tf = 76 yrs, Ts ~ 2000 yrs

Hbar = 1413.2;              % (m)
Lbar = 184.75e3;            % (m)
hg = 526.3;                 % thickness @ grounding line (m)
Sbar = 0.5/3.15e7;          % interior accumulation rate (m/s)
bx = -2e-3;                 % bed slope (prograde if negative)

n = 3;                      % creep exponent
m = 1/n;                    % sliding exponent
alpha = 2*n + 1;
gamma = n;
beta = (m+n+3)/(m+1);       % Grounding line flux exponent; see Schoof 2007, JGR
theta = 0.7;                % Buttressing parameter for Schoof 2007 flux condition 
rho_i = 917; rho_w = 1028;  % densities
lambda = rho_w/rho_i;       % water to ice density
g = 9.81;                   % gravity
A_glen = 4.22e-25;          % Nye-Glen coeff (Pa^-3 s^-1)
C = 7.624e6;                % Weertman coeff (Pa m^-1/3 s^1/3)
nu = (rho_i*g/C)^n;


% Alternate parameters - uncomment to use: (or build your own)
   %%% Glacier 2: Tf ~ 56, Ts ~ 1160 yrs
% Lbar = 212.022e3; Hbar = 1569.22; hg = 544.9; theta=0.75; b0 = 150; bx = -3e-3; Pbar = 0.6/3.15e7; 

   %%% Glacier 3: Tf ~ 144, Ts ~ 4600 yrs
% Lbar = 700.47e3; Hbar = 2813.56; hg = 673.2; theta=0.6; b0 = 100; bx = -1e-3; Pbar = 0.3/3.15e7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

St = 1 + beta*lambda*bx*Lbar/hg; % stability parameter
Tf_approx = (hg/Sbar)/(alpha + gamma + 1 - St); % fast timescale
Ts_approx = -(Hbar*hg)/(alpha*Tf_approx*Sbar^2*St); % slow timescale
Tf_approx = Tf_approx/3.15e7; % yrs
Ts_approx = Ts_approx/3.15e7; % yrs

% define grounding-line flux according to Schoof 2007:
omega = (A_glen*(rho_i*g)^(n+1)*(theta*(1-lambda^-1))^n*(4^n*C)^-1)^(1/(m+1));
Qg = omega*hg^beta; % grounding line flux

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% linearized couplings between length and thickness; see supplement of Robel et al. 2018 
ah = -Qg*alpha/(hg*Lbar);
al = Qg/Lbar^2*(1 + gamma*Hbar/hg + beta*lambda*bx*Lbar/hg*(1 - Hbar/hg));
bh = Qg*alpha/(Hbar*hg);
bl = Qg/hg*(beta*lambda*bx/hg - gamma/Lbar);

% re-define some combos to save space...
aa = 1/(1 - ah*dt);
ee = 1/(1 - bl*dt);
kk = 1/(1-al*bh*aa*ee*dt^2);

%%%%%%%%%%%%%%%% FORCINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% to generate fresh noise; also consider importing a file of the same
%%% anomalies to compare with another model, or by forcing type.

noiseS = randn(1,nts); noiseS = noiseS/std(noiseS); % interior SMB.
noiseO = randn(1,nts); noiseO = noiseO/std(noiseO); % grounding line flux coeff. (omega)

%%%%%%%%%%%%%%%% define SMB forcing
% tweak these:
sigS = 0.1;             % for noise (fractional)
deltaS = -0.0;          % for step/trend (fractional)
startyrS = 1000;        % year that trend or step starts
endyrS = startyrS+1;    % year that trend or step ends

% no need to tweak:
anom = linspace(0,1,endyrS-startyrS);
Sp = deltaS*Sbar*[zeros(1,startyrS),anom,ones(1,nts-endyrS)] + sigS*Sbar*noiseS; % step/trend + noise

%%%%%%%%%%%%%%%%% define GL flux forcing for Schoof parameterization:
% tweak these:
sigO = 0.0;             % for noise (fractional)
deltaO = 0.0;           % for step/trend (fractional)
startyrO = 1000;        % year that trend or step starts
endyrO = startyrO+1;    % year that trend or step ends

% no need to tweak:
anom = linspace(0,1,endyrO-startyrO); 
omegap = deltaO*omega*[zeros(1,startyrO),anom,ones(1,nts-endyrO)] + sigO*omega*noiseO; % step/trend + noise
Qgp = omegap*hg^beta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve vector AR equations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 2:nts
    % length:
    L_out(ii) = kk * ee * [L_out(ii-1) + aa*bh*dt*H_out(ii-1) + dt*(bh*dt*(Hbar/hg-1)/Lbar - 1/hg)*Qgp(ii) + bh*aa*dt^2*Sp(ii)];
    % thickness:
    H_out(ii) = kk * aa * [H_out(ii-1) + ee*al*dt*L_out(ii-1) + dt*((Hbar/hg-1)/Lbar - al*ee*dt/hg)*Qgp(ii) + dt*Sp(ii)];
end

%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%

figure(1)
plotlim = 1e4; % in years
plotlim = min(plotlim,nts);

subplot 211; hold on
plot([1:plotlim]/1e3,H_out(1:plotlim),'b')
xlabel('Time (kyr)')
ylabel('H anomaly (m)')

subplot 212; hold on
plot([1:plotlim]/1e3,L_out(1:plotlim),'b')
xlabel('Time (kyr)')
ylabel('L anomaly (m)')
