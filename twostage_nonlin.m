% %%% Nonlinear version of Alex Robel's two-equation marine-terminating glacier model (see Robel, 
% Roe, Haseloff JGR 2018). Several canned geometries and initial conditions are provided. 

% Idealized forcings are easy to play with, and are set up mainly to run trends, steps, or white 
% noise in Precipitation or Grounding line flux coefficient. These are defined as fractional 
% anomalies. Model outputs L and H, as well as fluxes. JEC Dec. 2019.

clear all;
%%%%%%%%%%%%% default parameters; may be changed below with initial conditions
nts = 11e3;                 % years to run
dt=1;                       % time step; typically 1 yr
year = 3.15e7;              % seconds
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

Sbar = 0.5;                 % avg. interior surface mass balance rate (m/yr)
b0 = -100;                    % ice divide bed height (m)
bx = -2e-3;                  % bed slope (prograde if negative)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Some specific parameters and corresponding initial steady-state conditions (xg and h)
%%% Just uncomment desired block.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 1) Glaciers from Christian et al. 2020, TC
   %%% Glacier 1: linearized response times of Tf ~ 76, Ts ~ 2000 yrs
    xg = 184.75e3; h = 1413.2;  beta = (m+n+3)/(m+1); theta=0.7; b0 = -100; bx = -2e-3; Sbar = 0.5; 

   %%% Glacier 2: Tf ~ 56, Ts ~ 1160 yrs
   % xg = 212.022e3; h = 1569.22;  beta = (m+n+3)/(m+1); theta=0.75; b0 = 150; bx = -3e-3; Sbar = 0.6; C = 7.624e6;

   %%% Glacier 3: Tf ~ 144, Ts ~ 4600 yrs
   % xg = 700.47e3; h = 2813.56; beta = (m+n+3)/(m+1); theta=0.6; b0 = 100; bx = -1e-3; Sbar = 0.3;
   
%%% 2) Example geometry from Robel et al. 2018:
%xg = 446038; h = 2174; b0 = -100; Sbar = 0.3; beta = (m+n+3)/(m+1); theta = 0.6; bx = -1e-3;

%%% 3) a test case w Haseloff et al. buttressing (J glac, 2018). Assumed
%%% strong buttressing by an ice shelf with length Ls and width Ws in
%%% meters. Forcing can come via changes to Ls.
% xg = 420129; beta = 4; h = 2100.855; b0 = -100; bx = -4e-3; Ls = 40e3; Ws = 7.5e3; Sbar = 0.3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% define grounding-line-flux coefficients, according to assumed flux rule (exponent)
% for Schoof parameterization (default):
omega_bar = (A_glen*(rho_i*g)^(n+1)*(theta*(1-lambda^-1))^n*(4^n*C)^-1)^(1/(m+1));
% for Haseloff parameterization:
if beta == 4    % ice shelf mass loss is dominated by calving: 
    omega_bar = (n/2)^n*(n+1)^-(n+1)*(rho_i*g*(1-1/lambda))^n * A_glen*Ls^(-n)*Ws^(n+1); 
end
omega = omega_bar*ones(1,nts);

%%%%%%%%%%%%%%%% FORCINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% to generate fresh noise; also consider importing a file of the same
%%% anomalies to compare with another model, or by forcing type.
noiseS = randn(1,nts); noiseS = noiseS/std(noiseS); % interior SMB
noiseO = randn(1,nts); noiseO = noiseO/std(noiseO); % grounding line flux coeff. (omega)

%%%%%%%%%%%%%%%% define SMB forcing
% tweak these:
sigS = 0.0;             % for noise (fractional)
deltaS = -0.2;          % for step/trend (fractional)
startyrS = 1000;        % year that trend or step starts
endyrS = startyrS+1;    % year that trend or step ends

% no need to tweak:
anom = linspace(0,1,endyrS-startyrS);
S = deltaS*Sbar*[zeros(1,startyrS),anom,ones(1,nts-endyrS)] + Sbar + sigS*Sbar*noiseS;
S = S/3.15e7;

%%%%%%%%%%%%%%%%% define GL flux forcing for Schoof parameterization:
% tweak these:
sigO = 0.0;             % for noise (fractional)
deltaO = 0.0;           % for step/trend (fractional)
startyrO = 1000;        % year that trend or step starts
endyrO = startyrO+1;    % year that trend or step ends

% no need to tweak:
anom = linspace(0,1,endyrO-startyrO);
omega = deltaO*omega_bar*[zeros(1,startyrO),anom,ones(1,nts-endyrO)] + omega_bar + sigO*omega_bar*noiseO;

%%%%%%%%%%%%%%%%% define ice-shelf forcing for Haseloff parameterization:
% n.b. there may be noise-induced drift if forcing Ls, even if no step change; see Robel et al. 2018
if beta == 4
    sigLs = 0.0; % fractional forcing
    sigO = 0.0; % to force omega directly
    deltaLs = -0.0; % for step/trend (fractional) 
    deltaO = 0.0; % to tweak grounding line flux directly
    startyrO = 1e3;
    endyrO = 1e3+1;
    
    anom = linspace(0,1,endyrO-startyrO);
    Lsp = deltaLs*Ls*[zeros(1,startyrO),anom,ones(1,nts-endyrO)] + Ls + sigLs*Ls*noiseO;  
    omega = (n/2)^n*(n+1)^-(n+1)*(rho_i*g*(1-1/lambda))^n * A_glen*Lsp.^(-n)*Ws^(n+1); 

    if sigO ~= 0 % that is, if forcing omega directly
        omega = deltaO*omega_bar*[zeros(1,startyrO),anom,ones(1,nts-endyrO)] + omega_bar + sigO*omega_bar*noiseO;
    end
    
    if deltaO ~= 0
        omega = deltaO*omega_bar*[zeros(1,startyrO),anom,ones(1,nts-endyrO)] + omega_bar + sigO*omega_bar*randn(1,nts);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% linearized response time estimates from Robel et al. 2018:
b = b0 + (bx*xg);              % bed elevation corresponding to initial condition
hg = -(rho_w/rho_i)*b;         % grounding line thickness
St = 1 + beta*lambda*bx*xg/hg; % stability parameter
Q = nu*h^alpha / xg^n ;        % interior flux
Q_g = omega_bar*(hg^beta);     % grounding line flux

Tf_approx = (hg/Sbar)/(alpha + gamma + 1 - St);     % approx. fast timescale
Ts_approx = -(h*hg)/(alpha*Tf_approx*Sbar^2*St);    % approx. slow timescale 
                                                                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% initialize arrays
L_out = zeros(1,nts);
H_out = zeros(1,nts);
S_out = zeros(1,nts);
Q_out = zeros(1,nts);
Qg_out = zeros(1,nts);

L1_out = zeros(1,nts); xg1 = xg;
LPp = zeros(1,nts) + xg;
Llp = zeros(1,nts) + xg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time loop: %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TIME LOOP
for t = 1:nts
    % update fluxes:
    b = b0 + (bx*xg);
    hg = -(rho_w/rho_i)*b;
    Q = nu*h^alpha / xg^n ;
    Q_g = omega(t)*(hg^beta);

    % update geometry:
    dh_dt = S(t) - (Q_g/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/hg;
    h = h + dh_dt*dt*year;
    xg = xg + dxg_dt*dt*year;
    
    % outputs
    L_out(t) = xg;
    H_out(t) = h;
    Qg_out(t) = Q_g;
    Q_out(t) = Q;
    S_out(t) = S(t)*xg;

end
%%%%%%%%%%%% PLOTTING

%%
figure(1); %clf
subplot 311; hold on
plot([1:nts]/1e3,L_out/1e3,'linewidth',1);
grid on
ylabel('L (km)')
xlabel('time (kyr)')
title('Length (ice divide to grounding line)')

subplot 312; hold on
plot([1:nts]/1e3,H_out,'linewidth',1); grid on
ylabel('H (m)')
xlabel('time (kyr)')
title('Interior thickness')

subplot 313; hold on 
% n.b. this is mainly for looking at step changes, will be very noisy for stochastic forcing
plot([1:nts]/1e3,S_out*3.15e7,'linewidth',1);
plot([1:nts]/1e3,Q_out*3.15e7,'linewidth',1)
plot([1:nts]/1e3,Qg_out*3.15e7,'linewidth',1); grid on
legend('S \times L','Q    (interior)','Qg  (grounding line)')
xlabel('time (kyr)')
ylabel('Flux (m^2 / yr)')
title('System fluxes')

