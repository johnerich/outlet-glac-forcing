%%% 2-D vector version of linearized 2-stage model glacier model (Robel et
%%% al., 2018 JGR-ES), implemented in Christian et al., 2020 TC. This
%%% script is a rough template for playing around with idealized
%%% nonstationary climate changes (See Christian et al. 2020 fig. 7)
%%% As discussed there, the
%%% focus is not the magnitude of length variations, temperature variations
%%% etc., but rather on how the response times filter anomalies from
%%% different sources. So take all magnitudes with a grain of salt, they
%%% are not optimized here in any way. 

%%% You can play around with the cartoon climate in the accompanying
%%% script, or change relative magnitudes of ocean (omega) and interior (P)
%%% forcing in this script. Or, load in your favorite paleo reconstruction
%%% for a more realistic timeseries.

%%% JEC May 2020

clear all; %close all;

% generates sequence of anomalies; cartoon deglaciation, LIA, and
% anthropogenic warming. You can mess around with parameters in that
% script. The default that the glacier is linearized w respect to the
% stable Holocene, but there is a long spinup period to allow it to come
% into equilibrium before the cartoon deglaciation signal.

cartoon_climate4paper;

dt = 1*3.15e7;  % typically 1 year (in seconds)
nts = length(lgm_lia_anth); % # time steps

H_out = zeros(1,nts);
L_out = zeros(1,nts);

%%%% Glacier Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% steady-state values from nonlinear model: medum ice stream for AGU: Tf = 76 yrs, Ts ~ 2000 yrs
Hbar = 1413.2;              % (m)
Lbar = 184.75e3;            % (m)
hg = 526.3;                 % thickness @ grounding line (m)
Pbar = 0.5/3.15e7;          % interior accumulation rate (m/s)
bx = -2e-3;                 % bed slope

n = 3;                      % creep exponent
m = 1/n;                    % sliding exponent
beta = (m+n+3)/(m+1);       % grounding line flux exponent (Schoof 2007)
alpha = 2*n + 1;
gamma = n;

g = 9.81;
rho_i = 917; rho_w = 1028;  % ice, water densities
lambda = rho_w/rho_i;       
theta=0.7;                  % buttressing parameter (between 0 and 1; see Schoof 2007)
A_glen = 4.22e-25;          % Nye-Glen coeff (Pa^-3 s^-1)
C = 7.624e6;                % Weertman coeff (Pa m^-1/3 s^1/3)

omega = (A_glen*(rho_i*g)^(n+1)*(theta*(1-lambda^-1))^n*(4^n*C)^-1)^(1/(m+1)); % grounding line flux coeff.
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

%%% forcing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paleo = lgm_lia_anth;

anomP = paleo;
anomO = paleo;
%noiseO = noise(1:nts); % minus so that flux anomalies have same sign

% forcing GL flux:
scaleO = 0.1;             % scaling for anomalies

omegap = scaleO*omega*anomO; 
Qgp = omegap*hg^beta;

% forcing in SMB.:
scaleP = 0.0;             % scaling for anomalies
Pp = scaleP*Pbar*anomP; 


% solve vector AR equations:
for ii = 2:nts
    % length:
    L_out(ii) = kk * ee * [L_out(ii-1) + aa*bh*dt*H_out(ii-1) + dt*(bh*dt*(Hbar/hg-1)/Lbar - 1/hg)*Qgp(ii) + bh*aa*dt^2*Pp(ii)];
    % thickness:
    H_out(ii) = kk * aa * [H_out(ii-1) + ee*al*dt*L_out(ii-1) + dt*((Hbar/hg-1)/Lbar - al*ee*dt/hg)*Qgp(ii) + dt*Pp(ii)];
end
     

%%
date = time + 1950;
lw = 1

figure(1); subplot 221; hold on; grid on 
plot(time/1e3,paleo,'r','linewidth',lw)
xlabel('Time (kyr BP)')
ylabel('T'' (normalized)')
xlim([-31 0.2])

subplot 223; hold on; grid on 
plot(time/1e3,L_out/1e3,'linewidth',lw)
xlabel('kyr BP');
ylabel('L'' (km)')
xlim([-31 0.2])

subplot 222; hold on; grid on 
plot(date,paleo,'r','linewidth',lw)
xlabel('yr CE')
ylabel('T'' (normalized)')
xlim([1000 2100])

subplot 224; hold on; grid on 
plot(date,L_out/1e3,'linewidth',lw)
xlabel('yr CE')
ylabel('L'' (km)')
xlim([1000 2100])