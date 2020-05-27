%%% script to generate (highly) idealized paleoclimate timeseries for
%%% investigating timescales of outlet glacier terminus changes. (see
%%% Christian et al., 2020, TC, fig 7). 
%%% JEC May 2020
clear all;

lgm_mag = -10; % relative magnitude of cartoon deglaciation
lgm_t0 = 11e3; % years BP
lgm_tau = 5e2; % width of error function used to smooth (yrs)

lia_mag = -1; % relative mag. of cartoon LIA
lia_t0 = 500; % LIA start; years BP
lia_t1 = 100; % LIA end;
lia_tau = 5e1; % width of error function used to smooth (yrs)

anth_mag = 4; % anthropogenic climate change magnitude
anth_t0 = 70; % anthropogenic onset (yrs BP)
anth_t1 = -150; % end of trend

present = 1950; % year CE
start = 3e4; % years BP
fin = anth_t1; % years in future (from 1950 or whenever present is set)

time = [-start:-fin];
date = time+present;

% individual anomalies:
lgm = lgm_mag*0.5*(1-erf((time+lgm_t0)/lgm_tau));
lia = lia_mag*0.5*(erf((time+lia_t0)/lia_tau)) + lia_mag*0.5*(-erf((time+lia_t1)/lia_tau));
anth = [zeros(1,start-anth_t0),linspace(0,anth_mag,anth_t0-anth_t1+1)];
% combined:
lgm_lia_anth = lgm+lia+anth;

% figure(1)
% plot(time/1e3,lgm_lia_anth,'linewidth',2); grid on
% xlabel('kyr BP');
% ylabel('Anomaly (arbitrary units)')