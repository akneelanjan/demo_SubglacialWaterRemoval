%% Parameters for 2 side canals and On-Target flowrate reduction in central Water Film
clear all; clc
%%
addpath('..\lib')
addpath('..\cbrewer2')
gpuD = gpuDevice;

% Cbrewer2 colorscales, perceptually uniform, color-blind friendly
surfacespeed = cbrewer2('seq','Blues',5);
basalstrength = cbrewer2('seq','Blues',5);
%%
dx = 10; %dx: nominal grid spacing along lateral direction
vert_scale = 1; % ratio of grid-spacing in vertical direction relative to the lateral direction 
dy = vert_scale*dx;
%%
alpha = 0.027; % for Water Film
Q = 'Q_0'; % Water film flux [m^3/s]

sed_l = 0.8e3; % Rock-sediment interface on the left
sed_r = 1.2e3; % Rock-sediment interface on the right
N_bedrock = 100e3; % [100 kPa] % Film Effective pressure at the Ice-Bedrock contact
Cf = 0.1; % Sliding Coefficient
As = 10e3; % [m yr^(-1) MPa^(-3)] % Hard-bed drag parameter
u_0 = As*((0.01)^3); % [m yr^-1] % Threshold to verify Coulomb-Slip regime

N_wf = 4e3; %[in kPa] % Sediment Cohesion

% Canal parameters
N_infi = 4e3; % far-field effective pressure away from canal [Pa]
N_canal = 50e3; % Canal effective pressure 

% N_canal = {50,100,150,200} [kPa] corresponding to flowrates:
% Q = {0.005,0.004,0.0025,0.002} [m^3/s] (baseline, 20% flowrate reduction, 50% flowrate reduction, 60% flowrate reduction)

a_canal = -0.118/1000; % [m/Pa] from Damsgaard et. al. 2017
b_canal = 4.6; % [m] from Damsgaard et. al. 2017
W_max = a_canal*N_infi + b_canal; % Maximum width of canal, [m], from Damsgaard et. al. 2017

X_canal_sed_l = 0.8e3+20; % Location of side canal on the left
X_canal_sed_r = 1.2e3-20; % Location of side canal on the right

% Water-film parameters
N_wfmod = N_wf; % Strength perturbation in Water Film [Pa]
% N_wfmod = N_wf = 4 kPa for Baseline
% N_wfmod = 10*N_wf/9 = 4444 Pa for 30% flowrate reduction
% N_wfmod = 4*N_wf/3 = 5333 Pa for 30% flowrate reduction
% N_wfmod = 2*N_wf = 8 kPa for 30% flowrate reduction

WFsuction1 = 1e3; % Location of flux reduction in the water film (center of the bed)
delta = 10; % [m] Spatial scale of the Gaussian distribution for the local rise in basal strength

% Basal-drag, tau_bed formulation as per Coulomb-Slip on hard-bedrock
tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock).*abs(u) ...
                 + heaviside(sed_r - x).*heaviside(x - sed_l).* ...
                 ((N_canal-N_wf)*(1-erf(0.2*abs(x-X_canal_sed_l)/W_max))+ ...
                  ((N_wfmod-N_wf)*exp(-abs(x-WFsuction1)/delta)+N_wf)+ ...
                  (N_canal-N_wf)*(1-erf(0.2*abs(x-X_canal_sed_r)/W_max))).*abs(u) ...
                 + heaviside(x - sed_r).*(Cf*N_bedrock).*abs(u) ...
                 ).*heaviside((900+dy/2) - y);
% with speed multiplied

basal_tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock) ...
                 + heaviside(sed_r - x).*heaviside(x - sed_l).* ...
                 ((N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_l)/W_max))+ ...
                  ((N_wfmod-N_wf)*exp(-abs(x-WFsuction1)/delta)+N_wf)+ ...
                  (N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_r)/W_max))) ...
                 + heaviside(x - sed_r).*(Cf*N_bedrock) ...
                 ).*heaviside((900+dy/2) - y);
% original (without speed multiplied)

% Ice surface speed plot legend color
surfaceSpeedColor = surfacespeed(5,:);

% Basal strength plot legend color
basalStrengthColor = basalstrength(5,:);

%% save parameters in a file
filename = "..\InputParameterFiles\WaterFilm_4000Pa_Baseline.mat";
save(filename)