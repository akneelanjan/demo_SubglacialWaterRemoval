%% Runs Anti-flow line on downsteam (grounding line) cross section. 
% Adjust bed properties for geometry and peak flow speed. Change driving force too.

clear all; clc;
addpath('lib')
addpath('cbrewer2')
%addpath 'D:\SIGMA Group\Ice Streams\Ice Preservation Cooper model\shearmargins_public-master-Neel-Copy\shearmargins_public-master\lib\'
gpuD = gpuDevice;

%% cbrewer2
icespeed = cbrewer2('seq','RdPu',100);
icespeedmod = icespeed(1:75,:);
icespeedmod2 = icespeed(25:75,:);
icespeedmod3 = icespeed(1:90,:);
surfacespeed = cbrewer2('seq','Blues',5);
basalstrength = cbrewer2('seq','Blues',5);
%% Initialization
% Physical parameters
p = 4/3;
rho = 900;
g = 9.8;
%alpha = 2.4E-3;
%alpha = 0.05; % 
%alpha = 0.04; % for R-channel
%alpha = 0.025; % for Canal, Linked Cavity

% bedrock = 50 kPa
%alpha = 0.052; % for Canal
%alpha = 0.051; % for Water film

% bedrock = 10 kPa
alpha = 0.0265; % for Canal
%alpha = 0.0275; % for Water film

f = rho*g*sin(alpha);

factor1 = 1;
%factor2 = 0.5;

% Domain parameters
dx = 10;                     %dx: nominal grid spacing
vert_scale = 1;
dy = vert_scale*dx;

x_max = 2E3*factor1;
y_max = 1000;
%% Modified Domain (idealized mountain glacier setting)
x_bed = [2:-0.05:0]';
y_hardrock2_down = 0.0625*[2:-0.05:1.2]'+0.775;
%y_hardrock2_up = -0.3*1.45+1.30;
y_sediment = 0.85*ones(size([1.15:-0.05:0.85]',1),1);
%y_hardrock1_down = 0.3*0.55+0.70;
y_hardrock1_up = -0.0625*[0.8:-0.05:0]'+0.90;
%y_disc = 0.05*x_disc.^2 -0.1*x_disc + 0.9;

y_bed = [y_hardrock2_down;y_sediment;y_hardrock1_up];

topedge = [0:0.05:2]';
pv5 = zeros(size(x_bed,1)+size(topedge,1)+1,2);
pv5(1:size(topedge,1),1) = 0:0.05:2;
pv5(1:size(topedge,1),2) = 1;
pv5(size(topedge,1)+1:end-1,1) = x_bed;
pv5(size(topedge,1)+1:end-1,2) = y_bed;
pv5(end,:) = [0,1];

figure;
plot(x_bed,y_bed);

pv6 = factor1*1E3.*pv5;
%% Mesh Generation
%pv = 40E3.*[0,1;2,1;2,.93;.3,.9;0,.91;-.75,.92;-.75,1;0,1];
%[xy,t] = distmesh2d(@dpoly,@huniform,dx,40E3.*[-.85,.75;2.1,1.1],pv,pv);

[xy,t] = distmesh2d(@dpoly,@huniform,dx,factor1*1E3.*[-0.1,0.84;2.1,1.1],pv6,pv6);
%xy(:,2) = factor2*((xy(:,2)-y_max) - min(xy(:,2)-y_max));

%% Past strength experiments
% Numerical experiment #4 (rock/sediment with channel)
% beta_var = 5E6;
% beta = 2.9E6;
% rock_trans = 10E3;
% sed_const = 19.4E3;
% sed_var = 8E3;
% sed_trans_l = -23E3;
% sed_trans_r = 71E3;
% ch_str = 24.5E3;
% ch_decay = 5.9E3;
% ch_loc = 10E3;

% Numerical experiment #6 (rock/sediment with channel)
% beta_var = 5E6;
% beta = 2.9E6;
% rock_trans = 0.5E3*factor1;
% sed_const = 19.4E3;
% sed_var = 8E3;
% sed_trans_l = 0.1E3*factor1;
% sed_trans_r = 1.6E3*factor1;
% ch_str = 24.5E3;
% ch_decay = 5E3;
% ch_loc = 0.5E3*factor1;

% tau_c =@(x,y,u) (...
%                  heaviside(sed_trans_l - x).*1E8.*abs(u) + ...
%                  heaviside(x - sed_trans_l).*heaviside(rock_trans - x).*(beta_var*(x/80E3) + beta).^2.*pow_abs(u,5/2) + ...
%                  heaviside(x - rock_trans).*(sed_var*((80E3-x)/80E3) + ...
%                                              ch_str*(exp(-abs(x-ch_loc)/ch_decay)) + ...
%                                              sed_const + ...
%                                              1E8*smearedHeavi(x - sed_trans_r,0)).*abs(u)... 
%                  ).*heaviside((450+dy/2) - y);

% tau_c =@(x,y,u) (...
%                  heaviside(sed_trans_l - x).*1E8.*abs(u) + ...
%                  heaviside(x - sed_trans_l).*heaviside(rock_trans - x).*(beta_var*(x/x_max) + beta).^2.*pow_abs(u,5/2) + ...
%                  heaviside(x - rock_trans).*(sed_var*((x_max-x)/x_max) + ...
%                                              ch_str*(exp(-abs(x-ch_loc)/ch_decay)) + ...
%                                              sed_const + ...
%                                              1E8*heaviside(x - sed_trans_r)).*abs(u)... 
%                  ).*heaviside((500+dy/2) - y);

% R-channel
%tau_c =@(x,y,u) N_ch*exp(-abs(x-ch_loc)/ch_decay).*pow_abs(u,(4/3))/(u_0^(1/3)).*heaviside((900+dy/2) - y);

% Constant N (water film)
%tau_c =@(x,y,u) 4450.*pow_p(u,(4/3))/(u_0^(1/3)).*heaviside((900+dy/2) - y);

% tau_c =@(x,y,u) (...
%                  heaviside(rock_trans - x).*133.*pow_abs(u,5/2) + ...
%                  heaviside(x - rock_trans).*(2E-5 + 1.1*heaviside(x - sed_trans_r)).*abs(u)...
%                 ).*heaviside((900+dy/2) - y);

%% Current hydrology-informed mixed bed strength experiments

%y_base = xy(b(xy(b,2)<905),2);

%u_0 = 1e15; % As = 1e3, C = 0.1, N_wf = 100e3 --> As*(CN_wf)^3 = 1e15
sed_l = 0.8e3;
sed_r = 1.2e3;
N_bedrock = 100e3; % [500 kPa]
Cf = 0.1;
As = 10e3; % [m yr^(-1) MPa^(-3)]
u_0 = As*((0.01)^3); % [m yr^-1]

% R-channels
% N_ch = 550e3; % [Pa]
% ch_decay = 0.5e2;
% ch_loc1 = 1e3;
% ch_loc2 = 1.5e3;
% 
% tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock).*abs(u) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).*((N_ch-N_wf)*exp(-abs(x-ch_loc1)/ch_decay)+N_wf).*abs(u) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock).*abs(u) ...
%                  ).*heaviside((900+dy/2) - y);


% Canal
N_infi = 4e3; % effective pressure far from canal [Pa]
N_canal = 50e3; % canal effective pressure corresponding to flowrates
% Q_w = [0.005, 0.004, 0.002] [m^3/s] {no suction, 20% suction, 60% suction}
a = -0.118/1000; % [m/Pa] from Damsgaard et. al. 2017
b = 4.6; % [m] from Damsgaard et. al. 2017
W_max = a*N_infi + b;
X_canal_center = 1e3;
X_canal_sed_l = 0.8e3+20;
X_canal_sed_r = 1.2e3-20;

% tau_c =@(x,y,u) N_canal*(1-erf(0.02*abs(x-X_canal1)/W_max)).*pow_abs(u,(4/3))/(u_0^(1/3)).*heaviside((900+dy/2) - y) ...
%     + 100*N_canal*(1-erf(0.02*abs(x-X_canal2)/W_max)).*pow_abs(u,(4/3))/(u_0^(1/3)).*heaviside((900+dy/2) - y);

% Velocity-strengthening on hard-bedrock
% tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_wf).*pow_abs(u,(4/3))/(u_0^(1/3)) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).*((N_canal-N_infi)*(1-erf(0.02*abs(x-X_canal)/W_max))+N_infi).*abs(u) ...
%                  + heaviside(x - sed_r).*(Cf*N_wf).*pow_abs(u,(4/3))/(u_0^(1/3)) ...
%                  ).*heaviside((900+dy/2) - y);

% Coulomb-Slip on hard-bedrock
% tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock).*abs(u) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).* ...
%                  ((N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_l)/W_max))+ ...
%                  (N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_center)/W_max))+ ...
%                  (N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_r)/W_max))+ ...
%                  N_infi).*abs(u) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock).*abs(u) ...
%                  ).*heaviside((900+dy/2) - y);

tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock).*abs(u) ...
                 + heaviside(sed_r - x).*heaviside(x - sed_l).* ...
                 ((4*N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_center)/W_max))+ ...
                 N_infi).*abs(u) ...
                 + heaviside(x - sed_r).*(Cf*N_bedrock).*abs(u) ...
                 ).*heaviside((900+dy/2) - y);


% Linked Cavity
% L1 = 1e3; %L2 = 1.5e3;
% a1 = 250e3; %a2 = 250e3;
% sigma = 20; % [m]
% 
% tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock).*abs(u) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).*((a1-N_wf)*exp(-((x-L1).^2)/sigma^2) + N_wf).*abs(u) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock).*abs(u) ...
%                  ).*heaviside((900+dy/2) - y);

% Rock WaterFilmSedimentWaterFilm Rock
N_wf = 4e3; %[in kPa]
% N_wfmod = 10*N_wf/9;
% WFsuction1 = 1e3; %WFsuction2 = 1.5e3;
% delta = 10; % [m]

% tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock).*abs(u) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).*(N_wf).*abs(u) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock).*abs(u) ...
%                  ).*heaviside((900+dy/2) - y);


% tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock).*abs(u) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).* ...
%                  ((N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_l)/W_max))+ ...
%                   N_wf+ ...
%                   (N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_r)/W_max))).*abs(u) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock).*abs(u) ...
%                  ).*heaviside((900+dy/2) - y);


%% Figure 2 Hydrology Schematic with Till Strength profiles
% figure;
% x_base = xy(b(xy(b,2)<905),1);
% y_base = xy(b(xy(b,2)<905),2);
% 
% % Water Film
% % till_tau_c = @(x,y) (heaviside(sed_l - x).*(Cf*N_bedrock) ...
% %                  + heaviside(sed_r - x).*heaviside(x - sed_l).* ...
% %                  ((N_canal-N_infi)*(1-erf(0.2*abs(x-820)/W_max))+ ...
% %                   N_wf+ ...
% %                   (N_canal-N_infi)*(1-erf(0.2*abs(x-1180)/W_max))) ...
% %                  + heaviside(x - sed_r).*(Cf*N_bedrock) ...
% %                  ).*heaviside((900+dy/2) - y);
% 
% % Linked Cavity
% %till_tau_c = @(x,y) ((a1-N_wf)*exp(-((x-L1).^2)/sigma^2) + N_wf).*heaviside((900+dy/2) - y);
% 
% % Canal
% %till_tau_c = @(x,y) ((N_canal-N_infi)*(1-erf(0.02*abs(x-X_canal)/W_max))+N_infi).*heaviside((900+dy/2) - y);
% till_tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).* ...
%                  ((N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_l)/W_max))+ ...
%                  (N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_center)/W_max))+ ...
%                  (N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_r)/W_max))+ ...
%                  N_infi) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock) ...
%                  ).*heaviside((900+dy/2) - y);
% 
% % R-Channel
% %till_tau_c = @(x,y) ((N_ch-N_wf)*exp(-abs(x-ch_loc1)/ch_decay)+N_wf).*heaviside((900+dy/2) - y);
% 
% till_tau_c_specific = till_tau_c(x_base,y_base);
% 
% plot(x_base,till_tau_c_specific/1000,'LineWidth',3,'DisplayName','Q = 0.005 m^3/s','Color',basalstrength(5,:)); hold on
% 
% axis([0,x_max,0,1.1*max(till_tau_c_specific/1000)])
% %axis([850,1150,0,10])
% setFontSize(14)
% xlabel("Lateral direction, y [m]", FontSize=16)
% ylabel("Basal yield strength [kPa]", FontSize=16)
%%
% baseArray_tau_c =@(x) (heaviside(sed_l - x).*(Cf*N_bedrock) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).* ...
%                  ((N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_center)/W_max))+ ...
%                  N_infi) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock) ...
%                  );
% x_baseArray = 0:10:x_max;
% baseArray_tau_c_specific = baseArray_tau_c(x_baseArray);
% 
% figure;
% plot(x_baseArray,baseArray_tau_c_specific/1000,'LineWidth',3,'DisplayName','Q = 0.005 m^3/s','Color',basalstrength(5,:)); hold on
% 
% axis([0,x_max,0,1.1*max(baseArray_tau_c_specific/1000)])
% %axis([850,1150,0,10])
% setFontSize(14)
% xlabel("Lateral direction, y [m]", FontSize=16)
% ylabel("Basal yield strength [kPa]", FontSize=16)

%%
nN = size(xy,1);                     %nN: number of nodes
nE = size(t,1);                      %nE: number of elements
b = unique(boundedges(xy,t));        %b:  boundary node numbers
nB = size(b,1);                      %nB: number of boundary nodes

%%
e = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
e = sort(e,2);
[foo,~,ifoo] = unique(e,'rows');
eB = foo(accumarray(ifoo,1) == 1,:); %eB: boundary edges

T = 25/y_max*(y_max-xy(:,2)) + 248;
T_old = 25/y_max*(y_max-xy(:,2)) + 248;
res = 1;
while(res > 1E-3)
T = .3*T + (1-.3)*T_old;

A = 3.5e-25; %A: preexponential constant[1/s*Pa^3]
E = 1; %E: enhancement factor[]
Q_h = 0*xy(:,1); %Q_h: melting point adjusted activation energy[J/mol]
R = 8.314; %R: gas constant[J/K*mol]
Tstar = 263.15; %Tstar: activation threshold[K] (-10C)
p0 = 7e-8; %p0: pressure heating coef[K/Pa]
P = rho*g*(y_max-xy(:,2)); %P: hydrostatic pressure[Pa]
T_h = T + p0*P; %T_h: melting point adjusted temperature[K]
Tstar_h = Tstar*(0*xy(:,1) + 1) + p0*P; %Tstar_h: adjusted activation threshold[K]
Q2 = 60e3; %Q2: lower activation energy[J/mol]
Q3 = 115e3; %Q3: higher activation energy[J/mol]
Q_h(heaviside(T_h-Tstar_h)==0) = Q2;
Q_h(heaviside(T_h-Tstar_h)==1) = Q3;
a = 1/2*(A*E*exp((-Q_h./R).*((1./T_h)-(1./Tstar_h)))).^(-1/3); %a = A^(-1/n)

%% Build System
A = sparse(nE,nN);
B = sparse(nE,nN);
D = zeros(nE,nN);
F = zeros(1,nN);

tau_area = zeros(nE,1);
for E = 1:nE  % integration over each element
  nodes = t(E,:);
  xyE = [ones(3,1),xy(nodes,:)];
  Area = abs(det(xyE))/2;
  tau_area(E) = Area;
  C = inv(xyE);
  A_E = C(2,:);
  B_E = C(3,:);
  D_E = 1/3*ones(1,3);
  F_E = f*Area/3*ones(1,3);
  A(E,nodes) = A(E,nodes) + A_E;
  B(E,nodes) = B(E,nodes) + B_E;
  D(E,nodes) = D(E,nodes) + D_E;
  F(nodes) = F(nodes) + F_E;
end

b_dx = zeros(size(b,1),1);
for N = 1:size(b,1) % integration over each boundary node
    edges = eB(sum(eB == b(N),2) == 1,:);
    for i = 1:size(edges,1)
        edge = edges(i,:);
        b_dx(N) = b_dx(N) + .5*sqrt(sum(diff(xy(edge,:),1).^2,2));
    end
end

%% Solve

cvx_begin
cvx_quiet true
    variables u(nN)
    obj = (1/p)*sum((D*a).*tau_area.*pow_pos(norms([A*u,B*u],2,2),p)) - F*u + sum(b_dx.*tau_c(xy(b,1),xy(b,2),u(b)));
    subject to
        u >= 0;
        u(xy(:,1) > x_max - dx/2) == 0;
        u(xy(:,1) < 0 + dx/2) == 0;
    minimize(obj)
cvx_end

u_x = A*u;
u_y = B*u;
mu = (D*a).*(sqrt((u_x).^2 + (u_y).^2)).^(1/2);
tau_E = sqrt((mu.*u_x).^2 + (mu.*u_y).^2); %tau_E: effective stress[Pa]
epsilon_E = sqrt((u_x/2).^2 + (u_y/2).^2); %epsilon_E: effective strain rate[1/s]
f_therm = 2*tau_E.*epsilon_E;

k1 = 9.828; %k1: conductivity preexponential[W/m*K]
k2 = 5.7; %k2: conductivity postexponential[1/K]
k = k1*exp(-k2*1e-3.*T);
tau_T =@(y) (55./(max(xy(:,2))-y.*heaviside((900+dy/2) - y))).*heaviside((900+dy/2) - y);

F_therm = zeros(1,nN);
for E = 1:nE  % integration over each element
  nodes = t(E,:);
  xyE = [ones(3,1),xy(nodes,:)];
  Area = abs(det(xyE))/2;
  F_therm_E = f_therm(E)*Area/3*ones(1,3);
  F_therm(nodes) = F_therm(nodes) + F_therm_E;
end

T_old = T;
cvx_begin
cvx_quiet true
    variables T(nN)
    obj = sum((D*k).*tau_area.*pow_pos(norms([A*T,B*T],2,2),2)) - 2E10*F_therm*T/40 - sum(b_dx.*k(b).*tau_T(xy(b,2)).*T(b));
    subject to
        T <= 273;
        T(xy(:,2) > max(xy(:,2))-dx) == 248;
    minimize(obj)
cvx_end

res = norm(T-T_old)/norm(T);
disp(res)

end
%%
u = u*pi*1E7; %conversion from [m/s] to [m/yr]
T = T - 273; %conversion from K to deg. C

%% Visualization
x_base = xy(b(xy(b,2)<905),1);
y_base = xy(b(xy(b,2)<905),2);
u_base = u(b(xy(b,2)<905));

%% 2D speed map
figure;
trisurf(t,xy(:,1),xy(:,2),u,u,'edgecolor','none','facecolor','interp'); hold on
cmap = colormap(icespeedmod3);
set(gca,'ColorScale','linear')
caxis([0 100])
%set(gca,'ColorScale','log')
%colormap(viridis)
clims = [0 100];
%colormap(flipud(cmap));
h = colorbar;
h.Ticks = [25 50 75 100];
h.TickLabels = [25 50 75 100];
h.Label.String = 'u [m/yr]';
h.Label.Position = [1.5 clims(2)+50];
h.Label.Rotation = 0;
h.Label.FontSize = 12;
%view([0 90]); %daspect([2 1 1])

% left edge no-slip marking
plot3(xy(b(u(b) < 1e-6 & xy(b,1)==0),1),xy(b(u(b) < 1e-6 & xy(b,1)==0),2),max(u)+0*b(u(b) < 1e-6 & xy(b,1)==0), ...
    'LineWidth',3,'Color',[0.4 0.4 0.4])

% right edge no-slip marking
plot3(xy(b(u(b) < 1e-6 & xy(b,1)==2000),1),xy(b(u(b) < 1e-6 & xy(b,1)==2000),2),max(u)+0*b(u(b) < 1e-6 & xy(b,1)==2000), ...
    'LineWidth',3,'Color',[0.4 0.4 0.4])

% bottom edge induced no-slip marking
plot3(xy(b(u(b) < 1e-6 & xy(b,2)==850),1),xy(b(u(b) < 1e-6 & xy(b,2)==850),2),max(u)+0*b(u(b) < 1e-6 & xy(b,2)==850), ...
    'LineWidth',3,'Color',[0.6 0.6 0.6])

view(2)
setFontSize(14)
xlabel("Lateral direction, y [m]", FontSize=16)
ylabel("Vertical direction, z [m]", FontSize=16)
%title("Antiplane ice flow speed map")
axis equal tight

%% Surface Ice Speed profile
figure;
x_surf = xy(xy(:,2) > max(xy(:,2))-dx/100,1);
u_surf = u(xy(:,2) > max(xy(:,2))-dx/100);
plot(x_surf,u_surf, ...
    '-pentagram','LineWidth',3,'DisplayName','Q = 0.002 m^3/s','Color',surfacespeed(2,:), ...
    'MarkerIndices', [1 ...
    1+0.5*(length(x_surf)-1)-40 ...
    1+0.5*(length(x_surf)-1) ...
    1+0.5*(length(x_surf)-1)+40 ...
    length(x_surf)]); hold on
axis([min(xy(:,1)),max(xy(:,1)),0,1.5*max(u)])
setFontSize(14)
xlabel("Lateral direction, y [m]", FontSize=16)
ylabel("Ice surface speed [m/yr]", FontSize=16)

%% 2D Temperature Contour
figure('Position',[10 10 1000 200])
load iceColorMap
colormap(iceColorMap)
% tricontf(xy(:,1),xy(:,2),t,T,...
%          [0,-5,-10,-15,-20,-25]);
trisurf(t,xy(:,1),xy(:,2),T,T,...
       'edgecolor','none','facecolor','interp');
axis off
view(2)
xlabel("Lateral direction [m]", FontSize=16)
ylabel("Vertical direction [m]", FontSize=16)


%% Basal Strength Profile formulations without speed-multiplication
% Coulomb-slip


% R-channel case
% basal_tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).*((N_ch-N_wf)*exp(-abs(x-ch_loc1)/ch_decay)+N_wf) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock) ...
%                  ).*heaviside((900+dy/2) - y);

% Linked Cavity case
% basal_tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).*((a1-N_wf)*exp(-((x-L1).^2)/sigma^2) + N_wf) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock) ...
%                  ).*heaviside((900+dy/2) - y);

% Water Film Case
% basal_tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).*(N_wf) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock) ...
%                  ).*heaviside((900+dy/2) - y);

% basal_tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).* ...
%                  ((N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_l)/W_max))+ ...
%                   N_wf+ ...
%                   (N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_r)/W_max))) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock) ...
%                  ).*heaviside((900+dy/2) - y);

% basal_tau_c =@(x,y,u) (heaviside(sed_l - x).*0.2.*rho.*g.*(y_max-y).*pow_p(u,(1/3))/(u_0^(1/3)) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).*((N_wfmod-N_wf)*exp(-abs(x-WFsuction1)/delta) + 2*N_wf) ...
%                  + heaviside(x - sed_r).*0.2.*rho.*g.*(y_max-y).*pow_p(u,(1/3))/(u_0^(1/3)) ...
%                  ).*heaviside((900+dy/2) - y);

% Canal case

% Velocity-strengthening
% basal_tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_wf).*pow_p(u,(1/3))/(u_0^(1/3)) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).*((N_canal-N_infi)*(1-erf(0.02*abs(x-X_canal_center)/W_max))+N_infi) ...
%                  + heaviside(x - sed_r).*(Cf*N_wf).*pow_p(u,(1/3))/(u_0^(1/3)) ...
%                  ).*heaviside((900+dy/2) - y);

% Coulomb-slip
% basal_tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).* ...
%                  ((N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_l)/W_max))+ ...
%                  (N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_center)/W_max))+ ...
%                  (N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_sed_r)/W_max))+ ...
%                  N_infi) ...
%                  + heaviside(x - sed_r).*(Cf*N_bedrock) ...
%                  ).*heaviside((900+dy/2) - y);

basal_tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_bedrock) ...
                 + heaviside(sed_r - x).*heaviside(x - sed_l).* ...
                 ((4*N_canal-N_infi)*(1-erf(0.2*abs(x-X_canal_center)/W_max))+ ...
                 N_infi) ...
                 + heaviside(x - sed_r).*(Cf*N_bedrock) ...
                 ).*heaviside((900+dy/2) - y);

% basal_tau_c =@(x,y,u) (heaviside(sed_l - x).*(Cf*N_wf) ...
%                  + heaviside(sed_r - x).*heaviside(x - sed_l).*((N_canal-N_infi)*(1-erf(0.02*abs(x-X_canal)/W_max))+N_infi) ...
%                  + heaviside(x - sed_r).*(Cf*N_wf) ...
%                  ).*heaviside((900+dy/2) - y);

basal_tau_c_specific = basal_tau_c(x_base,y_base,u_base);
%%
figure
%yyaxis left
plot(x_base,basal_tau_c_specific/1000, ...
    '-pentagram','LineWidth',3,'DisplayName','Q = 0.002 m^3/s','Color',basalstrength(2,:), ...
    'MarkerIndices', [1 1+0.5*(size(x_base,1)-1)-40 1+0.5*(size(x_base,1)-1) 1+0.5*(size(x_base,1)-1)+40 size(x_base,1)]); hold on

axis([min(x_base),max(x_base),0,1.1*max(basal_tau_c_specific/1000)])
setFontSize(14)
xlabel("Lateral direction, y [m]", FontSize=16)
ylabel("Basal drag [kPa]", FontSize=16)
%% Make sound when done
beep
pause(1)
beep
pause(1)
beep
%% Save to .mat file
filepath = "D:\SIGMA Group\Ice Streams\The Cryosphere Basal Anchoring Manuscript Submission\Mat files\";
filename = filepath+"HighRes_only1Canal_"+"200"+"kPa.mat";
save(filename)
