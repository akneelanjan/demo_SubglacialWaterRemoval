%% Formal Plotting Figures 3, 4, 5
clear all; clc;
addpath('lib')
addpath('cbrewer2')

gpuD = gpuDevice;
addpath matplotlib\

icespeed = cbrewer2('seq','RdPu',100);
icespeedmod3 = icespeed(1:90,:);

surfacespeed = cbrewer2('seq','Blues',5);
basalstrength = cbrewer2('seq','Blues',5);
%% Load results

load("ResultsMatFiles\HighRes_Canal_50kPa.mat","xy","t","b","x_surf","x_base"); % Load essential independent variables for plotting

%%% Load essential dependent variables for plotting
Canal50kPa = load("ResultsMatFiles\HighRes_Canal_50kPa.mat","u","u_surf","basal_tau_c_specific");
Canal100kPa = load("ResultsMatFiles\HighRes_Canal_100kPa.mat","u","u_surf","basal_tau_c_specific");
Canal150kPa = load("ResultsMatFiles\HighRes_Canal_150kPa.mat","u","u_surf","basal_tau_c_specific");

RChannel550kPa = load("ResultsMatFiles\HighRes_RChannel_550kPa.mat","u","u_surf","basal_tau_c_specific");
RChannel520kPa = load("ResultsMatFiles\HighRes_RChannel_520kPa.mat","u","u_surf","basal_tau_c_specific");

LinkedCavity250kPa = load("ResultsMatFiles\HighRes_LinkedCavity_250kPa.mat","u","u_surf","basal_tau_c_specific");
LinkedCavity500kPa = load("ResultsMatFiles\HighRes_LinkedCavity_500kPa.mat","u","u_surf","basal_tau_c_specific");

WaterFilm4000Pa = load("ResultsMatFiles\HighRes_WaterFilm_4000Pa.mat","u","u_surf","basal_tau_c_specific");
WaterFilm5333Pa = load("ResultsMatFiles\HighRes_WaterFilm_5333Pa.mat","u","u_surf","basal_tau_c_specific");

%% 2D ice speed map Fig. 3 (a) Canal 50 kPa
figure;
trisurf(t,xy(:,1),xy(:,2),Canal50kPa.u,Canal50kPa.u,'edgecolor','none','facecolor','interp'); hold on
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

u = Canal50kPa.u;
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

%% 2D ice speed map Fig. 3 (b) Canal 100 kPa
figure;
trisurf(t,xy(:,1),xy(:,2),Canal100kPa.u,Canal100kPa.u,'edgecolor','none','facecolor','interp'); hold on
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

u = Canal100kPa.u;
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

%% 2D ice speed map Fig. 3 (c) Canal 150 kPa
figure;
trisurf(t,xy(:,1),xy(:,2),Canal150kPa.u,Canal150kPa.u,'edgecolor','none','facecolor','interp'); hold on
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

u = Canal150kPa.u;
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

%% Ice Surface Speed profiles Fig. 3 (e) Canal
figure;
%x_surf = xy(xy(:,2) > max(xy(:,2))-dx/100,1);
%u_surf = u(xy(:,2) > max(xy(:,2))-dx/100);
plot(x_surf,Canal50kPa.u_surf, ...
    '-pentagram','LineWidth',3,'DisplayName','Q = 0.005 m^3/s','Color',surfacespeed(5,:), ...
    'MarkerIndices', [1 ...
    1+0.5*(length(x_surf)-1)-40 ...
    1+0.5*(length(x_surf)-1) ...
    1+0.5*(length(x_surf)-1)+40 ...
    length(x_surf)]); hold on

plot(x_surf,Canal100kPa.u_surf, ...
    '-o','LineWidth',3,'DisplayName','Q = 0.004 m^3/s','Color',surfacespeed(4,:), ...
    'MarkerIndices', [1 ...
    1+0.5*(length(x_surf)-1)-40 ...
    1+0.5*(length(x_surf)-1) ...
    1+0.5*(length(x_surf)-1)+40 ...
    length(x_surf)]); hold on

plot(x_surf,Canal150kPa.u_surf, ...
    '-square','LineWidth',3,'DisplayName','Q = 0.0025 m^3/s','Color',surfacespeed(3,:), ...
    'MarkerIndices', [1 ...
    1+0.5*(length(x_surf)-1)-40 ...
    1+0.5*(length(x_surf)-1) ...
    1+0.5*(length(x_surf)-1)+40 ...
    length(x_surf)]); hold on

axis([min(xy(:,1)),max(xy(:,1)),0,1.5*max(Canal50kPa.u_surf)])
setFontSize(14)
xlabel("Lateral direction, y [m]", FontSize=16)
ylabel("Ice surface speed [m/yr]", FontSize=16)
legend

%% Basal Drag Profiles Fig. 3 (d) Canal
figure
%yyaxis left
plot(x_base,Canal50kPa.basal_tau_c_specific/1000, ...
    '--pentagram','LineWidth',4,'DisplayName','Q = 0.005 m^3/s','Color',basalstrength(5,:), ...
    'MarkerIndices', [1 1+0.5*(size(x_base,1)-1)-40 1+0.5*(size(x_base,1)-1) 1+0.5*(size(x_base,1)-1)+40 size(x_base,1)]); hold on
plot(x_base,Canal100kPa.basal_tau_c_specific/1000, ...
    '--o','LineWidth',3,'DisplayName','Q = 0.004 m^3/s','Color',basalstrength(4,:), ...
    'MarkerIndices', [1 1+0.5*(size(x_base,1)-1)-40 1+0.5*(size(x_base,1)-1) 1+0.5*(size(x_base,1)-1)+40 size(x_base,1)]); hold on
plot(x_base,Canal150kPa.basal_tau_c_specific/1000, ...
    '-square','LineWidth',2,'DisplayName','Q = 0.0025 m^3/s','Color',basalstrength(3,:), ...
    'MarkerIndices', [1 1+0.5*(size(x_base,1)-1)-40 1+0.5*(size(x_base,1)-1) 1+0.5*(size(x_base,1)-1)+40 size(x_base,1)]); hold on

axis([min(x_base),max(x_base),0,1.1*max(Canal150kPa.basal_tau_c_specific/1000)])
setFontSize(14)
xlabel("Lateral direction, y [m]", FontSize=16)
ylabel("Basal drag [kPa]", FontSize=16)
legend

%% 2D ice speed map Fig. 4 (a) R-Channel 550 kPa
figure;
trisurf(t,xy(:,1),xy(:,2),RChannel550kPa.u,RChannel550kPa.u,'edgecolor','none','facecolor','interp'); hold on
cmap = colormap(icespeedmod3);
set(gca,'ColorScale','linear')
caxis([0 110])
%set(gca,'ColorScale','log')
%colormap(viridis)
clims = [0 110];
%colormap(flipud(cmap));
h = colorbar;
h.Ticks = [25 50 75 100];
h.TickLabels = [25 50 75 100];
h.Label.String = 'u [m/yr]';
h.Label.Position = [1.5 clims(2)+50];
h.Label.Rotation = 0;
h.Label.FontSize = 12;
%view([0 90]); %daspect([2 1 1])

u = RChannel550kPa.u;
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

%% 2D ice speed map Fig. 4 (b) R-Channel 520 kPa
figure;
trisurf(t,xy(:,1),xy(:,2),RChannel520kPa.u,RChannel520kPa.u,'edgecolor','none','facecolor','interp'); hold on
cmap = colormap(icespeedmod3);
set(gca,'ColorScale','linear')
caxis([0 110])
%set(gca,'ColorScale','log')
%colormap(viridis)
clims = [0 110];
%colormap(flipud(cmap));
h = colorbar;
h.Ticks = [25 50 75 100];
h.TickLabels = [25 50 75 100];
h.Label.String = 'u [m/yr]';
h.Label.Position = [1.5 clims(2)+50];
h.Label.Rotation = 0;
h.Label.FontSize = 12;
%view([0 90]); %daspect([2 1 1])

u = RChannel520kPa.u;
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

%% Ice Surface Speed profiles Fig. 4 (d) R-Channel
figure;
%x_surf = xy(xy(:,2) > max(xy(:,2))-dx/100,1);
%u_surf = u(xy(:,2) > max(xy(:,2))-dx/100);
plot(x_surf,RChannel550kPa.u_surf, ...
    '-diamond','LineWidth',3,'DisplayName','Q = 0.086 m^3/s','Color',surfacespeed(5,:), ...
    'MarkerIndices', [1 ...
    1+0.5*(length(x_surf)-1)-40 ...
    1+0.5*(length(x_surf)-1) ...
    1+0.5*(length(x_surf)-1)+40 ...
    length(x_surf)]); hold on
plot(x_surf,RChannel520kPa.u_surf, ...
    '-^','LineWidth',3,'DisplayName','Q = 0.043 m^3/s','Color',surfacespeed(4,:), ...
    'MarkerIndices', [1 ...
    1+0.5*(length(x_surf)-1)-40 ...
    1+0.5*(length(x_surf)-1) ...
    1+0.5*(length(x_surf)-1)+40 ...
    length(x_surf)]); hold on

axis([min(xy(:,1)),max(xy(:,1)),0,1.5*max(RChannel520kPa.u_surf)])
setFontSize(14)
xlabel("Lateral direction, y [m]", FontSize=16)
ylabel("Ice surface speed [m/yr]", FontSize=16)
legend

%% Basal Drag Profiles Fig. 4 (c) R-Channel
figure
%yyaxis left
plot(x_base,RChannel550kPa.basal_tau_c_specific/1000, ...
    '-diamond','LineWidth',3,'DisplayName','Q = 0.086 m^3/s','Color',basalstrength(5,:), ...
    'MarkerIndices', [1 1+0.5*(size(x_base,1)-1)-40 1+0.5*(size(x_base,1)-1) 1+0.5*(size(x_base,1)-1)+40 size(x_base,1)]); hold on
plot(x_base,RChannel520kPa.basal_tau_c_specific/1000, ...
    '-^','LineWidth',3,'DisplayName','Q = 0.043 m^3/s','Color',basalstrength(4,:), ...
    'MarkerIndices', [1 1+0.5*(size(x_base,1)-1)-40 1+0.5*(size(x_base,1)-1) 1+0.5*(size(x_base,1)-1)+40 size(x_base,1)]); hold on

axis([min(x_base),max(x_base),0,1.1*max(RChannel550kPa.basal_tau_c_specific/1000)])
setFontSize(14)
xlabel("Lateral direction, y [m]", FontSize=16)
ylabel("Basal drag [kPa]", FontSize=16)
legend

%% Ice Surface Speed profiles Fig. 5 (b) water film
figure;
%x_surf = xy(xy(:,2) > max(xy(:,2))-dx/100,1);
%u_surf = u(xy(:,2) > max(xy(:,2))-dx/100);
plot(x_surf,WaterFilm4000Pa.u_surf, ...
    '--','LineWidth',4,'DisplayName','Q = Q_0 m^3/s','Color',surfacespeed(5,:), ...
    'MarkerIndices', [1 ...
    1+0.5*(length(x_surf)-1)-40 ...
    1+0.5*(length(x_surf)-1) ...
    1+0.5*(length(x_surf)-1)+40 ...
    length(x_surf)]); hold on
plot(x_surf,WaterFilm5333Pa.u_surf, ...
    '-o','LineWidth',2,'DisplayName','Q = 0.5 Q_0 m^3/s','Color',surfacespeed(4,:), ...
    'MarkerIndices', [1 ...
    1+0.5*(length(x_surf)-1)-40 ...
    1+0.5*(length(x_surf)-1) ...
    1+0.5*(length(x_surf)-1)+40 ...
    length(x_surf)]); hold on

axis([min(xy(:,1)),max(xy(:,1)),0,1.5*max(WaterFilm4000Pa.u_surf)])
setFontSize(14)
xlabel("Lateral direction, y [m]", FontSize=16)
ylabel("Ice surface speed [m/yr]", FontSize=16)
legend

%% Basal Drag Profiles Fig. 5 (a) water film
figure
%yyaxis left
plot(x_base,WaterFilm4000Pa.basal_tau_c_specific/1000, ...
    '--','LineWidth',4,'DisplayName','Q = Q_0 m^3/s','Color',basalstrength(5,:), ...
    'MarkerIndices', [1 1+0.5*(size(x_base,1)-1)-40 1+0.5*(size(x_base,1)-1) 1+0.5*(size(x_base,1)-1)+40 size(x_base,1)]); hold on
plot(x_base,WaterFilm5333Pa.basal_tau_c_specific/1000, ...
    '-o','LineWidth',2,'DisplayName','Q = 0.5 Q_0 m^3/s','Color',basalstrength(4,:), ...
    'MarkerIndices', [1 1+0.5*(size(x_base,1)-1)-40 1+0.5*(size(x_base,1)-1) 1+0.5*(size(x_base,1)-1)+40 size(x_base,1)]); hold on

axis([min(x_base),max(x_base),0,1.1*max(WaterFilm4000Pa.basal_tau_c_specific/1000)])
setFontSize(14)
xlabel("Lateral direction, y [m]", FontSize=16)
ylabel("Basal drag [kPa]", FontSize=16)
legend

%% Ice Surface Speed profiles Fig. 5 (d) linked cavity
figure;
%x_surf = xy(xy(:,2) > max(xy(:,2))-dx/100,1);
%u_surf = u(xy(:,2) > max(xy(:,2))-dx/100);
plot(x_surf,LinkedCavity250kPa.u_surf, ...
    '-square','LineWidth',3,'DisplayName','Q = 0.010 m^3/s','Color',surfacespeed(5,:), ...
    'MarkerIndices', [1 ...
    1+0.5*(length(x_surf)-1)-40 ...
    1+0.5*(length(x_surf)-1) ...
    1+0.5*(length(x_surf)-1)+40 ...
    length(x_surf)]); hold on
plot(x_surf,LinkedCavity500kPa.u_surf, ...
    '-v','LineWidth',3,'DisplayName','Q = 0.005 m^3/s','Color',surfacespeed(4,:), ...
    'MarkerIndices', [1 ...
    1+0.5*(length(x_surf)-1)-40 ...
    1+0.5*(length(x_surf)-1) ...
    1+0.5*(length(x_surf)-1)+40 ...
    length(x_surf)]); hold on

axis([min(xy(:,1)),max(xy(:,1)),0,1.5*max(LinkedCavity250kPa.u_surf)])
setFontSize(14)
xlabel("Lateral direction, y [m]", FontSize=16)
ylabel("Ice surface speed [m/yr]", FontSize=16)
legend

%% Basal Drag Profiles Fig. 5 (c) linked cavity
figure
%yyaxis left
plot(x_base,LinkedCavity250kPa.basal_tau_c_specific/1000, ...
    '--square','LineWidth',3,'DisplayName','Q = 0.010 m^3/s','Color',basalstrength(5,:), ...
    'MarkerIndices', [1 1+0.5*(size(x_base,1)-1)-40 1+0.5*(size(x_base,1)-1) 1+0.5*(size(x_base,1)-1)+40 size(x_base,1)]); hold on
plot(x_base,LinkedCavity500kPa.basal_tau_c_specific/1000, ...
    '-v','LineWidth',3,'DisplayName','Q = 0.005 m^3/s','Color',basalstrength(4,:), ...
    'MarkerIndices', [1 1+0.5*(size(x_base,1)-1)-40 1+0.5*(size(x_base,1)-1) 1+0.5*(size(x_base,1)-1)+40 size(x_base,1)]); hold on

axis([min(x_base),max(x_base),0,1.1*max(LinkedCavity500kPa.basal_tau_c_specific/1000)])
setFontSize(14)
xlabel("Lateral direction, y [m]", FontSize=16)
ylabel("Basal drag [kPa]", FontSize=16)
legend