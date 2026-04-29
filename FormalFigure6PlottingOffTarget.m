%% Formal Plotting Figure 6
clear all; clc;
addpath('lib')
addpath('cbrewer2')

gpuD = gpuDevice;
addpath matplotlib\
%% Load x-axis
load("ResultsMatFiles\OffTarget_HighRes_Canal50kPaY1000.mat","x_base");

%% Load results for baseline, 20% on-target flowrate reduction, and 4 cases of off-target flux reduction

%%% AUC of Ice Surface Speed Profiles
load("ResultsMatFiles\OffTarget_HighRes_Canal50kPaY1000.mat","IceFlux50kPaY1000");

load("ResultsMatFiles\OffTarget_HighRes_Canal100kPaY1000.mat","IceFlux100kPaY1000");

load("ResultsMatFiles\OffTarget_Canal50kPa_Blip1333Pa_XSuction1010.mat","IceFluxXSuction1010");

load("ResultsMatFiles\OffTarget_Canal50kPa_Blip1333Pa_XSuction1020.mat","IceFluxXSuction1020");

load("ResultsMatFiles\OffTarget_Canal50kPa_Blip1333Pa_XSuction1050.mat","IceFluxXSuction1050");

load("ResultsMatFiles\OffTarget_Canal50kPa_Blip1333Pa_XSuction1100.mat","IceFluxXSuction1100");

%%% Basal Drag Profiles
load("ResultsMatFiles\OffTarget_HighRes_Canal50kPaY1000.mat",'basaltauc_Canal50kPaY1000');

load("ResultsMatFiles\OffTarget_HighRes_Canal100kPaY1000.mat",'basaltauc_Canal100kPaY1000');

load("ResultsMatFiles\OffTarget_Canal50kPa_Blip1333Pa_XSuction1010.mat",'basaltauc_Canal50kPaY1010');

load("ResultsMatFiles\OffTarget_Canal50kPa_Blip1333Pa_XSuction1020.mat",'basaltauc_Canal50kPaY1020');

load("ResultsMatFiles\OffTarget_Canal50kPa_Blip1333Pa_XSuction1050.mat",'basaltauc_Canal50kPaY1050');

load("ResultsMatFiles\OffTarget_Canal50kPa_Blip1333Pa_XSuction1100.mat",'basaltauc_Canal50kPaY1100');


%% Compute Ice Speed Reduction Efficiencies

y_labels = [1000, 1010, 1020, 1050, 1100];
reduction = [(IceFlux50kPaY1000 - IceFlux100kPaY1000), (IceFlux50kPaY1000 - IceFluxXSuction1010), (IceFlux50kPaY1000 - IceFluxXSuction1020), (IceFlux50kPaY1000 - IceFluxXSuction1050), (IceFlux50kPaY1000 - IceFluxXSuction1100)];

blue = cbrewer2('seq','Blues',5);
reds = cbrewer2('seq','Reds',6);

color1 = blue(4,:); color2 = reds(6,:); color3 = reds(5,:); color4 = reds(4,:); color5 = reds(3,:);

reductionEfficiency = reduction/max(reduction)*100;


%% Figure 6 Panel (b)
figure;
plot(y_labels,reductionEfficiency,'-k','MarkerSize',18,'HandleVisibility','off'); hold on
scatter(y_labels(1),reductionEfficiency(1),24,color1,"filled",'DisplayName','y_{off} = 0 m'); hold on
scatter(y_labels(2),reductionEfficiency(2),36,color2,"filled","square",'DisplayName','y_{off} = 10 m'); hold on
scatter(y_labels(3),reductionEfficiency(3),36,color3,"filled","^",'DisplayName','y_{off} = 20 m'); hold on
scatter(y_labels(4),reductionEfficiency(4),36,color4,"filled","o",'DisplayName','y_{off} = 50 m'); hold on
scatter(y_labels(5),reductionEfficiency(5),36,color5,"filled","diamond",'DisplayName','y_{off} = 100 m'); hold on
setFontSize(14)
xlabel("Lateral direction [m], y", FontSize=16)
ylabel("Ice speed reduction, %", FontSize=16)
axis tight
axis([1000,1100,0,100])
legend

%% Figure 6 Panel (a)
figure
%yyaxis left
plot(x_base,basaltauc_Canal50kPaY1000/1000,'-pentagram','LineWidth',3,'DisplayName','Q = 0.005 m^3/s','Color',blue(5,:)); hold on
plot(x_base,basaltauc_Canal100kPaY1000/1000,'--pentagram','LineWidth',1,'DisplayName','Q = 0.004 m^3/s','Color',blue(4,:)); hold on
plot(x_base,basaltauc_Canal50kPaY1010/1000,'--+','LineWidth',1,'DisplayName','y_{off} = 10 m','Color',reds(6,:)); hold on
plot(x_base,basaltauc_Canal50kPaY1020/1000,'-^','LineWidth',1,'DisplayName','y_{off} = 20 m','Color',reds(5,:)); hold on
plot(x_base,basaltauc_Canal50kPaY1050/1000,'-o','LineWidth',1,'DisplayName','y_{off} = 50 m','Color',reds(4,:)); hold on
plot(x_base,basaltauc_Canal50kPaY1100/1000,'-diamond','LineWidth',1,'DisplayName','y_{off} = 100 m','Color',reds(3,:)); hold on

axis([min(x_base),max(x_base),0,1.1*max(basaltauc_Canal100kPaY1000/1000)])
setFontSize(14)
xlabel("Lateral direction [m]", FontSize=16)
ylabel("Basal drag profile [kPa]", FontSize=16)
legend