clear all
close all
cacheDir = './cache'; %fullfile(strrep(userpath,pathsep,''), 'cache');
forceRecalc = false; 
% define units and constants objects
u = units;
const = constants;

%% Parallel Computing
% matlabpool open 3
% matlabpool open 8


%% Set simulation parameters
time            = (-2:0.1:1000)*u.ps;           % the time we want to simulate
dbf             = 0;                        % Debye-Waller factor
initTemp        = 300*u.K;                  % initial temperature of the sample
fluence         = 10*u.mJ/u.cm^2;           % fluence we pump the sample
pulseWidth      = 0.1*u.ps;
heatDiffusion   = true;                     % enable heat diffusion

%% Initialize the Sample
disp('Initialize sample ...');

%% Initialize the atoms of the Sample Structure
O = atomBase('O');
Mg = atomBase('Mg');
Fe = atomBase('Fe');
Pt = atomBase('Pt');
Al = atomBase('Al');

%% Initialize the unitCells of the Sample

% Sample structure
% 2nm Al / 10nm FePt / MgO
% 8nm FePt /MgO

% lattice constants in Angstrom
cFePt       = 3.71      *u.ang;
aFePt       = 3.85      *u.ang;
cMgOsub     = 4.212     *u.ang;
aAl         = 4.0494    *u.ang;
%Thin film thicknesses
%dFePt = 10*u.nm ;
dFePt = 10*u.nm ;
nFePt = round(dFePt/cFePt) ;

dAl = 2*u.nm;
nAl = round(dAl/aAl);


% sound velocities [nm/ps]
svMgO       = 9.120     *u.nm/u.ps; % taken from Reichmann et al, Geophys. Res. Lett. 27, 799 (2000).
svFePt      = 3.64      *u.nm/u.ps; % 2.2nm according to 
svAl        = 5.0       *u.nm/u.ps; % according to Wikipedia 

MAl = 26.98;
% %%% Al layer 
% 
propAl.aAxis           = aAl;                  % aAxis
propAl.bAxis           = aAl;                  % bAxis
propAl.debWalFac       = dbf;                    % Debye-Waller factor
propAl.soundVel        = svAl;                 % sound velocity
propAl.optPenDepth     = 15*u.nm;                % value assumed by Kimling and Kimling
propAl.thermCond       = 237*u.W/(u.m *u.K);     % heat conductivity Wikipedia
propAl.linThermExp     = 2.31e-5;                % linear thermal expansion Wikipedia
propAl.heatCapacity    = 897*u.J/(u.kg *u.K) ; %Wikipedia 

AlUC = unitCell('AlUC','AlUC',aAl,propAl);
AlUC.addAtom(Al,0);AlUC.addAtom(Al,0)
AlUC.addAtom(Al,0.5);AlUC.addAtom(Al,0.5);
AlUC.optRefIndex = [2.7673, 8.3543];


% %%% FePt layer (001)
propFePt.aAxis           = aFePt;                  % aAxis
propFePt.bAxis           = aFePt;                  % bAxis
propFePt.debWalFac       = dbf;                    % Debye-Waller factor
propFePt.soundVel        = svFePt;                 % sound velocity
propFePt.optPenDepth     = 15*u.nm;                % value assumed by Kimling and Kimling
propFePt.thermCond       = 5.72*u.W/(u.m *u.K);    % heat conductivity CHECK
propFePt.linThermExp     = 3.4e-5;                % linear thermal expansion CHECK
% propFePt.heatCapacity    = @(T)((107.7 + 2.65e-2.*T -5.19e5./T.^2)/(6.0200e+023*3.9303e-025));% heat capacity [J / kg K]
propFePt.heatCapacity    = (3.46+0.033+0.17)*1e6 *u.J/((u.m)^3 *u.K)/ 15154*u.kg/(u.m)^3; % Kimling and Kimling et al. 

FePt = unitCell('FePt','FePt',cFePt,propFePt);
FePt.addAtom(Pt,0);FePt.addAtom(Pt,0)
FePt.addAtom(Fe,0.5);FePt.addAtom(Fe,0.5);
FePt.optRefIndex = [3.33, 2.635];

% FePt.mass = FePt.mass / FePt.aAxis/ FePt.bAxis * (1*u.ang)^2;
% FePt.aAxis = 1*u.ang;
% FePt.bAxis = 1*u.ang;
% FePt.area =  FePt.aAxis * FePt.bAxis;
% FePt.volume = FePt.cAxis * FePt.area;
% FePt.soundVel = svFePt; % one need to redefine sound velocity here since it adapts the spring constant to the modified unit cell mass
% 

%%% MgO substrate

propMgOsub.aAxis           = cMgOsub;                  % aAxis
propMgOsub.bAxis           = cMgOsub;                  % bAxis
propMgOsub.debWalFac       = dbf;                      % Debye-Waller factor
propMgOsub.soundVel        = svMgO;                    % sound velocity
propMgOsub.optPenDepth     = Inf;                      % optical penetration depth
propMgOsub.thermCond       = 45*u.W/(u.m *u.K);        % heat conductivity
propMgOsub.linThermExp     = 1.4e-5;                   % linear thermal expansion
propMgOsub.heatCapacity    = 37.8/0.0403044;           % heat capacity [J / kg K]


MgOsub = unitCell('MgOsub', 'MgOsub', cMgOsub, propMgOsub);
MgOsub.addAtom(Mg,0);
MgOsub.addAtom(O,0);
MgOsub.addAtom(Mg,0);
MgOsub.addAtom(O,0);
MgOsub.addAtom(Mg,0.5);
MgOsub.addAtom(O,0.5);
MgOsub.addAtom(Mg,0.5);
MgOsub.addAtom(O,0.5);

MgOsub.optRefIndex = [1.7276 ,    0.0];

% MgOsub.mass = MgOsub.mass / MgOsub.aAxis/ MgOsub.bAxis * (1*u.ang)^2;
% MgOsub.aAxis = 1*u.ang;
% MgOsub.bAxis = 1*u.ang;
% MgOsub.area =  MgOsub.aAxis * MgOsub.bAxis;
% MgOsub.volume = MgOsub.cAxis * MgOsub.area;
% MgOsub.soundVel = svMgO; % one need to redefine sound velocity here since it adapts the spring constant to the modified unit cell mass
% 

%% Initialize the Sample Structure

S = structure('FePt-thinFilm');
S.addSubStructure(AlUC,nAl);
S.addSubStructure(FePt,nFePt);
S.addSubStructure(MgOsub,400);  % sound velocity in DSO ~ 20 u.c./ps


% add a static substrate to the sample structure
substrate = structure('MgOsubstrate');
substrate.addSubStructure(MgOsub,1000000)% 
S.addSubstrate(substrate);
S.display()
%S.visualize()

distances = S.getDistancesOfUnitCells(); % these are the distances of each unitCell from the surface





%% Calculating the excitation profile
angle        = 44;
wavelength   = 800e-9;
stepsOptical = 250;

interfaces = S.getDistancesOfInterfaces;
dAir  = 3*u.nm; 
dAl   = interfaces(2) - interfaces(1);
dFePt = interfaces(3) - interfaces(2);
dMgO  = interfaces(4) - interfaces(3);

layers800 =[1.000  ,    0.0,       1e-6,     1.0,    dAir;                                                         ... %common sense values for air
            AlUC.optRefIndex(1) ,      AlUC.optRefIndex(2),  AlUC.density,    propAl.heatCapacity,          dAl;     ... %refractive Index.info
            FePt.optRefIndex(1)   ,  FePt.optRefIndex(2),    FePt.density,    propFePt.heatCapacity,      dFePt;     ... %Z. H. Cen, B. X. Xu, J. F. Hu, J. M. Li, K. M. Cher, Y. T. Toh, K. D. Ye, and J. Zhang, "Optical property study of FePt-C nanocomposite thin film for heat-assisted magnetic recording," Opt. Express 21, 9906-9914 (2013) 
            MgOsub.optRefIndex(1) ,  MgOsub.optRefIndex(2),  MgOsub.density,  propMgOsub.heatCapacity,    dMgO];         %refractive Index.info
    


[zs8,Ints8,dAs8,dTs8,R8,T8] = multilayerAbsorption(layers800, angle, wavelength,stepsOptical);  
[z8,Int8,dA8,dT8]           = combineLayers(zs8, Ints8, dAs8, dTs8);

% 


% Heat diffusion
% Initialize the Heat Diffusion Simulation
%Provide the sample structure for the heat simulation.
H = heat(S,forceRecalc,heatDiffusion);
H.setCacheDir(cacheDir); % set the cache directory

z8cut = dA8(z8-dAir>0)-dAir;
dA8cut = dA8(z8-dAir>=0);
H.dAlphadz = interp1(z8-dAir,dA8,S.getDistancesOfUnitCells());
H.dAlphadz(1) = dA8cut(2);
fontSize = 16;
figure(2);
plot(z8/u.nm-dAir/u.nm, dA8*1e-9*100, 'LineWidth',2, 'Color','red')
hold on;
plot(distances/u.nm,H.getAbsorptionProfile*1e-9*100,'-o','LineWidth',2, 'Color','blue')
plot(distances/u.nm,H.dAlphadz*1e-9*100,'-s','LineWidth',2, 'Color','green')

legendtext = {'Transfer Matrix Algorithm Profile','Lambert Beer Excitation','Test of interpolation'};
legend(legendtext,'FontSize',14)
title('FePt thin film absorption', 'FontSize', fontSize);
xlabel('Position (nm)', 'FontSize', fontSize)
ylabel('Absorption (%.nm$^{-1}$)', 'FontSize', fontSize)
xlim([0,13])
ylim([0,10])


excitation = fluence;    % fluence
tic
[tempMap, deltaTempMap] = H.getTempMap(time,excitation,initTemp);
toc


%% Plot the Results of the Heat Simulations


TFePt = mean(tempMap(:,S.getAllPositionsPerUniqueUnitCell{2},1),2);
TAl  = mean(tempMap(:,S.getAllPositionsPerUniqueUnitCell{1},1),2);



figure(1)%(101)
subplot(2,1,1)
plot(time/u.ps, TFePt , '-b', 'LineWidth', 2);
hold on
plot(time/u.ps, TAl, '-r', 'LineWidth', 2);


xlim([time(1) time(end)]/u.ps)
title('FePt granular film ');
xlabel('Delay [ps]');
ylabel('Temperature [K]');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');
legend('T_{FePt}','T_{Al}', 'FontSize',10);

xmin = 0;
xmax = 40;
subplot(2,1,2)
kk = surf(time/u.ps,distances/u.nm,tempMap(:,:,1)');
set(kk, 'LineStyle', 'none');
ylabel('z (nm)'); xlabel('Delay (ps)');
axis([time(1)/u.ps time(end)/u.ps xmin xmax] )
box; colorbar; colormap jet;
set(gca,'FontSize',14);

box on; grid on;
hold off


%% Checking the initial temperature profile
figure(3)
index0 = finderb(0*u.ps,time);
index1 = finderb(1*u.ps,time);
plot(distances/u.nm,tempMap(index0,:,1),'-r','LineWidth',2)
hold on
plot(distances/u.nm,tempMap(index1,:,1),'-g','LineWidth',2)
xmin = 0;
xmax = 20;

title('FePt granular film ');
xlabel('distance (nm)');
ylabel('Temperature (K)');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');
legend('0 ps','1 ps', 'FontSize',10);
xlim([xmin,xmax]);
% 
% for f = 1:length(fluences)
%     excitation = fluences(f)*0.25;    % fluence
%     tic
%     [tempMap, deltaTempMap] = H.getTempMap(time,excitation,initTemp);
%     toc
% 
%     TFePt = mean(tempMap(:,S.getAllPositionsPerUniqueUnitCell{1},1),2);
%     TMgO  = mean(tempMap(:,S.getAllPositionsPerUniqueUnitCell{2},1),2);
%     
%     % finding the time at which the temperature drops below a threshold
%     index20 = finderA(TFePt,486);
%     index50 = finderA(TFePt,658);
%     index80 = finderA(TFePt,724);
%     times(f,2) = time(index20)/u.ps;
%     times(f,3) = time(index50)/u.ps;
%     times(f,4) = time(index80)/u.ps;
%     temperatures(:,f+1) = TFePt;
%     
%     for t = 1:length(temperatures)
%         if TFePt(t) > Tc
%             magnetization(t,f+1) = 0;
%         else
%             magnetization(t,f+1) = dataM(finderA(dataM(:,1),TFePt(t)),2);            
%         end;
%     end;    
% end
% ExportMatrixToFile(times,'feptTimesThinFilm.dat',{'%t(ps)', 'T(K)'},10)
% ExportMatrixToFile(temperatures,'feptTemperaturesThinFilm.dat',{'%t(ps)', 'T(K)'},10)
% ExportMatrixToFile(magnetization,'feptMagnetizationThinFilm.dat',{'%t(ps)', 'M(t) (%)'},10)
% %%
% figure()
% plot(times(:,1)/(u.mJ/u.cm^2),times(:,2),'s-b','LineWidth',2)
% hold on
% plot(times(:,1)/(u.mJ/u.cm^2),times(:,3),'s-y', 'LineWidth',2)
% plot(times(:,1)/(u.mJ/u.cm^2),times(:,4),'s-r','LineWidth',2)
% leg = legend('t20', 't50', 't80');
% set(leg,'Location','NorthWest')
% xlabel('Incident fluence (mJ/cm²)')
% ylabel('delay (ps)')
% 
% 
% %% Temperature plot
% figure()
% hold on
% cMap =  jet(length(fluences));
% for i = 1:length(fluences)
%     plot(temperatures(:,1),temperatures(:,i+1),'color', cMap(i,:),'LineWidth',2)
% end
% xlabel('delay (ps)')
% ylabel('Temperature (K)')
% xlim([-5,300])
% %% Magnetization plot
% figure()
% hold on
% cMap =  jet(length(fluences));
% for i = 1:length(fluences)
%     plot(temperatures(:,1),magnetization(:,i+1),'color', cMap(i,:),'LineWidth',2)
% end
% xlabel('delay (ps)')
% ylabel('magnetization (%)')
% xlim([-5,300])