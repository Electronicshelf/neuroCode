SPACE_UNITS = 'µm';
TIME_UNITS = 's';

N_PARTICLES = 178;
N_TIME_STEPS = 34;
N_DIM = 2; % 2D

% Typical values taken from studies of proteins diffusing in membranes:
% Diffusion coefficient
D  = 1e-3; % µm^2/s
% Time step between acquisition; fast acquisition!
dT = 1; % s,

% Area size, just used to disperse particles in 2D. Has no impact on
% analysis.
SIZE = 2; % µm
 ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
 % This works:

ma = ma.addAll(tracks);

% Indeed:
disp(ma)
ma.plotTracks;
ma.labelPlotTracks;
figure
%% plot Specific points 
plot(tracks{4,1}(:,2),tracks{4,1}(:,3))
plot(tracks{4,1}(:,2),tracks{4,1}(:,3))
%ax2 = axes();
%ma.plotTracks(ax2,tracks{4,1}(2:3));
%% plot
figure
ma.plotMSD;

cla
ma.plotMeanMSD(gca,true)

[fo, gof] = ma.fitMeanMSD;
plot(fo)
ma.labelPlotMSD;
legend off