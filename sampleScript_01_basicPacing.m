% This script runs and plots 5 beats of the endocardial cell at 1 Hz.
%
clear

param.model = @model_TWorld;
param.bcl = 1000; % basic cycle length (1000 ms ~ 1 Hz)

X0 = getStartingState('TW_endo');

options = [];

% Settings below simulate 5 beats, showing all of them. 
beats = 5;
ignoreFirst = 0;

% Parameters below would simulate 100 beats, keeping data for only the last one.
% beats = 100;
% ignoreFirst = beats - 1;

% simulate model
[time1Hz, X1Hz] = modelRunner(X0, options, param, beats, ignoreFirst);
% and extract currents
currents1Hz = getCurrentsStructure(time1Hz, X1Hz, param, 0);

%% Plotting of key variables
figure(1);
plot(currents1Hz.time, currents1Hz.V);
xlabel('Time (ms)'); ylabel('Membrane potential');

figure(2);
plot(currents1Hz.time, currents1Hz.Ca_i * 1e6); % 1e6 to convert to nM
xlabel('Time (ms)'); ylabel('[Ca]_i (nM)');

figure(3);
plot(currents1Hz.time, currents1Hz.Land_Ta);
xlabel('Time (ms)'); ylabel('Active tension');

