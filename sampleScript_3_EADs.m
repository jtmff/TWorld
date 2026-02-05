% This script shows how EADs can be evoked.
% 
clear

param.model = @model_TWorld;
param.bcl = 4000; % slow pacing
param.IKr_Multiplier = 0.15; % reduced hERG
param.cao = 2; % elevated extracellular calcium

X0 = getStartingState('TW_endo');

options = [];
beats = 100;
ignoreFirst = 99;

% simulate model
[time1Hz, X1Hz] = modelRunner(X0, options, param, beats, ignoreFirst);
currents1Hz = getCurrentsStructure(time1Hz, X1Hz, param, 0);

% and plot outputs
figure(1);
plot(currents1Hz.time, currents1Hz.V);
xlabel('Time (ms)'); ylabel('Membrane potential');
xlim([0 2000]);
