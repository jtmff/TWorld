% This script shows how to simulate effect of beta-adrenergic signalling, setting phoshorylation levels
% of targets through cAMP signalling.
%
% You may notice there is a small difference between the outputs here and of the sample script 7.
% This arises from a slight difference in phosphatase activity with full BARS activation.
% If you set param.Whole_cell_PP1 = 0.11 in sampleScript_7, you will get a virtually perfect
% agreement.
clear

param.model = @model_TWorld;
param.bcl = 1000; % basic cycle length (1000 ms ~ 1 Hz)
param.runSignallingPathway = 1; % setting signalling on
param.ISO = 1; % 1 uM isoproterenol

% using a starting state that includes signalling state variables
X0 = getStartingState('TW_endo_sig');

options = [];
beats = 250;
ignoreFirst = 249;

tic
% simulate model
[time1Hz, X1Hz] = modelRunner(X0, options, param, beats, ignoreFirst);
% and extract currents
currents1Hz = getCurrentsStructure(time1Hz, X1Hz, param, 0);
toc

%% Plotting of key variables
figure(1); clf;
plot(currents1Hz.time, currents1Hz.V);
xlabel('Time (ms)'); ylabel('Membrane potential');

figure(2); clf;
plot(currents1Hz.time, currents1Hz.Ca_i * 1e6); % 1e6 to convert to nM
xlabel('Time (ms)'); ylabel('[Ca]_i (nM)');

figure(3); clf;
plot(currents1Hz.time, currents1Hz.Land_Ta);
xlabel('Time (ms)'); ylabel('Active tension');
