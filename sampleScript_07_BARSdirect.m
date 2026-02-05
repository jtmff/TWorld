% This script shows how to simulate effect of beta-adrenergic signalling, setting phoshorylation levels
% of targets directly.
%
clear

param.model = @model_TWorld;
param.bcl = 1000; % basic cycle length (1000 ms ~ 1 Hz)

% We set all the phosphorylation fractions to 1 here. You can try setting some to 0 (or any value
% between 0 and 1) to see how the behaviour of the cell changes.
fINa_PKA = 1;
fICaL_PKA = 1;
fINaK_PKA = 1;
fIKs_PKA = 1;
fPLB_PKA = 1;
fTnI_PKA = 1;
fMyBPC_PKA = 1;

param.PKA_P = [fINa_PKA, fICaL_PKA, fINaK_PKA, fIKs_PKA, fPLB_PKA, fTnI_PKA, fMyBPC_PKA];

X0 = getStartingState('TW_endo');

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

