% This script shows how SERCA reduction promotes alternans at slower pacing rates.
% The simulations show no alternans at bcl of 400 ms in control model, but alternans in
% SERCA-reduced one.
%
clear

% Define baseline model - all perturbed parameters are to be set to 1 here.
param.model = @model_TWorld;
param.bcl = 400; 
param.Jup_Multiplier = 1;

parameters(1:2) = param;
parameters(2).Jup_Multiplier = 0.65;

X0 = getStartingState('TW_endo');

options = [];

% We want to simulate 250 beats, plotting the last four.
beats = 250;
ignoreFirst = 246;

parfor iModel = 1:2 % change to parfor if parallel computing toolbox is available
    [time1Hz{iModel}, X1Hz{iModel}] = modelRunner(X0, options, parameters(iModel), beats, ignoreFirst);
    currents1Hz{iModel} = getCurrentsStructure(time1Hz{iModel}, X1Hz{iModel}, parameters(iModel), 0);
end


%% Plotting of key variables
figure(1); clf;
hold on
for iModel = 1:numel(currents1Hz)
    plot(currents1Hz{iModel}.time, currents1Hz{iModel}.V);
end
hold off
xlabel('Time (ms)'); ylabel('Membrane potential');
legend('Control', 'Perturbed')

figure(2); clf;
hold on
for iModel = 1:numel(currents1Hz)
    plot(currents1Hz{iModel}.time, currents1Hz{iModel}.Ca_i * 1e6);
end
hold off
xlabel('Time (ms)'); ylabel('[Ca]_i (nM)');
legend('Control', 'Low SERCA')
