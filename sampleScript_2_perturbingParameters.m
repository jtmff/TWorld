% This script shows how to compare the baseline and perturbed model (half INa and IKr).
% 
clear

% Define baseline model - all perturbed parameters are to be set to 1 here.
param.model = @model_TWorld;
param.bcl = 1000; 
% When perturbing any multipliers, the default model should have these defined as 1.
% (otherwise they are autofilled as empty and lead to the simulation crashing)
param.INa_Multiplier = 1; 
param.IKr_Multiplier = 1;

parameters(1:2) = param;
parameters(2).INa_Multiplier = 0.5;
parameters(2).IKr_Multiplier = 0.5;

X0 = getStartingState('TW_endo');

options = [];

% This runs the model for 100 beats, not storing the first 99, hence storing and plotting the last beat
beats = 100;
ignoreFirst = 99;

for iModel = 1:2 % change to parfor if parallel computing toolbox is available
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
    plot(currents1Hz{iModel}.time, currents1Hz{iModel}.Ca_i * 1e6); % 1e6 to convert to nM
end
hold off
xlabel('Time (ms)'); ylabel('[Ca]_i (nM)');
legend('Control', 'Perturbed')


