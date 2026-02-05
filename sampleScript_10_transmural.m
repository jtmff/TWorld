% This script runs and plots a beat of each of three cell types
clear

param.model = @model_TWorld;
param.bcl = 1000; % basic cycle length (1000 ms ~ 1 Hz)

startingState{1} = getStartingState('TW_endo');
startingState{2} = getStartingState('TW_epi');
startingState{3} = getStartingState('TW_mid');

params(1:3) = param;
params(1).cellType = 0;
params(2).cellType = 1;
params(3).cellType = 2;

options = [];
beats = 100;
ignoreFirst = 99;

parfor i = 1:3
    tic
    % simulate model
    [time1Hz{i}, X1Hz{i}] = modelRunner(startingState{i}, options, params(i), beats, ignoreFirst);
    % and extract currents
    currents1Hz{i} = getCurrentsStructure(time1Hz{i}, X1Hz{i}, params(i), 0);
    toc
end

%% Plotting of key variables
figure(1); clf;
hold on
for i = 1:3
    plot(currents1Hz{i}.time, currents1Hz{i}.V);
end
hold off
legend('Endo','Epi','Mid');
xlabel('Time (ms)'); ylabel('Membrane potential');

figure(2); clf;
hold on
for i = 1:3
    plot(currents1Hz{i}.time, currents1Hz{i}.Ca_i*1e6);
end
hold off
legend('Endo','Epi','Mid');
xlabel('Time (ms)'); ylabel('[Ca]_i (nM)');

figure(3); clf;
hold on
for i = 1:3
    plot(currents1Hz{i}.time, currents1Hz{i}.Land_Ta);
end
hold off
legend('Endo','Epi','Mid');
xlabel('Time (ms)'); ylabel('Active tension');

