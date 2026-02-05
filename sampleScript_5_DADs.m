% A script demonstrating DAD evocation following rapid prepacing at high Cao and beta-AR stimulation
clear

param.model = @model_TWorld;
param.bcl = 400; % fast pre-pacing
param.PKA_P = ones(1,8); % full phosphorylation of PKA (sympathetic) targets
param.cao = 3.25; % high extra-cellular calcium

X0 = getStartingState('TW_endo');

options = [];

% Prepacing the model first using regular pacing at 400 ms bcl
beats = 200;
ignoreFirst = 0;
[time, X] = modelRunner(X0, options, param, beats, ignoreFirst);

% And then adding one beat with long bcl, to create a period of quiescence
param.bcl = 10000;
beats = 1;
ignoreFirst = 0;
[time_DAD, X_DAD] = modelRunner(X{end}(end,:), options, param, beats, ignoreFirst);
currents_DAD = getCurrentsStructure(time_DAD, X_DAD, param, 0);


%% Plotting of key variables
figure(1); clf;
plot(currents_DAD.time, currents_DAD.V);
xlabel('Time (ms)'); ylabel('Membrane potential');

figure(2); clf;
plot(currents_DAD.time, currents_DAD.Ca_i * 1e3); % 1e3 to bring to uM range.
xlabel('Time (ms)'); ylabel('[Ca]_i (uM)');

% there is only one stimulated beat - the second one is spontaneous.