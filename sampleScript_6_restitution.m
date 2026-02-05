% A script demonstrating the extraction of S1-S2 restitution curve.
% 
clear

param.model = @model_TWorld;
param.bcl = 1000;
param.PKA_P = zeros(1,8);

X0 = getStartingState('TW_endo');

% APD of S1 and S2 beats is taken at a fixed threshold. Otherwise, APD measurements are confounded
% by differences in peaks of S1 and S2 (or diffferent S2) beats.
threshold = -70; 

% The model is first pre-paced for 100 beats
options = [];
beats = 100;
ignoreFirst = beats-1;
[time1Hz, X1Hz] = modelRunner(X0, options, param, beats, ignoreFirst);

newX0 = X1Hz{1}(end, :);

% For a range of S2 intervals, two additional beats are simulated (the last S1 and the S2 one)
S2interval = [200:5:param.bcl]; % here S2 interval is 200, 205,. 210, ..., 1000

for i = 1:length(S2interval)
    parameterS1S2_1000.bcl = S2interval(i);
    [timeS1S2_1000_2{i,1}, XS1S2_1000_2{i,1}] = modelRunner(newX0, options, parameterS1S2_1000, 2, 0); % 2 beats - one as S1, one as S2.
    currentsS1S2{i} = getCurrentsStructure(timeS1S2_1000_2{i,1}, XS1S2_1000_2{i,1}, parameterS1S2_1000, 0);
end

apdS1S2_1000 = nan(size(S2interval));
diS1S2_1000 = nan(size(S2interval));

for i = 1:length(S2interval)
    % For each coupling interval, we get diastolic interval
    % (exists if and only if there are three segments under
    % APD90 threshold and then it's the middle one's length) and
    % APD of the 2nd spike
    timeS2beats = currentsS1S2{i}.time; % [timeS1S2_1000{i,1}{1};S2interval(i)+timeS1S2_1000{i,2}{1}];
    vS2beats = currentsS1S2{i}.V; %[XS1S2_1000{i,1}{1}(:,1);XS1S2_1000{i,2}{1}(:,1)];


    s = bwconncomp(vS2beats > threshold);
    if (s.NumObjects > 2)
        disp('too many objects found');
    end

    if s.NumObjects == 1 % no 2nd beat was evoked
        apdS1(i) = NaN;
        apdS2(i) = NaN;
        di(i) = 0;
    else
        try
            t1 = timeS2beats(s.PixelIdxList{1});
            t2 = timeS2beats(s.PixelIdxList{2});
            v1 = vS2beats(s.PixelIdxList{1});
            v2 = vS2beats(s.PixelIdxList{2});
        catch
            error(['Crashed in ' num2str(i)]);
        end

        if max(v2) < 0 % 2nd beat has insufficient amplitude (peak < 0 mV), we don't consider this an AP
            apdS1(i) = t1(end) - t1(1);
            apdS2(i) = NaN;
            di(i) = 0;
        else
            t1Start = timeS2beats(s.PixelIdxList{1}(1));
            t1Pre = timeS2beats(s.PixelIdxList{1}(1)-1);
            t1End = timeS2beats(s.PixelIdxList{1}(end));
            t1Post = timeS2beats(s.PixelIdxList{1}(end)+1);
            v1Start = vS2beats(s.PixelIdxList{1}(1));
            v1Pre = vS2beats(s.PixelIdxList{1}(1)-1);
            v1End = vS2beats(s.PixelIdxList{1}(end));
            v1Post = vS2beats(s.PixelIdxList{1}(end)+1);
            partPreStart = (threshold - v1Start)/(v1Pre - v1Start);
            partPostEnd = (v1End - threshold )/(v1End - v1Post);
            apdS1(i) = t1(end) - t1(1) + (t1Start - t1Pre) * partPreStart + (t1Post - t1End) * partPostEnd;

            t2Start = timeS2beats(s.PixelIdxList{2}(1));
            t2Pre = timeS2beats(s.PixelIdxList{2}(1)-1);
            t2End = timeS2beats(s.PixelIdxList{2}(end));
            t2Post = timeS2beats(s.PixelIdxList{2}(end)+1);
            v2Start = vS2beats(s.PixelIdxList{2}(1));
            v2Pre = vS2beats(s.PixelIdxList{2}(1)-1);
            v2End = vS2beats(s.PixelIdxList{2}(end));
            v2Post = vS2beats(s.PixelIdxList{2}(end)+1);
            partPreStart = (threshold - v2Start)/(v2Pre - v2Start);
            partPostEnd = (v2End - threshold )/(v2End - v2Post);
            apdS2(i) = t2(end) - t2(1) + (t2Start - t2Pre) * partPreStart + (t2Post - t2End) * partPostEnd;

            di(i) = t2(1) - t1(end) - (t2Start - t2Pre) * partPreStart - (t2Post - t2End) * partPostEnd;
        end
    end

end

apdS2_control = apdS2;
dvdt1000 = diff(apdS2)./diff(S2interval);
peakDVDT_control = max(dvdt1000);
apdControl = apdS1(end);
di_control = di;

figure(1); clf;
yyaxis left
plot(di, apdS2);
ylabel('APD (ms)')
yyaxis right
plot(di(1:end-1), dvdt1000);
ylabel('Restitution slope');
xlabel('Diastolic interval (ms)');

% (peak slope shows a tiny difference compared to the paper due to a small change in starting state)