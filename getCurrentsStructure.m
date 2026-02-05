%     Cardiac model T-World
%     Copyright (C) 2026 Jakub Tomek. Contact: jakub.tomek.mff@gmail.com
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

function currents = getCurrentsStructure(time, X, parameters, ignoreFirstSpikes)
% A function which computes currents and extract other key state variables,
% collecting these in an output structure.
%
% INPUTS:
% time - a cell array of time vectors (one cell per beat), as returned by
% modelRunner
%
% X - a cell array of matrices of state variables over time (one cell per
% beat), as returned by modelRunner
%
% parameters - passed to the model and specify, e.g. scaling of
% conductances or which implementation of a current is used. Should be the
% same as passed to modelRunner
%
% ignoreFirstSpikes - specifies how many action potentials are to be ignored out of beats recorded.
%
% OUTPUTS:
% currents - a structure where its fields contain the current recordings,
% as well as key state variables (V, Cai, Cass, etc.). Time is also stored here, to allow
% convenient plotting/processing using only the currents structure (e.g. "plot(currents.time, currents.V)"). These variables are not
% anymore separated into separate cells per beat, but are in a single
% vector containing all the beats.

%% Parameters are set here
parameters = getDefaultParameters(parameters);
% Cell type (endo/epi/mid)
cellType = parameters.cellType;

% Extracellular concentrations
nao = parameters.nao;
cao = parameters.cao;
ko = parameters.ko;

% Localization of ICaL and NCX: the fraction in junctional subspace
ICaL_fraction_junc = parameters.ICaL_fraction_junc;
INaCa_fraction_junc = parameters.INaCa_fraction_junc;

% Current multipliers
INa_Multiplier = parameters.INa_Multiplier;
ICaL_Multiplier = parameters.ICaL_Multiplier;
Itos_Multiplier = parameters.Itos_Multiplier;
Itof_Multiplier = parameters.Itof_Multiplier;
INaL_Multiplier = parameters.INaL_Multiplier;
IKr_Multiplier = parameters.IKr_Multiplier;
IKs_Multiplier = parameters.IKs_Multiplier;
IK1_Multiplier = parameters.IK1_Multiplier;
IKb_Multiplier = parameters.IKb_Multiplier;
INaCa_Multiplier = parameters.INaCa_Multiplier;
INaK_Multiplier = parameters.INaK_Multiplier;
INab_Multiplier = parameters.INab_Multiplier;
ICab_Multiplier = parameters.ICab_Multiplier;
IpCa_Multiplier = parameters.IpCa_Multiplier;
IClCa_Multiplier =  parameters.IClCa_Multiplier;
IClb_Multiplier =  parameters.IClb_Multiplier;
Jrel_Multiplier = parameters.Jrel_Multiplier;
Jup_Multiplier = parameters.Jup_Multiplier;
sex = parameters.sex;
PKA_P = parameters.PKA_P;
Whole_cell_PP1 = parameters.Whole_cell_PP1;
runSignallingPathway = parameters.runSignallingPathway;
ISO = parameters.ISO;

% Pre-calculating signalling constants if signalling is active
if (runSignallingPathway == 1)
    constantsSig = getConstantsPKASignalling(ISO);
else
    constantsSig = [];
end

% An array of extra parameters (if user-defined, this can be a cell array
% as well), passing other parameters not defined otherwise in this script
extraParams = parameters.extraParams;

% There may be parameters defining clamp behaviour
vcParameters = parameters.vcParameters;
apClamp = parameters.apClamp;

% Stimulus parameters
stimAmp = parameters.stimAmp;
stimDur = parameters.stimDur;

model = parameters.model;

%% preallocation of currents and other recorded variables
nPoints = length(cell2mat(time((ignoreFirstSpikes+1):end)));

currents.time = zeros(nPoints, 1);
currents.V = zeros(nPoints, 1);
currents.Ca_i = zeros(nPoints, 1);
currents.Ca_junc = zeros(nPoints, 1);
currents.Ca_sl = zeros(nPoints, 1);
currents.Ca_SR = zeros(nPoints, 1);
currents.Na_i = zeros(nPoints, 1);

currents.INa = zeros(nPoints, 1);
currents.INaL = zeros(nPoints, 1);
currents.Ito = zeros(nPoints, 1);
currents.ICaL = zeros(nPoints, 1);
currents.IKr = zeros(nPoints, 1);
currents.IKs = zeros(nPoints, 1);
currents.IK1 = zeros(nPoints, 1);
currents.INaCa_junc = zeros(nPoints, 1);
currents.INaCa_sl = zeros(nPoints, 1);
currents.INaCa = zeros(nPoints, 1);
currents.INaK = zeros(nPoints, 1);
currents.IKb = zeros(nPoints, 1);
currents.INab = zeros(nPoints, 1);
currents.ICab = zeros(nPoints, 1);
currents.IpCa = zeros(nPoints, 1);
currents.Jrel = zeros(nPoints, 1);
currents.Jleak = zeros(nPoints, 1);
currents.IClCa = zeros(nPoints, 1);
currents.IClCa_junc = zeros(nPoints, 1);
currents.IClCa_sl = zeros(nPoints, 1);
currents.IClb = zeros(nPoints, 1);
currents.ICaLtot = zeros(nPoints, 1);
currents.Jup = zeros(nPoints, 1);
currents.ICaL_junc = zeros(nPoints, 1);
currents.ICaL_sl = zeros(nPoints, 1);
currents.Jrel_ICaLdep = zeros(nPoints, 1);
currents.Itos = zeros(nPoints, 1);
currents.Itof = zeros(nPoints, 1);
currents.Land_Ta = zeros(nPoints, 1);
currents.fINa_PKA = zeros(nPoints, 1);
currents.fICaL_PKA = zeros(nPoints, 1);
currents.fINaK_PKA = zeros(nPoints, 1);
currents.fIKs_PKA = zeros(nPoints, 1);
currents.fPLB_PKA = zeros(nPoints, 1);
currents.fTnI_PKA = zeros(nPoints, 1);
currents.fMyBPC_PKA = zeros(nPoints, 1);

%% Using model file to get variables of interest.
i = 1;
for iBeat = 1:size(time, 1)
    if (iBeat <= ignoreFirstSpikes)
        continue
    end
    % currents.time = [currents.time; (iBeat-1)*parameters.bcl + time{iBeat}];

    for j=1:size(X{iBeat},1)
        IsJs=model(time{iBeat}(j),X{iBeat}(j,:)',0, cellType, ICaL_Multiplier, ...
            INa_Multiplier, Itos_Multiplier, Itof_Multiplier, INaL_Multiplier, IKr_Multiplier, IKs_Multiplier, IK1_Multiplier, IKb_Multiplier,INaCa_Multiplier,...
            INaK_Multiplier, INab_Multiplier, ICab_Multiplier, IpCa_Multiplier, IClCa_Multiplier, IClb_Multiplier, Jrel_Multiplier,Jup_Multiplier,nao,cao,ko,ICaL_fraction_junc,INaCa_fraction_junc, PKA_P, Whole_cell_PP1, runSignallingPathway, constantsSig, sex, stimAmp, stimDur, vcParameters, apClamp, extraParams);

        currents.time(i) = (iBeat-1)*parameters.bcl + time{iBeat}(j);
        currents.V(i) = X{iBeat}(j,1);
        currents.Ca_i(i) = X{iBeat}(j,8);
        currents.Ca_junc(i) = X{iBeat}(j,6);
        currents.Ca_sl(i) = X{iBeat}(j,7);
        currents.Ca_SR(i) = X{iBeat}(j,10);
        currents.Na_i(i) = X{iBeat}(j,4);

        currents.INa(i)=IsJs(1);
        currents.INaL(i)=IsJs(2);
        currents.Ito(i)=IsJs(3);
        currents.ICaL(i)=IsJs(4);
        currents.IKr(i)=IsJs(5);
        currents.IKs(i)=IsJs(6);
        currents.IK1(i)=IsJs(7);
        currents.INaCa_junc(i)=IsJs(8);
        currents.INaCa_sl(i)=IsJs(9);
        currents.INaCa(i) = currents.INaCa_sl(i) + currents.INaCa_junc(i);
        currents.INaK(i)=IsJs(10);
        currents.IKb(i)=IsJs(11);
        currents.INab(i)=IsJs(12);
        currents.ICab(i)=IsJs(13);
        currents.IpCa(i)=IsJs(14);
        currents.Jrel(i) = IsJs(15);
        currents.Jleak(i) = IsJs(16);
        currents.IClCa(i) = IsJs(17);
        currents.IClCa_junc(i) = IsJs(18);
        currents.IClCa_sl(i) = IsJs(19);
        currents.IClb(i) = IsJs(20);
        currents.ICaLtot(i) = IsJs(21);
        currents.Jup(i) = IsJs(22);
        currents.ICaL_junc(i) = IsJs(23);
        currents.ICaL_sl(i) = IsJs(24);
        currents.Jrel_ICaLdep(i) = IsJs(25);
        currents.Itos(i) = IsJs(26);
        currents.Itof(i) = IsJs(27);
        currents.Land_Ta(i) = IsJs(28);
        currents.fINa_PKA(i) = IsJs(29);
        currents.fICaL_PKA(i) = IsJs(30);
        currents.fINaK_PKA(i) = IsJs(31);
        currents.fIKs_PKA(i) = IsJs(32);
        currents.fPLB_PKA(i) = IsJs(33);
        currents.fTnI_PKA(i) = IsJs(34);
        currents.fMyBPC_PKA(i) = IsJs(35);

        i = i + 1; % incrementing the index where everything is stored
    end
end
end

