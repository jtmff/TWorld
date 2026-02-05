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

function [parameters] = getDefaultParameters(parameters)
%GETDEFAULTPARAMETERS Populates the structure of parameters with defaults for those that are
%undefined.

% Cell type (endo/epi/mid)
if (~isfield(parameters,'cellType')) parameters.cellType = 0;end

% If the number of simulated beat is to be printed out.
if (~isfield(parameters,'verbose')) parameters.verbose = 0; end

% Extracellular concentrations
if (~isfield(parameters,'nao')) parameters.nao = 140; end
if (~isfield(parameters,'cao')) parameters.cao = 1.8; end
if (~isfield(parameters,'ko')) parameters.ko = 5; end


% Localization of ICaL and NCX: the fraction in junctional subspace
if (~isfield(parameters,'ICaL_fraction_junc')) parameters.ICaL_fraction_junc = 0.8; end
if (~isfield(parameters,'INaCa_fraction_junc')) parameters.INaCa_fraction_junc = 0.31; end

% Current multipliers
if (~isfield(parameters,'INa_Multiplier')) parameters.INa_Multiplier = 1; end
if (~isfield(parameters,'INaL_Multiplier')) parameters.INaL_Multiplier = 1; end
if (~isfield(parameters,'ICaL_Multiplier')) parameters.ICaL_Multiplier = 1; end
if (~isfield(parameters,'Itos_Multiplier')) parameters.Itos_Multiplier = 1; end
if (~isfield(parameters,'Itof_Multiplier')) parameters.Itof_Multiplier = 1; end
if (~isfield(parameters,'IKr_Multiplier')) parameters.IKr_Multiplier = 1; end
if (~isfield(parameters,'IKs_Multiplier')) parameters.IKs_Multiplier = 1; end
if (~isfield(parameters,'IK1_Multiplier')) parameters.IK1_Multiplier = 1; end
if (~isfield(parameters,'IKb_Multiplier')) parameters.IKb_Multiplier = 1; end
if (~isfield(parameters,'INaCa_Multiplier')) parameters.INaCa_Multiplier = 1; end
if (~isfield(parameters,'INaK_Multiplier')) parameters.INaK_Multiplier = 1; end
if (~isfield(parameters,'INab_Multiplier')) parameters.INab_Multiplier = 1; end
if (~isfield(parameters,'ICab_Multiplier')) parameters.ICab_Multiplier = 1; end
if (~isfield(parameters,'IpCa_Multiplier')) parameters.IpCa_Multiplier = 1; end
if (~isfield(parameters,'IClCa_Multiplier')) parameters.IClCa_Multiplier = 1; end
if (~isfield(parameters,'IClb_Multiplier')) parameters.IClb_Multiplier = 1; end
if (~isfield(parameters,'Jrel_Multiplier')) parameters.Jrel_Multiplier = 1; end
if (~isfield(parameters,'Jup_Multiplier')) parameters.Jup_Multiplier = 1; end

if (~isfield(parameters,'sex')) parameters.sex = 'U'; end % Unisex

% PKA signalling
if (~isfield(parameters,'PKA_P')) parameters.PKA_P = [0, 0, 0, 0, 0, 0, 0]; end 
if (~isfield(parameters,'Whole_cell_PP1')) parameters.Whole_cell_PP1 = 0.13698; end 
if (~isfield(parameters,'runSignallingPathway')) parameters.runSignallingPathway = 0; end 
if (~isfield(parameters,'ISO')) parameters.ISO = 0; end % Nonspecific

% Whole_cell_PP1
% PKA_P

% An array of extra parameters (if user-defined, this can be a cell array
% as well), passing other parameters not defined otherwise in this script
if (~isfield(parameters,'extraParams')) parameters.extraParams = []; end

% There may be parameters defining clamp behaviour
if (~isfield(parameters,'vcParameters')) parameters.vcParameters = []; end
if (~isfield(parameters,'apClamp')) parameters.apClamp = []; end

% Stimulus parameters
if (~isfield(parameters,'stimAmp')) parameters.stimAmp = -53; end
if (~isfield(parameters,'stimDur')) parameters.stimDur = 1; end

% Model
if (~isfield(parameters,'model')) parameters.model = @model_TWorld; end


end

