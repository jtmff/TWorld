# T-World

T-World is a a human ventricular cardiomyocyte model for simulating membrane voltage, calcium handling, ionic currents, contraction, and beta-adrenergic signalling. The repository includes the core model (Matlab and CellML), helper utilities for parameter handling and post-processing, and example scripts covering common use cases such as baseline pacing, parameter perturbation, EADs, DADs, alternans, restitution, beta-adrenergic stimulation, sex differences, and transmural comparisons.

Publications describing the model will be included here once published.

## Online simulator

You can try out the model using the online interactive simulator here:

**[T-World online simulator](https://t-world-simulator-multipage-production.up.railway.app/)**

## Repository contents

### Core files
The overall logic is similar to our prior model ToR-ORd.

- `model_TWorld.m` — main T-World model implementation (ordinary differential equations).
- `modelRunner.m` — runs the model for one or more beats using a parameter structure and chosen starting state.
- `getDefaultParameters.m` — provides the parameter structure with default values.
- `getStartingState.m` — returns initial conditions for supported model variants, including endocardial, epicardial, midmyocardial, and signalling-including states.
- `getCurrentsStructure.m` — reconstructs and collects currents and key state variables from a simulation into a convenient output structure for plotting and analysis.
- `getConstantsPKASignalling.m` — provides constants used for beta-adrenergic / PKA signalling simulations.
- `DataReporter.m` — helper class used to extract APD in simulations (other functions than getAPD are not to be used for TWorld).

### Example scripts

- `sampleScript_01_basicPacing.m` — basic pacing and plotting.
- `sampleScript_02_perturbingParameters.m` — comparison of control and perturbed parameter sets.
- `sampleScript_03_EADs.m` — example of early afterdepolarization (EAD) induction.
- `sampleScript_04_alternans.m` — example of alternans induced by reduced SERCA uptake.
- `sampleScript_05_DADs.m` — example of delayed afterdepolarization (DAD) induction.
- `sampleScript_06_restitution.m` — S1–S2 restitution protocol.
- `sampleScript_07_BARSdirect.m` — beta-adrenergic effects using directly assigned phosphorylation fractions.
- `sampleScript_08_BARSsignalling.m` — beta-adrenergic effects through the signalling pathway.
- `sampleScript_09_MFcomparison.m` — comparison of female and male model version.
- `sampleScript_10_transmural.m` — comparison of endocardial, epicardial, and midmyocardial cell types.

## Getting started

1. Open MATLAB navigate to the TWorld folder.
2. Choose a starting state with `getStartingState`.
3. Define a structure of simulation parameters, typically including at least:
   - `param.model = @model_TWorld`
   - `param.bcl = 1000`
4. Run the simulation with `modelRunner`.
5. Extract currents and state traces with `getCurrentsStructure`.
6. Plot variables such as membrane potential, intracellular calcium, or active tension.

## Minimal example

```matlab
clear

param.model = @model_TWorld;
param.bcl = 1000; % basic cycle length in ms (1000 ms = 1 Hz)

X0 = getStartingState('TW_endo');
options = [];

beats = 5;
ignoreFirst = 0;

[time, X] = modelRunner(X0, options, param, beats, ignoreFirst);
currents = getCurrentsStructure(time, X, param, 0);

figure;
plot(currents.time, currents.V);
xlabel('Time (ms)');
ylabel('Membrane potential');
```

## Common outputs

The `currents` structure returned by `getCurrentsStructure` can be used directly for plotting and analysis. Commonly used fields include:

- `currents.time`
- `currents.V`
- `currents.Ca_i`
- `currents.Land_Ta`

For example:

```matlab
figure;
plot(currents.time, currents.Ca_i * 1e6);
xlabel('Time (ms)');
ylabel('[Ca]_i (nM)');

figure;
plot(currents.time, currents.Land_Ta);
xlabel('Time (ms)');
ylabel('Active tension');
```

## Model capabilities

T-World supports a range of simulation scenarios (illustrated in sample scripts) through its parameter structure, including:

- ventricular cell-type selection (`endo`, `epi`, `mid`)
- scaling of major ionic currents and fluxes through multipliers (INa_Multiplier for fast sodium current, ICaL_Multiplier for L-type calcium current, Jup_Multiplier for SERCA, etc.)
- changes in extracellular ion concentrations (`nao`, `cao`, `ko`)
- sex-specific parameterizations (`sex = 'F'` or `sex = 'M'`)
- direct control of beta-adrenergic / PKA phosphorylation targets
- simulation of the signalling pathway via cAMP/PKA dynamics
- pacing protocols ranging from regular stimulation to restitution and quiescent intervals

## License

This project is licensed under the **GNU General Public License v3.0**. See the `LICENSE` file for details.

## Contact

Copyright (C) 2026 Jakub Tomek.
Contact: jakub.tomek.mff@gmail.com
