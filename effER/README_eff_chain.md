# Acc x Eff production chain

Run these commands from **effER/**. Every call works on one tree in one year, independently of the others.
Supported names are exact and case-sensitive:

- trees: **ntmix_X3872**, **ntmix_PSI2S**, **ntKp**, **ntphi**, **ntKstar**
- years: **ppRef** and the PbPb years (**PbPb23**, **PbPb24**, ...)
- map dimensions: **0D**, **1D**, **2D**
- efficiency weights: **raw** and **Prediction** (the classifier score), or another selection variable of the
  validation weights; never a variable that is an axis of the map (**Bpt**, **By**), see section 1
- reading methods: **sPlot**, **mWindow**

There are no aliases, fallback samples, or name normalization.

## Inputs

Only binned fit products enter this chain, from **../fitER/ROOTfiles/<YEAR>/** for every year: the single-year
fit (`roofitB.C`) or, for a PbPb year, the per-year files of the combined fit (`combinedANA`), which have the same
names and replace them.

~~~text
../fitER/ROOTfiles/<YEAR>/fitResults_<TREE>_Bpt_<YEAR>.root
../fitER/ROOTfiles/<YEAR>/fitResults_<TREE>_<VAR>_<YEAR>.root
../fitER/ROOTfiles/<YEAR>/mcFitResults_<TREE>_<VAR>_<YEAR>.root
~~~

The first file supplies **inputMC**, the reconstructed selection, and the Bpt analysis bins to map production.
The second supplies fitted data yields, stored data, snapshots, and models to map reading. The third is the
MC-only closure input.

**inputMC** is a scored MC sample (selectionER), which holds only the reconstructed tree. The generator tree
**ntGen** is read from the flat file named by its metadata (**input_file**). A fit made on a flat MC file
(no metadata) breaks the map production there.

Validation supplies the weights of each tree and year:

~~~text
../plotER/Validation/WEIGHTS/<DATA TREE>_<YEAR>_<PARTICLE>_weight.root
    ntmix_<YEAR>_X3872, ntmix_<YEAR>_PSI2S, ntKp_<YEAR>_Bp, ntphi_<YEAR>_Bs, ntKstar_<YEAR>_B0
~~~

The map code reads **hWeight_<WEIGHT>** and its stored **weightExpression_<WEIGHT>**. A validation weight
modifies only the selected reconstructed numerator. Generator and acceptance counts are never variable-reweighted.
Every generated, accepted, and selected MC count uses **pThatreweight**.

The weight therefore tests how well the MC models the selection: the data/MC ratio of a selection variable, the
classifier score **Prediction** first. A variable that is an axis of the map (**Bpt**, **By**) must not be used: in a
fine map bin the numerator-only weight just scales the efficiency by the data/MC ratio of the spectrum, so the
corrected yield moves by its inverse (B+ PbPb: up to -23% at 30-50 GeV), which measures the spectrum difference, not
an efficiency mismodeling. The data distribution inside each analysis bin is already taken into account by the sPlot
reading of the 2D map.

## Settings of each tree and year

`GetAccEffSetup` (`aux/uti.h`) holds them; a tree or year it does not list breaks the lookup:

| tree | year | generator cut | track acceptance | map pT bins (number) | map abs(y) bins | mass window |
|---|---|---|---|---|---|---|
| ntmix_X3872 | ppRef | 7.5 < pT < 50, abs(y) < 2.4 | 2 tracks, pT > 0.5 | 43, 0.5-5 GeV wide | 0, 0.6, 1.2, 1.6, 2.0, 2.4 | +-20 MeV |
| ntmix_X3872 | PbPb | 15 < pT < 50, abs(y) < 1.6 | 2 tracks, pT > 0.9 | 15 | 0, 0.8, 1.6 | +-20 MeV |
| ntmix_PSI2S | as X(3872) | | | | | +-15 MeV |
| ntKp | all | 7.5 < pT < 50, abs(y) < 2.4 | 1 track, pT > 0.5 (ppRef), 0.9 (PbPb) | 14 | 0, 0.5, 1.0, 1.5, 2.0, 2.4 | +-60 MeV |
| ntKstar | all | as ntKp | 2 tracks | 14 | as ntKp | +-40 MeV |
| ntphi | all | as ntKp | 2 tracks | 14 | as ntKp | +-25 MeV |

Muon acceptance (all trees): the soft-muon acceptance in pT and abs(eta) < 2.4. Tracks: abs(eta) < 2.4 and the pT
threshold of the reconstruction. The B0 and Bs mass windows are the Punzi windows of selectionER (about 2 sigma).
The B+ window is +-60 MeV, wider than its +-40 MeV Punzi window: the window drops the mass tails, whose candidates
have another pT, |y| mix, so the MC closure of the reading at 7.5-10 GeV is 4.4% off with +-40 MeV, 2.5% with
+-50 MeV and 1.4% with +-60 MeV (within 0.5% above 10 GeV). Wider windows close better in MC but bring more
background into the data.

## Outputs

Everything of one tree in one year is in **output/<TREE>/<YEAR>/**:

~~~text
output/<TREE>/<YEAR>/ROOTs/            maps, corrected yields, comparisons, closure
output/<TREE>/<YEAR>/2Dmaps/           2D map PDFs
output/<TREE>/<YEAR>/ACCxEFF_plots/    <1/(Acc x Eff)> per analysis bin
output/<TREE>/<YEAR>/systematicFILES/  comparison PDFs
output/<TREE>/<YEAR>/closure/          closure PDFs
~~~

## 1. Produce the maps

~~~cpp
accXeff_MAPS(TREE, YEAR, DIMENSION, WEIGHT)
~~~

- **2D**: fine pT-versus-abs(y) map from MC counts; ROOT and PDF output.
- **1D**: pT axis of those counts, integrated over the abs(y) range; ROOT output only.
- **0D**: one direct **Nsel/Nacc** value in every Bpt analysis bin; ROOT output only.

~~~bash
# for every tree and year, e.g. ntKp PbPb24
root -l -b -q 'accXeff_MAPS.C("ntKp","PbPb24","0D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntKp","PbPb24","1D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntKp","PbPb24","2D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntKp","PbPb24","2D","Prediction")'
~~~

~~~text
output/<TREE>/<YEAR>/ROOTs/<TREE>_<YEAR>_<DIMENSION>map_ACCxEFF_<WEIGHT>.root
output/<TREE>/<YEAR>/ROOTs/<TREE>_<YEAR>_<DIMENSION>map_ACC.root
~~~

## 2. Read the maps in each analysis bin

~~~cpp
accXeff_READ(TREE, YEAR, VAR, METHOD, DIMENSION, WEIGHT)
~~~

**sPlot** reloads **nominalPars_binN** from the binned fit of the year, floats its yields (signal, combinatorial
background and, for ntKp, the partially reconstructed B) for the sPlot covariance, and computes:

~~~text
corrected yield = sum_i signal_sWeight_i / AccEff_i
~~~

**mWindow** averages **1/AccEff** over all selected data candidates in the resonance window (table above), then
multiplies it by the binned fitted signal yield.

~~~bash
root -l -b -q 'accXeff_READ.C("ntKp","PbPb24","Bpt","sPlot","0D","raw")'
root -l -b -q 'accXeff_READ.C("ntKp","PbPb24","Bpt","sPlot","1D","raw")'
root -l -b -q 'accXeff_READ.C("ntKp","PbPb24","Bpt","sPlot","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntKp","PbPb24","Bpt","mWindow","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntKp","PbPb24","Bpt","sPlot","2D","Prediction")'
~~~

The raw 2D sPlot case is nominal, but it uses the same code path as every other combination.

~~~text
output/<TREE>/<YEAR>/ROOTs/<TREE>_<YEAR>_<VAR>_<DIMENSION>_<WEIGHT>_<METHOD>_CorrectedYields.root
~~~

## 3. Produce the two systematic comparisons

The comparison receives explicit **{DIMENSION, WEIGHT, METHOD}** triples:

~~~cpp
accXeff_COMPARISONS(TREE, YEAR, VAR, CASES, TAG, SAVE_ROOT)
~~~

The raw 2D sPlot case is the reference. A PDF is always written. With **SAVE_ROOT=true**, the comparison ROOT
file stores both `<1/(Acc x Eff)>` and the corrected yield for every named case, together with `referenceCase` and
`case<N>` metadata. `resultER` uses the corrected yields to propagate paired X(3872) and psi(2S) variations through
the final ratio. With **SAVE_ROOT=false**, any existing systematic ROOT artifact remains untouched, so an exploratory
plot cannot replace a resultER input.

~~~bash
root -l -b -q 'accXeff_COMPARISONS.C("ntKp","PbPb24","Bpt",{{"0D","raw","sPlot"},{"1D","raw","sPlot"},{"2D","raw","sPlot"},{"2D","raw","mWindow"}},"METHODS_comparison",true)'
root -l -b -q 'accXeff_COMPARISONS.C("ntKp","PbPb24","Bpt",{{"2D","raw","sPlot"},{"2D","Prediction","sPlot"}},"DATA_MC_AGREEMENT_reweight_comparison",true)'
~~~

~~~text
output/<TREE>/<YEAR>/systematicFILES/METHODS_comparison_<TREE>_<YEAR>_<VAR>.pdf
output/<TREE>/<YEAR>/ROOTs/METHODS_comparison_<TREE>_<YEAR>_<VAR>.root
output/<TREE>/<YEAR>/systematicFILES/DATA_MC_AGREEMENT_reweight_comparison_<TREE>_<YEAR>_<VAR>.pdf
output/<TREE>/<YEAR>/ROOTs/DATA_MC_AGREEMENT_reweight_comparison_<TREE>_<YEAR>_<VAR>.root
~~~

resultER reads only these exact tagged ROOT filenames. Other tags do not enter uncertainty propagation.

## 4. Run the MC-only closure

**Closure_methods.C** never opens data. It reads the raw maps and the **mcFitResults** companion of the year.
All candidates in each stored MC dataset are selected and generator matched already.

~~~bash
root -l -b -q 'Closure_methods.C("ntKp","PbPb24","Bpt")'
~~~

The single closure PDF contains:

- raw 0D using all selected signal MC
- raw 1D using all selected signal MC
- raw 2D using all selected signal MC
- raw 2D using the fitted-MC-mean mass window
- raw 2D using the saved MC model and sPlot

All five cases use pThat weights. No validation-derived variable weight is used in closure.

~~~text
output/<TREE>/<YEAR>/closure/closure_<TREE>_<YEAR>_<VAR>.pdf
output/<TREE>/<YEAR>/ROOTs/closure_<TREE>_<YEAR>_<VAR>.root
~~~

## 5. resultER contract

~~~text
../effER/output/ntmix_X3872/<YEAR>/ROOTs/ntmix_X3872_<YEAR>_<VAR>_2D_raw_sPlot_CorrectedYields.root
../effER/output/ntmix_PSI2S/<YEAR>/ROOTs/ntmix_PSI2S_<YEAR>_<VAR>_2D_raw_sPlot_CorrectedYields.root
~~~

Run the ratio-level uncertainty producer from `resultER/`:

~~~bash
root -l -b -q 'ntmix_UNCpropagator.C("ppRef","Bpt")'
root -l -b -q 'ntmix_UNCpropagator.C("ppRef","nChargedTracks")'
~~~

For each Acc x Eff source, **ntmix_UNCpropagator.C** matches the same stored case for X(3872) and psi(2S), forms the
varied corrected-yield ratio, and assigns the largest absolute variation from the nominal ratio in each analysis bin.
It combines the two particle fit-model uncertainties as independent because the fit artifacts currently contain only
the leading per-particle uncertainty. The fit-ratio, data-MC-ratio, and method-ratio sources are then added in
quadrature. The result is:

~~~text
resultER/output_ntmix/systematicFILES/ntmix_ratioUnc_X3872_over_PSI2S_<YEAR>_<VAR>.root
~~~

`Ratio_CorrectedYields.C` reads this ratio-level total. The procedure is independent of whether `<VAR>` is pT,
multiplicity, rapidity, or centrality; the binning comes from the stored corrected-yield histograms.
