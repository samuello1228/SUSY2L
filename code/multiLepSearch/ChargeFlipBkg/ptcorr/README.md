# README: ptCorr
Last updated 1 March 2017

Part of the [ChargeFlipBkg](../) code package developed by Gabriel Gallardo for the UM/HK EWK SUSY SS 2L analysis.

Compute the absolute pt/energy difference between a reconstructed electron and it's corresponding original-electron-from-Z in MC truth.  
Runs on the converted MC NTuples produced by ConvertNTuples.cxx. ([Documentation](../convert/README.md))

There are two main files in this sub-package:
- [makeDEhistos.C](#makedehistosc)
- [drawDEhistos.py](#drawdehistospy)

Execute both to obtain the necessary histogram for the files in [ChargeFlipBkg/estimate](../estimate/)

Also included is [EnergyCorr.C](#energycorrc), which computes only the energy difference (as opposed to both energy and pt) on UM/HK SUSY MC Ntuples

## makeDEhistos.C
Finds \Delta E and \Delta p_T of reconstructed electrons where \Delta p_T is defined as:  
\Delta p_T = (reconstructed p_T) - (p_T of original electron from Z)  
\Delta E is defined similarly.

### Execution
`root -l -b -q "makeDEhistos.C(\"outputDir\", \"inFileList.txt\")"`
- `outputDir`: Specify output directory
- `inFileList.txt`: Specify full path and filename of input converted NTuples

*Note:* Compilation is not necessary as the converted NTuples are small enough such that this script runs in a sufficiently reasonable amount of time. In addition, the necessary #include statements are not included.

### Options
- `bool signalOnly`: `true` for Signal, `false` for LooseBaseline
- `bool applyPRW`: `true` to apply MCPileupWeight, `false` otherwise

## drawDEhistos.py
- Draws histograms created by `makeDEhistos.C`
- Projects 3D histograms to 2D histograms, among which `hDPTflipped_pxy` will be used by the scripts in [ChargeFlipBkg/estimate](../estimate/)

### Execution
`./drawDEhistos.py dEhistos.root`  
OR  
`python drawDEhistos.py dEhistos.root`

## EnergyCorr.C
### Execution:
- Input the file list of UM/HK NTuples as the argument to the loadData() function called inside EnergyCorr()
- Execute `root -l -b- q EnergyCorr.C`
