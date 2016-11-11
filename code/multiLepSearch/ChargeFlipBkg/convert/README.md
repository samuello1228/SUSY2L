# README: ConvertNTuples 
Last update: 7 September 2016

Part of the [ChargeFlipBkg](../) code package developed by Gabriel Gallardo for the UM/HK EWK SUSY SS 2L analysis. 

Written by Gabriel Gallardo, Jul 2016.

Takes UM/HK SUSY NTuples as input and creates NTuples of an [acceptable format](#output-tree-structure) for the [EGamma chargeMisID tool](../chargeMisID/EG/) and the [pt correction tool](../ptcorr/)

## Execution

### STEP 1: List files in text file 
Do not specify full path, just the filenames
```sh
find /path/to/root/files > convert/myFiles.txt
```

### STEP 2: Check that the makefile is setup properly
Example entry:
```sh
mc: $(PROG) myFiles.txt
    ./$(PROG) myFiles.txt /full/path/to/destination/
```

### STEP 3: Convert NTuples
Following the previous example:
`make mc`

## Output Tree structure
```c++
float MCEvtWeight;
float MCPileupWeight; 
float Zcand_M; 
unsigned long int trigCode;
unsigned int cuts;

// Electron 1
int elCand1_charge; 
float elCand1_pt; 
float elCand1_cl_eta; 
float elCand1_phi;
float elCand1_E;
int elCand1_ID;

// Electron 2
int elCand2_charge; 
float elCand2_pt; 
float elCand2_cl_eta; 
float elCand2_phi;
float elCand2_E;
int elCand2_ID;

// Variables for cuts (SUSY)
int elCand1_flag;
int elCand2_flag; 

// Variables for cuts (both)
float elCand1_d0significance;  // = d0 / sig_d0 
float elCand2_d0significance; 

// Truth pt
float elCand1_truthPt;
float elCand1_truthE;
float elCand2_truthPt;
float elCand2_truthE;

// Original electron pt
float elCand1_origPt;
float elCand1_origE;
int elCand1_origCharge;
float elCand2_origPt;
float elCand2_origE;
int elCand2_origCharge;

float elCand1_dRwOrig;
float elCand2_dRwOrig;
```