#!/bin/bash

LooseRates="../results/0805-MC-loose-80100/output.root"
LooseDPT="../results/0808-MC-EG-loose/dEhistos.root"
SignalRates="../results/0805-MC-signal-80100/output.root"
SignalDPT="../results/0808-MC-EG-signal/dEhistos.root"

tag="0810-ZeeMC-"
inFileList="../common/inFileList-MC.txt"

# folder=$tag"loose-old-pwt"
# root -l -b -q "saveSSprediction.C+(\"$folder\", \"$inFileList\", \"$LooseRates\", \"$LooseDPT\", \"el=l\")"
# decorations=$folder"/misIDdecorations.root"
# root -l -b -q "extractSSprediction.C+(\"$folder\", \"$inFileList\", \"$decorations\", \"el=l\")"

folder=$tag"loose-new-pwt-oldDpt"
root -l -b -q "getSSPrediction.C+(\"$folder\", \"$inFileList\", \"$LooseRates\", \"$LooseDPT\", \"el=l\")"


# folder=$tag"signal-old-pwt"
# root -l -b -q "saveSSprediction.C+(\"$folder\", \"$inFileList\", \"$SignalRates\", \"$SignalDPT\", \"el=s\")"
# decorations=$folder"/misIDdecorations.root"
# root -l -b -q "extractSSprediction.C+(\"$folder\", \"$inFileList\", \"$decorations\", \"el=s\")"

folder=$tag"signal-new-pwt-oldDpt"
root -l -b -q "getSSPrediction.C+(\"$folder\", \"$inFileList\", \"$SignalRates\", \"$SignalDPT\", \"el=s\")"

