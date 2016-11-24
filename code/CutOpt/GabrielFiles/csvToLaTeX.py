#!/usr/bin/env python

# import system, os
from collections import namedtuple
import csv
from operator import attrgetter

try:
    from itertools import imap
except ImportError:  # Python 3
    imap = map

data = []

with open("checks.csv") as file:
	reader = csv.reader(file)
	entry = namedtuple("entry", next(reader))

	for row in imap(entry._make, reader):
		data.append(row)

dm020 = []
dm050 = []
dm100 = []
dm200 = []

massPoints = (dm020, dm050, dm100, dm200)

for row in data:
    if row.dm == "20":
        dm020.append(row)
    elif row.dm == "50":
        dm050.append(row)
    elif row.dm == "100":
        dm100.append(row)
    elif row.dm == "+":
        dm200.append(row)
    else:
        print "Unknown mass point"

for p in massPoints:
    p.sort(key=lambda x:(x.Flavor, x.ISR, int(x.NTrees), int(x.NodeSize)))
    # p.sort(key=attrgetter("NodeSize","NTrees", "Flavor"))

for row in dm020:
    print "%s, %s, %s, %s" % (row.Flavor, row.ISR, row.NTrees, row.NodeSize)

def beginTable(dm, chan):
    chan = int(chan)
    if int(chan/10)==1:
        isr="ISR"
    elif int(chan/10)==0:
        isr="non ISR"
    elif int(chan/10)==2:
        isr="combined"

    if chan%3==0:
        flv="ee"
    elif chan%10==1:
        flv="e\mu"
    elif chan%10==2:
        flv="\mu\mu"

    latex = "\\begin{frame}{$\Delta m=%d$ GeV, %s $%s$ channel}\n" % (dm, isr, flv) 
    latex += "\\begin{table}\n"
    latex += "\\begin{tabular}{c c c c c c c c}\n"
    latex += "& MinNodeSize & NTrees & $KS_{sig}$ & $KS_{bkg}$ & $\chi^2_{sig}$ & $\chi^2_{bkg}$ & $\sigma_{max}$ \\ \n"
    latex += "\hline\n"
    print latex

def printRow(row):
    
beginTable(20,1)

