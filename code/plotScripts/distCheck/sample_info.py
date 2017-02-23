#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory

funlist=[]

class sampleInfo:
    def __init__(self, filename=None, tag='default'):
        self.trees = {}
        self.tree = None
        self.Unit_fb = False
        if filename: self.loadFile(filename, tag)
    def loadAll(self):
#         dir1 = '/home/dzhang/work/bsmSearch/ewSUSY/analysis/r20_evtSelect/RootCoreBin/data/SUSYTools/mc15_13TeV/'
#         dir1 = '../../RootCoreBin/data/SUSYTools/mc15_13TeV/'
        dir1 = '/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/SUSYTools/mc15_13TeV/'
        self.loadFile(dir1+'Backgrounds.txt')
        self.loadFile(dir1+'MGPy8EG_A14N23LO_C1N2_WZ_XX_YY.txt')
        self.loadFile(dir1+'MGPy8EG_A14N23LO_C1N2_Slep_XX_YY.txt')
        self.loadFile('../../multiLepSearch/script/SigXsec.txt')

    def setDefaultTree(tag):
        self.tree = self.trees[tag]

    def loadFile(self, filename, tag='default'):
        if self.trees.has_key(tag):
            tree = self.trees[tag]
            tree.ReadFile(filename)
        else:
            tree = TTree()
            self.trees[tag] = tree
            tree.ReadFile(filename, self.getHeader(filename))

    def getHeader(self, filename):
        with open(filename) as f:
            first_line = f.readline()
            return '' if first_line.find(':')!=-1 else 'id/I:name/C:xsec/F:kfac/F:eff/F:relunc/F'
        #id/I:name/C:xsec/F:kfac/F:eff/F:relunc/F

    def getXSec(self, query, tag=None):
        xsec = 0
        tr1 = self.trees[tag] if tag else self.tree
        if not tr1: tr1 = self.trees['default']

        tr1.Draw('>>slist', query)
        el1 = gDirectory.Get('slist')
        for i in range(el1.GetN()):
            tr1.Show(el1.GetEntry(i))
            tr1.GetEntry(el1.GetEntry(i))
            print i, tr1.id, tr1.name, tr1.xsec, tr1.kfac, tr1.eff
            xsec += tr1.xsec*tr1.kfac*tr1.eff
            if self.Unit_fb:
                xsec *= 1000.
        return xsec

def test():
    sI1 = sampleInfo()
    sI1.loadAll()
#     sI1 = sampleInfo('/home/dzhang/work/bsmSearch/ewSUSY/analysis/r20_evtSelect/RootCoreBin/data/SUSYTools/mc15_13TeV/Backgrounds.txt')
#     sI1.loadFile('/home/dzhang/work/bsmSearch/ewSUSY/analysis/r20_evtSelect/RootCoreBin/data/SUSYTools/mc15_13TeV/MGPy8EG_A14N23LO_C1N2_WZ_XX_YY.txt')
    print sI1.getXSec('id==410000')
    print sI1.getXSec('id==392205')
    print sI1.getXSec('id==392205&&name=="127"')

#     waitRootCmd()
funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    for fun in funlist: print fun()
