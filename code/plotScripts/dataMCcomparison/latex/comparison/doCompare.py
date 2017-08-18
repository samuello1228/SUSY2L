#!/usr/bin/env python
import re, os
from ROOT import *
from subprocess import call
### get the files
dirAP = re.compile('.*FigDirA{(.*)} *$')
# dirAP = re.compile('^\\def\\FigDirA\{(.*)\}$')
dirBP = re.compile('.*FigDirB{(.*)} *$')
outNameP = re.compile('%outName:(.*)$')

sDirectly = False
if gROOT.IsBatch(): sDirectly = True


def checkDir(x):
    plotdir='~/links/eosOther/cloShared/save'
#     mkdir1 = 'if [ ! -d {0:s} ]; then mkdir {0:s}; cd {0:s}; ln -s {1:s}/plot*/*.eps .; ln -s {1:s}/*.tex .; cd ..; fi'.format(x.replace('.','p'),plotdir+'/'+x[6:])
#     mkdir1 = 'if [ ! -d {0:s} ]; then mkdir {0:s}; cd {0:s}; ln -s {1:s}/plot*/*.eps .; ln -s {1:s}/*.tex .; cd ..; fi'.format(x,plotdir+'/'+x[6:].replace('2p8','2.8'))
    mkdir1 = 'if [ ! -d {0:s} ]; then mkdir {0:s}; cd {0:s}; ln -s {1:s}/plot*/*.eps .; ln -s {1:s}/*.tex .; cd ..; fi'.format(x,plotdir+'/'+x[6:])
#     print mkdir1
    call(mkdir1, shell=True)

def getComparePlots(A1, A2):
    from rootUtil import useAtlasStyle
    useAtlasStyle()

    from  compareOpt import ComparerX

    cx1 = ComparerX('tes1')
    cx1.sDir = './dCompare/'
    cx1.sTag = ('plots_'+A1+'_VS_plots_'+A2+'_').replace('.','p')
    cx1.autoSave = sDirectly
    print cx1.sDir
#     return

    print A1,A2
    ### git version
#     cx1.dir0 = './'
    cx1.dir0 = '../data/optimization/'
    cx1.rx1 = ('significance_'+A1,A1)
    cx1.rx2 = ('significance_'+A2,A2)
    cx1.opt = 'opt'
    
#     cx1.compareChans()
    cx1.compare2X()

def main():
    # %outName: test1.pdf
    # \title[Muon $|\eta|<2.8$]{SS Analysis}
    # \def\FigDirA{plots_2.8_1D_avg}
    # \def\FigDirB{plots_2.8_1D_loose_avg}
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-P', "--producePlots", help="if produce the plots", action='store_true', default=False)
    parser.add_option('-b', "--batch", help="if produce the plots", action='store_true', default=False)
    parser.add_option("--oCmd", help='opening command', default='okular')
    parser.add_option('-N', "--noOpen", help="don't open pdf file when it's compiled", action='store_true', default=False)
    (options, args) = parser.parse_args()


    pts = [('dA',dirAP), ('dB',dirBP), ('oN',outNameP)]
    vx = {}

    with open('metadata.tex','read') as f1:
        lines = f1.readlines()

    started = False
    for line in lines:
        if line[:5] == "%HERE":
            started = True
            continue
        if started:
            if line[:3] == '%--':
                started = False
                break
            
            for pt in pts:
                m = pt[1].match(line)
                if m:
                    vx[pt[0]] = m.group(1)


#     getComparePlots(vx['dA'][6:].replace('2p8','2.8'), vx['dB'][6:].replace('2p8','2.8'))
    if options.producePlots:
        getComparePlots(vx['dA'][6:], vx['dB'][6:])
#     return
    ### get plot dirs and compile
    checkDir(vx['dA'])
    checkDir(vx['dB'])
    compileCmd = 'pdflatex compare1; mv compare1.pdf '+vx['oN']
    call(compileCmd, shell=True)

    autoOpen = True
    openCmd = ['okular', vx['oN']]
    if autoOpen:
#     if not options.noOpen:
#         openCmd = [options.oCmd, pdfFile]
        if which(openCmd[0]) is None:
            print 'command', openCmd[0], 'does not exist'
            return
        if call(['pgrep', '-f', ' '.join(openCmd)])==0:
            print 'file', openCmd[1], 'probably is already opened with', openCmd[0], ', find/check/refersh your windows'
            print 'run `',' '.join(openCmd),'` to debug'
            return
        if not os.path.isfile(openCmd[1]):
            print 'file', openCmd[1], ' is not found. Check early (file production) processes'
            return
        print 'Opening file', openCmd[1], 'with', openCmd[0]
        Popen(openCmd)

def which(cmd):
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, cmd)):
                return os.path.join(path, cmd)
    return None




if __name__ == '__main__':
    main()
