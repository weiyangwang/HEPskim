from ROOT import TCanvas, TPad, TF1, TFile, TPaveLabel, TPaveText, TH1F, TH1I, TLorentzVector
from ROOT import gROOT, gDirectory, gBenchmark

import time, math
start = time.clock()
#ki = int(input('Enter a number: '))
#string1 = './mc'+`ki`+'_trees.root'
#string1 = 'trees_'+`ki`+'.root'
string1 = 'mc_skimTree.root'
#string2 = 'CFF_outputFileMC'+`ki`+'.root'
string2 = 'mc_histograms.root'
myfile1 = TFile(string1)
mychain1 = gDirectory.Get( 'mytree/tree' )
entries1 = mychain1.GetEntries()
print entries1

hm0vtxx      = TH1F( 'hm0vtxx', 'Vtx x position', 100, -1, 1 )
hm0vtxy      = TH1F( 'hm0vtxy', 'Vtx y position', 100, -1, 1 )
hm0vtxz      = TH1F( 'hm0vtxz', 'Vtx z position', 100, -10, 10 )

hm0pt        = TH1F( 'hm0pt', 'pT of charged tracks', 100, 0, 20 )
hm0eta       = TH1F( 'hm0eta', '#eta of charged tracks', 100, -2.2, 2.2 )
#hm0rap       = TH1F( 'hm0rap', 'Rapidity of charged trackes', 200, -20, 20)
hm0phi       = TH1F( 'hm0phi', '#phi of charged tracks', 50, -3.2, 3.2 )
hm0moch      = TH1F( 'hm0moch', 'number of charged particles', 400, 0, 400 )
hm0dz        = TH1F( 'hm0dz', 'dz', 200, -50, 50)
hm0d0        = TH1F( 'hm0d0', 'd0', 200, -40, 40)
hm0chi2n     = TH1F( 'hm0chi2n', '#chi^{2}', 200, 0, 60)
hm0highPurity= TH1F( 'hm0highPurity', 'highPurity', 5, 0, 5)
hm0algo      = TH1F( 'hm0algo', 'algo', 20, 0, 20)
hm0nValidHits= TH1F( 'hm0nValidHits', 'Valid hits', 60, 0, 60)
hm0nLostHits = TH1F( 'hm0nLostHits', 'Lost hits', 20, 0, 20)
hm0pterr     = TH1F( 'hm0pterr' , 'ptError', 100, 0, 5)
hm0d0err     = TH1F( 'hm0d0err' , 'd0Error', 100, -5, 5)
hm0dzerr     = TH1F( 'hm0dzerr' , 'dzError',  100, -5, 5)
hm0ptun      = TH1F( 'hm0ptun' , '#sigma_{pT}/pT', 100, 0, 1)
hm0d0byerr   = TH1F( 'hm0d0byerr' ,  'd0/#sigma_{d0}', 100, -5, 5)
hm0dzbyerr   = TH1F( 'hm0dzbyerr' , 'dz/#sigma_{dz}', 100, -5, 5)
hm0dzC       = TH1F( 'hm0dzC', 'dz (calculated)', 200, -50, 50)
hm0d0C       = TH1F( 'hm0d0C', 'd0 (calculated)', 200, -40, 40)

hm1pt        = TH1F( 'hm1pt', 'pT of charged tracks', 100, 0, 20 )
hm1eta       = TH1F( 'hm1eta', '#eta of charged tracks', 100, -2.2, 2.2 )
#hm1rap       = TH1F( 'hm1rap', 'Rapidity of charged trackes', 200, -20, 20)
hm1phi       = TH1F( 'hm1phi', '#phi of charged tracks', 50, -3.2, 3.2 )
hm1moch      = TH1F( 'hm1moch', 'number of charged particles', 200, 0, 200 )
hm1dz        = TH1F( 'hm1dz', 'dz', 200, -50, 50)
hm1d0        = TH1F( 'hm1d0', 'd0', 200, -40, 40)
hm1chi2n     = TH1F( 'hm1chi2n', '#chi^{2}', 200, 0, 60)
hm1highPurity= TH1F( 'hm1highPurity', 'highPurity', 5, 0, 5)
hm1algo      = TH1F( 'hm1algo', 'algo', 20, 0, 20)
hm1nValidHits= TH1F( 'hm1nValidHits', 'Valid hits', 60, 0, 60)
hm1nLostHits = TH1F( 'hm1nLostHits', 'Lost hits', 20, 0, 20)
hm1pterr     = TH1F( 'hm1pterr' , 'ptError', 100, 0, 1)
hm1d0err     = TH1F( 'hm1d0err' , 'd0Error', 100, -5, 5)
hm1dzerr     = TH1F( 'hm1dzerr' , 'dzError',  100, -5, 5)
hm1ptun      = TH1F( 'hm1ptun' , '#sigma_{pT}/pT', 100, 0, 1)
hm1d0byerr   = TH1F( 'hm1d0byerr' ,  'd0/#sigma_{d0}', 100, -5, 5)
hm1dzbyerr   = TH1F( 'hm1dzbyerr' , 'dz/#sigma_{dz}', 100, -5, 5)
hm1dzC       = TH1F( 'hm1dzC', 'dz (calculated)', 200, -50, 50)
hm1d0C       = TH1F( 'hm1d0C', 'd0 (calculated)', 200, -40, 40)

hm2pt        = TH1F( 'hm2pt', 'pT of charged tracks', 100, 0, 20 )
hm2eta       = TH1F( 'hm2eta', '#eta of charged tracks', 100, -2.2, 2.2 )
#hm2rap       = TH1F( 'hm2rap', 'Rapidity of charged trackes', 200, -20, 20)
hm2phi       = TH1F( 'hm2phi', '#phi of charged tracks', 50, -3.2, 3.2 )
hm2moch      = TH1F( 'hm2moch', 'number of charged particles', 20, 0, 20 )
hm2dz        = TH1F( 'hm2dz', 'dz', 200, -50, 50)
hm2d0        = TH1F( 'hm2d0', 'd0', 200, -40, 40)
hm2chi2n     = TH1F( 'hm2chi2n', '#chi^{2}', 200, 0, 60)
hm2highPurity= TH1F( 'hm2highPurity', 'highPurity', 5, 0, 5)
hm2algo      = TH1F( 'hm2algo', 'algo', 20, 0, 20)
hm2nValidHits= TH1F( 'hm2nValidHits', 'Valid hits', 60, 0, 60)
hm2nLostHits = TH1F( 'hm2nLostHits', 'Lost hits', 20, 0, 20)
hm2pterr     = TH1F( 'hm2pterr' , 'ptError', 100, 0, 1)
hm2d0err     = TH1F( 'hm2d0err' , 'd0Error', 100, -5, 5)
hm2dzerr     = TH1F( 'hm2dzerr' , 'dzError',  100, -5, 5)
hm2ptun      = TH1F( 'hm2ptun' , '#sigma_{pT}/pT', 100, 0, 1)
hm2d0byerr   = TH1F( 'hm2d0byerr' ,  'd0/#sigma_{d0}', 100, -5, 5)
hm2dzbyerr   = TH1F( 'hm2dzbyerr' , 'dz/#sigma_{dz}', 100, -5, 5)
hm2dzC       = TH1F( 'hm2dzC', 'dz (calculated)', 200, -50, 50)
hm2d0C       = TH1F( 'hm2d0C', 'd0 (calculated)', 200, -40, 40)


for jentry in xrange( entries1 ):
    ientry = mychain1.LoadTree( jentry )
    if ientry < 0:
        break

    nb = mychain1.GetEntry( jentry )
    if nb <= 0:
        continue

    nEvent = int(mychain1.trPt.size())
    if nEvent <= 0:
        continue
    
    if mychain1.Lumi < 90: continue

    hm0moch.Fill(mychain1.trPt.size())
    moch0 = 0
    moch1 = 0
    moch2 = 0
    numgoodvtx = 0
    min = -1
    if mychain1.vtxisFake.size() < mychain1.vtxisFakeBS.size(): min = mychain1.vtxisFake.size()
    else : min = mychain1.vtxisFakeBS.size()
    for i1 in xrange(0, min):
        vtxrho = math.sqrt(mychain1.vtxx.at(i1)*mychain1.vtxx.at(i1) + mychain1.vtxy.at(i1)*mychain1.vtxy.at(i1)) # count only good primary vertices
        if not mychain1.vtxisFake.at(i1) and abs(mychain1.vtxz.at(i1) - mychain1.vtxzBS.at(i1)) <= 10 and mychain1.vtxndof.at(i1) > 4 and vtxrho <= 2:
            numgoodvtx+=1
            if numgoodvtx==1:
                vtx_x = mychain1.vtxx.at(i1)
                vtx_y = mychain1.vtxy.at(i1)
                vtx_z = mychain1.vtxz.at(i1)
                hm0vtxx.Fill( vtx_x )
                hm0vtxy.Fill( vtx_y )
                hm0vtxz.Fill( vtx_z )

            #print vtxrho, mychain1.vtxisFake.at(i1), mychain1.vtxndof.at(i1), "sssssssss"
    #if numgoodvtx != 1: continue # and numgoodvtx != 1: print numgoodvtx
       
    if(mychain1.trPt.size() != mychain1.d0err.size()): print "ssssssssss"
    #if numgoodvtx > 1: print "aaaaaaaaaaaaaaaaaaaaaa"
    if numgoodvtx < 1: continue

    for i in xrange(0, mychain1.trPt.size()):
        pt=mychain1.trPt.at(i)
        eta = mychain1.trEta.at(i)
        phi = mychain1.trPhi.at(i)
        p4 = TLorentzVector()
        p4.SetPtEtaPhiM(pt, eta, phi, 0)
        tr_d0Err = mychain1.d0err.at(i)
        tr_dzErr = mychain1.dzerr.at(i)
        tr_ptErr = mychain1.ptErr.at(i)
        tr_x   = mychain1.vx.at(i)
        tr_y   = mychain1.vy.at(i)
        tr_z   = mychain1.vz.at(i)
        trpx = p4.Px()
        trpy = p4.Py()
        trpz = p4.Pz()
        tr_d0  = (- (tr_x-vtx_x) * trpy + (tr_y - vtx_y) * trpx)/pt
        tr_dz  = (tr_z-vtx_z) - ((tr_x-vtx_x)*trpx + (tr_y-vtx_y)*trpy)/pt * trpz/pt

        if pt <= 0.5 or abs(eta)>=2:
            continue
        
        moch0 += 1
        ptun1 = mychain1.ptErr.at(i)/pt
        d0byerr1 = mychain1.d0.at(i)/tr_d0Err
        dzbyerr1 = mychain1.dz.at(i)/tr_dzErr
        chi2n1 = mychain1.chi2n.at(i)

        hm0pt.Fill(pt)
        hm0eta.Fill(eta)
        hm0phi.Fill(phi)
        hm0dz.Fill(mychain1.dz.at(i))
        hm0d0.Fill(mychain1.d0.at(i))
        hm0chi2n.Fill(mychain1.chi2n.at(i)) 
        hm0highPurity.Fill(mychain1.highPurity.at(i))
        hm0algo.Fill(mychain1.algo.at(i))
        hm0nValidHits.Fill(mychain1.nValidHits.at(i))
        hm0nLostHits.Fill(mychain1.nLostHits.at(i))
        hm0pterr.Fill(mychain1.ptErr.at(i))
        hm0d0err.Fill(mychain1.d0.at(i))
        hm0dzerr.Fill(mychain1.dz.at(i))
        hm0ptun.Fill(ptun1)
        hm0d0byerr.Fill(d0byerr1)
        hm0dzbyerr.Fill(dzbyerr1)
        hm0dzC.Fill(tr_dz)
        hm0d0C.Fill(tr_d0)

        if abs(d0byerr1)<3 and abs(dzbyerr1)<3 and mychain1.highPurity.at(i) != 0:# and chi2n1 > 10:
            hm2ptun.Fill(ptun1)
            hm2pterr.Fill(mychain1.ptErr.at(i))
        if abs(ptun1)<0.05 and abs(d0byerr1)<3 and abs(dzbyerr1)<3:# and chi2n1 > 10:
            hm2highPurity.Fill(mychain1.highPurity.at(i))
        if abs(ptun1)<0.05 and abs(dzbyerr1)<3 and mychain1.highPurity.at(i) != 0:# and chi2n1 > 10:
            hm2d0err.Fill(mychain1.d0.at(i))
            hm2d0byerr.Fill(d0byerr1)
            hm2d0C.Fill(tr_d0)
        if abs(ptun1)<0.05 and abs(d0byerr1)<3 and mychain1.highPurity.at(i) != 0:# and chi2n1 > 10:
            hm2dzerr.Fill(mychain1.dzerr.at(i))
            hm2dzbyerr.Fill(dzbyerr1)
            hm2dzC.Fill(tr_dz)

        if abs(ptun1)<0.05 and abs(d0byerr1)<3 and abs(dzbyerr1)<3 and mychain1.highPurity.at(i) != 0:# and chi2n1 > 10:
            #print chi2n1, i, "          ddddddddddd"
            hm1pt.Fill(pt)
            hm1eta.Fill(eta)
            hm1phi.Fill(phi)
            hm1dz.Fill(mychain1.dz.at(i))
            hm1d0.Fill(mychain1.d0.at(i))
            hm1chi2n.Fill(mychain1.chi2n.at(i))
            hm1highPurity.Fill(mychain1.highPurity.at(i))
            hm1algo.Fill(mychain1.algo.at(i))
            hm1nValidHits.Fill(mychain1.nValidHits.at(i))
            hm1nLostHits.Fill(mychain1.nLostHits.at(i))
            hm1pterr.Fill(mychain1.ptErr.at(i))
            hm1d0err.Fill(mychain1.d0err.at(i))
            hm1dzerr.Fill(mychain1.dzerr.at(i))
            hm1ptun.Fill(ptun1)
            hm1d0byerr.Fill(d0byerr1)
            hm1dzbyerr.Fill(dzbyerr1)
            hm1dzC.Fill(tr_dz)
            hm1d0C.Fill(tr_d0)
            moch1 += 1

    if jentry%1000==1:
        print jentry,"events got processed ..."

    hm1moch.Fill(moch0)
    hm2moch.Fill(moch1)

hm1pt.SetTitle( 'pT of charged tracks' )
hm1pt.GetXaxis().SetTitle( 'pT (in GeV)' )
hm1eta.GetXaxis().SetTitle( '#eta' )
hm1phi.GetXaxis().SetTitle( '#phi' )
hm1moch.GetXaxis().SetTitle( 'Number of charged tracks' )
#h1mopho.GetXaxis().SetTitle( 'Number of photons' )
hm1dz.GetXaxis().SetTitle( 'dz (in cm)' )
hm1d0.GetXaxis().SetTitle( 'd0 (in cm)' )
hm1chi2n.GetXaxis().SetTitle( '{#chi}^2/ndof of track fit' )
hm1highPurity.GetXaxis().SetTitle( 'Track High-purity flag' )
hm1algo.GetXaxis().SetTitle( 'Iterations' )
hm1nValidHits.GetXaxis().SetTitle( '#' )
hm1nLostHits.GetXaxis().SetTitle( '#' )
hm1pterr.GetXaxis().SetTitle( '#sigma_{pT}' )
hm1d0err.GetXaxis().SetTitle( '#sigma_{dz}' )
hm1dzerr.GetXaxis().SetTitle( '#sigma_{dz}' )
hm1ptun.GetXaxis().SetTitle( '#sigma_{pT}/pT' )
hm1d0byerr.GetXaxis().SetTitle( 'd0/#sigma_{d0}' )
hm1dzbyerr.GetXaxis().SetTitle( 'dz/#sigma_{dz}' )
hm1dzC.GetXaxis().SetTitle( 'dz (in cm)' )
hm1d0C.GetXaxis().SetTitle( 'd0 (in cm)' )
hm2ptun.GetXaxis().SetTitle( '#sigma_{pT}/pT' )
hm2pterr.GetXaxis().SetTitle( '#sigma_{pT} (in GeV)' )
hm2highPurity.GetXaxis().SetTitle( 'Track High-purity flag' )
hm2d0err.GetXaxis().SetTitle( '#sigma_{dz}' )
hm2d0byerr.GetXaxis().SetTitle( 'dz/#sigma_{dz}' )
hm2dzerr.GetXaxis().SetTitle( '#sigma_{dz}' )
hm2dzbyerr.GetXaxis().SetTitle( 'dz/#sigma_{dz}' )
hm2dzC.GetXaxis().SetTitle( 'dz (in cm)' )
hm2d0C.GetXaxis().SetTitle( 'd0 (in cm)' )
hm2moch.GetXaxis().SetTitle( 'Number of charged tracks' )

hm1pt.GetYaxis().SetTitle( '# of events' )
hm1eta.GetYaxis().SetTitle( '# of events' )
hm1phi.GetYaxis().SetTitle( '# of events' )
hm1moch.GetYaxis().SetTitle( '# of events' )
hm1dz.GetYaxis().SetTitle( '# of events' )
hm1d0.GetYaxis().SetTitle( '# of events' )
hm1chi2n.GetYaxis().SetTitle( '# of events' )
hm1highPurity.GetYaxis().SetTitle( '# of events' )
hm1algo.GetYaxis().SetTitle( '# of events' )
hm1nValidHits.GetYaxis().SetTitle( '# of events' )
hm1nLostHits.GetYaxis().SetTitle( '# of events' )
hm1pterr.GetYaxis().SetTitle( '# of events' )
hm1d0err.GetYaxis().SetTitle( '# of events' )
hm1dzerr.GetYaxis().SetTitle( '# of events' )
hm1ptun.GetYaxis().SetTitle( '# of events' )
hm1d0byerr.GetYaxis().SetTitle( '# of events' )
hm1dzbyerr.GetYaxis().SetTitle( '# of events' )
hm1dzC.GetYaxis().SetTitle( '# of events' )
hm1d0C.GetYaxis().SetTitle( '# of events' )
hm2ptun.GetYaxis().SetTitle( '# of events' )
hm2pterr.GetYaxis().SetTitle( '# of events' )
hm2highPurity.GetYaxis().SetTitle( '# of events' )
hm2d0err.GetYaxis().SetTitle( '# of events' )
hm2d0byerr.GetYaxis().SetTitle( '# of events' )
hm2dzerr.GetYaxis().SetTitle( '# of events' )
hm2dzbyerr.GetYaxis().SetTitle( '# of events' )
hm2moch.GetYaxis().SetTitle( '# of events' )

hm0pt.GetXaxis().SetTitle( 'pT of charged tracks' )
hm0eta.GetXaxis().SetTitle( '#eta' )
hm0phi.GetXaxis().SetTitle( '#phi' )
hm0moch.GetXaxis().SetTitle( 'Number of charged tracks' )
#h1mopho.GetXaxis().SetTitle( 'Number of photons' )
hm0dz.GetXaxis().SetTitle( 'dz (in cm)' )
hm0d0.GetXaxis().SetTitle( 'd0 (in cm)' )
hm0chi2n.GetXaxis().SetTitle( '{#chi}^2/ndof of track fit' )
hm0highPurity.GetXaxis().SetTitle( 'Track High-purity flag' )
hm0algo.GetXaxis().SetTitle( 'Iterations' )
hm0nValidHits.GetXaxis().SetTitle( '#' )
hm0nLostHits.GetXaxis().SetTitle( '#' )
hm0pterr.GetXaxis().SetTitle( '#sigma_{pT}' )
hm0d0err.GetXaxis().SetTitle( '#sigma_{dz}' )
hm0dzerr.GetXaxis().SetTitle( '#sigma_{dz}' )
hm0ptun.GetXaxis().SetTitle( '#sigma_{pT}/pT' )
hm0d0byerr.GetXaxis().SetTitle( 'dz/#sigma_{dz}' )
hm0dzbyerr.GetXaxis().SetTitle( 'dz/#sigma_{dz}' )
hm0vtxx.GetXaxis().SetTitle( 'vtx x (in cm)' )
hm0vtxy.GetXaxis().SetTitle( 'vtx y (in cm)' )
hm0vtxz.GetXaxis().SetTitle( 'vtx y (in cm)' )
hm0dzC.GetXaxis().SetTitle( 'dz (in cm)' )
hm0d0C.GetXaxis().SetTitle( 'd0 (in cm)' )

hm0pt.GetYaxis().SetTitle( '# of events' )
hm0eta.GetYaxis().SetTitle( '# of events' )
hm0phi.GetYaxis().SetTitle( '# of events' )
hm0moch.GetYaxis().SetTitle( '# of events' )
hm0dz.GetYaxis().SetTitle( '# of events' )
hm0d0.GetYaxis().SetTitle( '# of events' )
hm0chi2n.GetYaxis().SetTitle( '# of events' )
hm0highPurity.GetYaxis().SetTitle( '# of events' )
hm0algo.GetYaxis().SetTitle( '# of events' )
hm0nValidHits.GetYaxis().SetTitle( '# of events' )
hm0nLostHits.GetYaxis().SetTitle( '# of events' )
hm0pterr.GetYaxis().SetTitle( '# of events' )
hm0d0err.GetYaxis().SetTitle( '# of events' )
hm0dzerr.GetYaxis().SetTitle( '# of events' )
hm0ptun.GetYaxis().SetTitle( '# of events' )
hm0d0byerr.GetYaxis().SetTitle( '# of events' )
hm0dzbyerr.GetYaxis().SetTitle( '# of events' )
hm0vtxx.GetYaxis().SetTitle( '# of events' )
hm0vtxy.GetYaxis().SetTitle( '# of events' )
hm0vtxz.GetYaxis().SetTitle( '# of events' )
hm0dzC.GetYaxis().SetTitle( '# of events' )
hm0d0C.GetYaxis().SetTitle( '# of events' )

outfile1 = TFile(string2, "recreate")
hm0vtxx.Write()
hm0vtxy.Write()
hm0vtxz.Write()

hm0pt.Write()
hm0eta.Write()
hm0phi.Write()
hm0moch.Write()
hm0dz.Write()
hm0d0.Write()
hm0chi2n.Write()
hm0highPurity.Write()
hm0algo.Write()
hm0nValidHits.Write()
hm0nLostHits.Write()
hm0pterr.Write()
hm0d0err.Write()
hm0dzerr.Write()
hm0ptun.Write()
hm0d0byerr.Write()
hm0dzbyerr.Write()
hm0dzC.Write()
hm0d0C.Write()
hm1pt.Write()
hm1eta.Write()
hm1phi.Write()
hm1moch.Write()
hm1dz.Write()
hm1d0.Write()
hm1chi2n.Write()
hm1highPurity.Write()
hm1algo.Write()
hm1nValidHits.Write()
hm1nLostHits.Write()
hm1pterr.Write()
hm1d0err.Write()
hm1dzerr.Write()
hm1ptun.Write()
hm1d0byerr.Write()
hm1dzbyerr.Write()
hm1dzC.Write()
hm1d0C.Write()
hm2ptun.Write()
hm2pterr.Write()
hm2highPurity.Write()
hm2d0err.Write()
hm2d0byerr.Write()
hm2dzerr.Write()
hm2dzbyerr.Write()
hm2dzC.Write()
hm2d0C.Write()
hm2moch.Write()
outfile1.Close()

end = time.clock()
print "%.2g seconds are wasted as execution time....." %(end - start)

