import ROOT


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


fileRecov=ROOT.TFile.Open("zMuMu_PierreTag_sum8gt20.root", "READ")
fileTrue=ROOT.TFile.Open("zMuMu_AOD_2018B.root", "READ")

recoTree=fileRecov.Get("ntupler/selected")
trueTree=fileTrue.Get("ntupler/selected")

xtalEnergies=[] # (trueEn, recoEn)
histo=ROOT.TH2F("xtal_recovEnVsTrueEn", "xtal_recovEnVsTrueEn", 500, 0, 500, 500, 0, 500)
ratiohisto=ROOT.TH2F("ratiohisto", "E/Etrue vs Etrue", 500,0, 500, 50,  0, 5) 
hTrue=ROOT.TH2F("truePattern_recovEnVsTrueEn", "truePattern_recovEnVsTrueEn", 100, 0, 500, 100, 0, 500)
hFalse=ROOT.TH2F("falsePattern_recovEnVsTrueEn", "falsePattern_recovEnVsTrueEn", 100, 0, 500, 100, 0, 500)
count=0
oldevt=-999
for row in recoTree:
    if count%100==000: print 'number of events processed:', count
    count+=1
    for line  in trueTree:
        evtReco=str(row.runNumber)+':'+str(row.lumiBlock)+':'+str(row.eventNumber)
        evtTrue=str(line.runNumber)+':'+str(line.lumiBlock)+':'+str(line.eventNumber)
        if evtReco!=evtTrue: continue
        if evtTrue==oldevt: continue
        oldevt=evtTrue
        #print evtTrue
        for n,ixtalID in enumerate(row.xtalRawId):
            for m,jxtalID in enumerate(line.xtalRawId):
                if ixtalID!=jxtalID: continue
                if line.xtalEn[m]==0.: continue
                if row.xtalEn[n]==0.: continue
                xtalEnergies.append((line.xtalEn[m], row.xtalEn[n]))
                histo.Fill(line.xtalEn[m], row.xtalEn[n])
                ratiohisto.Fill(line.xtalEn[m], row.xtalEn[n]/line.xtalEn[m])

                ## ho tutte le info della 3x3. Posso verificare se e' il rec hit piu' energetico o no

                neigh=[row.vNeigh0[n], row.vNeigh1[n], row.vNeigh2[n], row.vNeigh3[n], row.vNeigh5[n], row.vNeigh6[n], row.vNeigh7[n], row.vNeigh8[n]]
                
                isGoodPattern=False
                if row.vNeigh3[n]==max(neigh) or row.vNeigh5[n]==max(neigh): isGoodPattern=True
                #print n,m, ixtalID, line.xtalEn[m], row.xtalEn[n], isGoodPattern, neigh
                if isGoodPattern:
                    hTrue.Fill(line.xtalEn[m], row.xtalEn[n])
                else:
                    hFalse.Fill(line.xtalEn[m], row.xtalEn[n])
        #print '\n'


print len(xtalEnergies)

canvas=ROOT.TCanvas("c", "c", 600, 600)
histo.Draw("COLZ")
histo.GetXaxis().SetRangeUser(0, 50)
histo.GetYaxis().SetRangeUser(0, 50)
histo.GetXaxis().SetTitle("xtal true energy [GeV]")
histo.GetYaxis().SetTitle("xtal recov energy [GeV]")

canvas.SaveAs("plots/true_vs_Etrue.png")
canvas.SaveAs("plots/true_vs_Etrue.pdf")

canvas.Clear()
ratiohisto.Draw("COLZ")

ratiohisto.GetXaxis().SetTitle("xtal true energy [GeV]")
ratiohisto.GetYaxis().SetTitle("xtal E/Etrue")
canvas.SaveAs("plots/EoverEtrue_vs_Etrue.png")
canvas.SaveAs("plots/EoverEtrue_vs_Etrue.pdf")

ratiohisto.GetXaxis().SetRangeUser(0, 50)
canvas.SaveAs("plots/EoverEtrue_vs_Etrue_zoommed.png")
canvas.SaveAs("plots/EoverEtrue_vs_Etrue_zoommed.pdf")



outfile=ROOT.TFile.Open("zMuMu_2018B_Comparison.root", "RECREATE")
outfile.cd()
histo.Write()
hTrue.Write()
hFalse.Write()
ratiohisto.Write()
outfile.Close()
print 'output file close'
