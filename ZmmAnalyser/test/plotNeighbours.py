import ROOT

ROOT.gROOT.SetBatch(1)

file0=ROOT.TFile.Open("output_Neighbours_20GeV.root","READ")

canvas= ROOT.TCanvas("c","c", 600, 600)

canvas.Draw()

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetHistLineWidth(2)
ROOT.gROOT.ForceStyle()

xtalEn=file0.Get("xtalEn")
xtalEn_vs_sum8=file0.Get("xtalEn_vs_sum8")

hSum8=file0.Get("hSum8")
hAve8=file0.Get("hAve8")
hE1overSum8=file0.Get("hE1overSum8")
hE2overSum8=file0.Get("hE2overSum8")
hE3overSum8=file0.Get("hE3overSum8")
hE4overSum8=file0.Get("hE4overSum8")
hE6overSum8=file0.Get("hE6overSum8")
hE7overSum8=file0.Get("hE7overSum8")
hE8overSum8=file0.Get("hE8overSum8")
hE9overSum8=file0.Get("hE9overSum8")

firstNeigh=file0.Get("firstNeigh")
secondNeigh=file0.Get("secondNeigh")
thirdNeigh=file0.Get("thirdNeigh")
fourthNeigh=file0.Get("fourthNeigh")

firstEnFrac=file0.Get("firstEnFrac")
secondEnFrac=file0.Get("secondEnFrac")
thirdEnFrac=file0.Get("thirdEnFrac")
fourthEnFrac=file0.Get("fourthEnFrac")
fifthEnFrac=file0.Get("fifthEnFrac")
sixthEnFrac=file0.Get("sixthEnFrac")
seventhEnFrac=file0.Get("seventhEnFrac")
eighthEnFrac=file0.Get("eighthEnFrac")

secondfirstEnAsymmetry_vs_sum8 =file0.Get("secondfirstEnAsymmetry_vs_sum8") 
thirdfirstEnAsymmetry_vs_sum8  =file0.Get("thirdfirstEnAsymmetry_vs_sum8")  
fourthfirstEnAsymmetry_vs_sum8 =file0.Get("fourthfirstEnAsymmetry_vs_sum8") 
fifthfirstEnAsymmetry_vs_sum8  =file0.Get("fifthfirstEnAsymmetry_vs_sum8")  
sixthfirstEnAsymmetry_vs_sum8  =file0.Get("sixthfirstEnAsymmetry_vs_sum8")  
seventhfirstEnAsymmetry_vs_sum8=file0.Get("seventhfirstEnAsymmetry_vs_sum8")
eighthfirstEnAsymmetry_vs_sum8 =file0.Get("eighthfirstEnAsymmetry_vs_sum8") 

secondEnFrac_vs_firstEnFrac=file0.Get("secondEnFrac_vs_firstEnFrac")
thirdEnFrac_vs_firstEnFrac =file0.Get("thirdEnFrac_vs_firstEnFrac")


nXtal90=file0.Get("nXtal90")
nXtal80=file0.Get("nXtal80")
nXtal70=file0.Get("nXtal70")
nXtal50=file0.Get("nXtal50")
nXtal30=file0.Get("nXtal30")

firstNeighIndex=file0.Get("firstNeighIndex")
secondNeighIndex=file0.Get("secondNeighIndex")

second_vs_firstNeighIndex=file0.Get("second_vs_firstNeighIndex")
first_vs_sum8 =file0.Get("first_vs_sum8")
second_vs_sum8=file0.Get("second_vs_sum8")
third_vs_sum8 =file0.Get("third_vs_sum8")

hE1oS_vs_sum8=file0.Get("hE1oS_vs_sum8")
hE2oS_vs_sum8=file0.Get("hE2oS_vs_sum8")
hE3oS_vs_sum8=file0.Get("hE3oS_vs_sum8")
hE4oS_vs_sum8=file0.Get("hE4oS_vs_sum8")
hE6oS_vs_sum8=file0.Get("hE6oS_vs_sum8")
hE7oS_vs_sum8=file0.Get("hE7oS_vs_sum8")
hE8oS_vs_sum8=file0.Get("hE8oS_vs_sum8")
hE9oS_vs_sum8=file0.Get("hE9oS_vs_sum8")

cornerOverCross=file0.Get("cornerOverCross")
cornerOverCross_vs_sum8=file0.Get("cornerOverCross_vs_sum8")
minTwoMaxTwoRatio=file0.Get("minTwoMaxTwoRatio")
minTwoMaxTwoRatio_vs_sum8=file0.Get("minTwoMaxTwoRatio_vs_sum8")
minTwoMaxTwoRatio_vs_invSum8=file0.Get("minTwoMaxTwoRatio_vs_invSum8")
maxTwoMinTwoRatio_vs_sum8=file0.Get("maxTwoMinTwoRatio_vs_sum8")
cornerOverCross_vs_invSum8=file0.Get("cornerOverCross_vs_invSum8")

firstNeighIndex_vs_cornerOverCross=file0.Get("firstNeighIndex_vs_cornerOverCross")
secondNeighIndex_vs_cornerOverCross=file0.Get("secondNeighIndex_vs_cornerOverCross")

firstNeighIndex_vs_xtalEn=file0.Get("firstNeighIndex_vs_xtalEn")
secondNeighIndex_vs_xtalEn=file0.Get("secondNeighIndex_vs_xtalEn")
maxTwoMinTwoRatio_vs_xtalEn=file0.Get("maxTwoMinTwoRatio_vs_xtalEn")

second_vs_first=file0.Get("second_vs_first")
third_vs_first=file0.Get("third_vs_first")

leftOverSum8 =file0.Get("leftOverSum8")
rightOverSum8=file0.Get("rightOverSum8")
leftRightOverSum8=file0.Get("leftRightOverSum8")
left_vs_sum8 =file0.Get("left_vs_sum8")
right_vs_sum8=file0.Get("right_vs_sum8")
leftRight_vs_sum8=file0.Get("leftRight_vs_sum8")

leftRightOverSum8_vs_sum8=file0.Get("leftRightOverSum8_vs_sum8")

upOverSum8 =file0.Get("upOverSum8")
downOverSum8=file0.Get("downOverSum8")
upDownOverSum8=file0.Get("upDownOverSum8")
up_vs_sum8 =file0.Get("up_vs_sum8")
down_vs_sum8=file0.Get("down_vs_sum8")
upDown_vs_sum8=file0.Get("upDown_vs_sum8")
upDownOverSum8_vs_sum8=file0.Get("upDownOverSum8_vs_sum8")


canvas.SetRightMargin(0.15)
xtalEn_vs_sum8.Scale(1./xtalEn_vs_sum8.Integral())
xtalEn_vs_sum8.Draw("COLZ")
xtalEn_vs_sum8.GetYaxis().SetTitle("recovered crystal energy [GeV]")
xtalEn_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
canvas.SaveAs("plots/xtalEn_vs_sum8.png")
canvas.SaveAs("plots/xtalEn_vs_sum8.pdf")

canvas.Clear()
xtalEn_vs_sum8.Draw("COLZ")
xtalEn_vs_sum8.GetYaxis().SetTitle("recovered crystal energy [GeV]")
xtalEn_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
xtalEn_vs_sum8.GetYaxis().SetRangeUser(0, 300)
xtalEn_vs_sum8.GetXaxis().SetRangeUser(0, 300)
canvas.SaveAs("plots/xtalEn_vs_sum8_zoomed.png")
canvas.SaveAs("plots/xtalEn_vs_sum8_zoomed.pdf")


canvas.Clear()
canvas.SetLogy(1)

xtalEn.Scale(1./xtalEn.Integral())
xtalEn.Draw("HIST")
xtalEn.GetXaxis().SetTitle("recovered crystal energy [GeV]")
canvas.SaveAs("plots/xtalEn.png")
canvas.SaveAs("plots/xtalEn.pdf")

canvas.Clear()
hSum8.Scale(1./hSum8.Integral())
hSum8.Draw("HIST")
hSum8.GetXaxis().SetTitle("sum of the neighbours energy [GeV]")
canvas.SaveAs("plots/sum8.png")
canvas.SaveAs("plots/sum8.pdf")

canvas.Clear()

hAve8.Scale(1./hAve8.Integral())
hAve8.Draw("HIST")
hAve8.GetXaxis().SetTitle("average energy of the neighbours [GeV]")
canvas.SaveAs("plots/ave8.png")
canvas.SaveAs("plots/ave8.pdf")

canvas.Clear()

hE1overSum8.Scale(1./hE1overSum8.Integral())
hE1overSum8.Draw("HIST")
hE1overSum8.GetXaxis().SetTitle("E1/sum8")
canvas.SaveAs("plots/logE1overSum8.png")
canvas.SaveAs("plots/logE1overSum8.pdf")

canvas.Clear()

hE2overSum8.Scale(1./hE2overSum8.Integral())
hE2overSum8.Draw("HIST")
hE2overSum8.GetXaxis().SetTitle("E2/sum8")
canvas.SaveAs("plots/logE2overSum8.png")
canvas.SaveAs("plots/logE2overSum8.pdf")

canvas.Clear()
hE3overSum8.Scale(1./hE3overSum8.Integral())
hE3overSum8.Draw("HIST")
hE3overSum8.GetXaxis().SetTitle("E3/sum8")
canvas.SaveAs("plots/logE3overSum8.png")
canvas.SaveAs("plots/logE3overSum8.pdf")

canvas.Clear()
hE4overSum8.Scale(1./hE4overSum8.Integral())
hE4overSum8.Draw("HIST")
hE4overSum8.GetXaxis().SetTitle("E4/sum8")
canvas.SaveAs("plots/logE4overSum8.png")
canvas.SaveAs("plots/logE4overSum8.pdf")

canvas.Clear()

hE6overSum8.Scale(1./hE6overSum8.Integral())
hE6overSum8.Draw("HIST")
hE6overSum8.GetXaxis().SetTitle("E6/sum8")
canvas.SaveAs("plots/logE6overSum8.png")
canvas.SaveAs("plots/logE6overSum8.pdf")

canvas.Clear()
hE7overSum8.Scale(1./hE7overSum8.Integral())
hE7overSum8.Draw("HIST")
hE7overSum8.GetXaxis().SetTitle("E7/sum8")
canvas.SaveAs("plots/logE7overSum8.png")
canvas.SaveAs("plots/logE7overSum8.pdf")

canvas.Clear()
hE8overSum8.Scale(1./hE8overSum8.Integral())
hE8overSum8.Draw("HIST")
hE8overSum8.GetXaxis().SetTitle("E8/sum8")
canvas.SaveAs("plots/logE8overSum8.png")
canvas.SaveAs("plots/logE8overSum8.pdf")

canvas.Clear()

hE9overSum8.Scale(1./hE9overSum8.Integral())
hE9overSum8.Draw("HIST")
hE9overSum8.GetXaxis().SetTitle("E9/sum8")
canvas.SaveAs("plots/logE9overSum8.png")
canvas.SaveAs("plots/logE9overSum8.pdf")

canvas.Clear()
firstNeigh.Scale(1./firstNeigh.Integral())
firstNeigh.Draw("HIST")
firstNeigh.GetXaxis().SetTitle("Most energetic neighbour crystal")
canvas.SaveAs("plots/firstNeigh.png")
canvas.SaveAs("plots/firstNeigh.pdf")

canvas.Clear()
secondNeigh.Scale(1./secondNeigh.Integral())
secondNeigh.Draw("HIST")
secondNeigh.GetXaxis().SetTitle("Second most energetic neighbour crystal")
canvas.SaveAs("plots/secondNeigh.png")
canvas.SaveAs("plots/secondNeigh.pdf")

canvas.Clear()
thirdNeigh.Scale(1./thirdNeigh.Integral())
thirdNeigh.Draw("HIST")
thirdNeigh.GetXaxis().SetTitle("Third most energetic neighbour crystal")
canvas.SaveAs("plots/thirdNeigh.png")
canvas.SaveAs("plots/thirdNeigh.pdf")

canvas.Clear()
canvas.SetLogy(0)

hE1overSum8.Draw("HIST")
hE1overSum8.GetXaxis().SetTitle("E1/sum8")
canvas.SaveAs("plots/E1overSum8.png")
canvas.SaveAs("plots/E1overSum8.pdf")

canvas.Clear()

hE2overSum8.Draw("HIST")
hE2overSum8.GetXaxis().SetTitle("E2/sum8")
canvas.SaveAs("plots/E2overSum8.png")
canvas.SaveAs("plots/E2overSum8.pdf")

canvas.Clear()

hE3overSum8.Draw("HIST")
hE3overSum8.GetXaxis().SetTitle("E3/sum8")
canvas.SaveAs("plots/E3overSum8.png")
canvas.SaveAs("plots/E3overSum8.pdf")

canvas.Clear()

hE4overSum8.Draw("HIST")
hE4overSum8.GetXaxis().SetTitle("E4/sum8")
canvas.SaveAs("plots/E4overSum8.png")
canvas.SaveAs("plots/E4overSum8.pdf")

canvas.Clear()

hE6overSum8.Draw("HIST")
hE6overSum8.GetXaxis().SetTitle("E6/sum8")
canvas.SaveAs("plots/E6overSum8.png")
canvas.SaveAs("plots/E6overSum8.pdf")

canvas.Clear()

hE7overSum8.Draw("HIST")
hE7overSum8.GetXaxis().SetTitle("E7/sum8")
canvas.SaveAs("plots/E7overSum8.png")
canvas.SaveAs("plots/E7overSum8.pdf")

canvas.Clear()

hE8overSum8.Draw("HIST")
hE8overSum8.GetXaxis().SetTitle("E8/sum8")
canvas.SaveAs("plots/E8overSum8.png")
canvas.SaveAs("plots/E8overSum8.pdf")

canvas.Clear()

hE8overSum8.Draw("HIST")
hE8overSum8.GetXaxis().SetTitle("E9/sum8")
canvas.SaveAs("plots/E9overSum8.png")
canvas.SaveAs("plots/E9overSum8.pdf")

canvas.Clear()

fourthNeigh.Scale(1./fourthNeigh.Integral())
fourthNeigh.Draw("HIST")
fourthNeigh.GetXaxis().SetTitle("Fourth most energetic neighbour crystal")
canvas.SaveAs("plots/fourthNeigh.png")
canvas.SaveAs("plots/fourthNeigh.pdf")

canvas.Clear()
nXtal90.Scale(1./nXtal90.Integral())
nXtal90.Draw("HIST")
nXtal90.GetXaxis().SetTitle("number of crystal containing 90% sum8")
canvas.SaveAs("plots/nXtal90.png")
canvas.SaveAs("plots/nXtal90.pdf")

canvas.Clear()
nXtal80.Scale(1./nXtal80.Integral())
nXtal80.Draw("HIST")
nXtal80.GetXaxis().SetTitle("number of crystal containing 80% sum8")
canvas.SaveAs("plots/nXtal80.png")
canvas.SaveAs("plots/nXtal80.pdf")

canvas.Clear()
nXtal70.Scale(1./nXtal70.Integral())
nXtal70.Draw("HIST")
nXtal70.GetXaxis().SetTitle("number of crystal containing 70% sum8")
canvas.SaveAs("plots/nXtal70.png")
canvas.SaveAs("plots/nXtal70.pdf")

canvas.Clear()
nXtal50.Scale(1./nXtal50.Integral())
nXtal50.Draw("HIST")
nXtal50.GetXaxis().SetTitle("number of crystal containing 50% sum8")
canvas.SaveAs("plots/nXtal50.png")
canvas.SaveAs("plots/nXtal50.pdf")

canvas.Clear()

nXtal30.Scale(1./nXtal30.Integral())
nXtal30.Draw("HIST")
nXtal30.GetXaxis().SetTitle("number of crystal containing 30% sum8")
canvas.SaveAs("plots/nXtal30.png")
canvas.SaveAs("plots/nXtal30.pdf")

canvas.Clear()
canvas.SetLogy(0)
firstNeighIndex.Scale(1./firstNeighIndex.Integral())
firstNeighIndex.Draw("HIST")
firstNeighIndex.GetXaxis().SetTitle("index of the most energetic neighbour xtal")
canvas.SaveAs("plots/firstNeighIndex.png")
canvas.SaveAs("plots/firstNeighIndex.pdf")
canvas.Clear()
secondNeighIndex.Scale(1./secondNeighIndex.Integral())
secondNeighIndex.Draw("HIST")
secondNeighIndex.GetXaxis().SetTitle("index of the second most energetic neighbour xtal")
canvas.SaveAs("plots/secondNeighIndex.png")
canvas.SaveAs("plots/secondNeighIndex.pdf")

canvas.Clear()
second_vs_firstNeighIndex.Scale(1./second_vs_firstNeighIndex.Integral())
second_vs_firstNeighIndex.Draw("COLZTEXT")
second_vs_firstNeighIndex.GetXaxis().SetTitle("index of the most energetic neighbour xtal")
second_vs_firstNeighIndex.GetYaxis().SetTitle("index of the second most energetic neighbour xtal")
canvas.SaveAs("plots/second_vs_firstNeighIndex.png")
canvas.SaveAs("plots/second_vs_firstNeighIndex.pdf")

canvas.Clear()
firstNeighIndex_vs_xtalEn.Scale(1./firstNeighIndex_vs_xtalEn.Integral())
firstNeighIndex_vs_xtalEn.Draw("COLZ")
firstNeighIndex_vs_xtalEn.GetXaxis().SetTitle("recovered xtal energy [GeV]")
firstNeighIndex_vs_xtalEn.GetYaxis().SetTitle("index of the most energetic crystal")
canvas.SaveAs("plots/firstNeighIndex_vs_xtalEn.png")
canvas.SaveAs("plots/firstNeighIndex_vs_xtalEn.pdf")


canvas.Clear()
secondNeighIndex_vs_xtalEn.Scale(1./secondNeighIndex_vs_xtalEn.Integral())
secondNeighIndex_vs_xtalEn.Draw("COLZ")
secondNeighIndex_vs_xtalEn.GetXaxis().SetTitle("recovered xtal energy [GeV]")
secondNeighIndex_vs_xtalEn.GetYaxis().SetTitle("index of the most energetic crystal")
canvas.SaveAs("plots/secondNeighIndex_vs_xtalEn.png")
canvas.SaveAs("plots/secondNeighIndex_vs_xtalEn.pdf")

canvas.Clear()
maxTwoMinTwoRatio_vs_xtalEn.Scale(1./maxTwoMinTwoRatio_vs_xtalEn.Integral())
maxTwoMinTwoRatio_vs_xtalEn.Draw("COLZ")
maxTwoMinTwoRatio_vs_xtalEn.GetXaxis().SetTitle("recovered xtal energy [GeV]")
maxTwoMinTwoRatio_vs_xtalEn.GetYaxis().SetTitle("two max energies / two min energies")
canvas.SaveAs("plots/maxTwoMinTwoRatio_vs_xtalEn.png")
canvas.SaveAs("plots/maxTwoMinTwoRatio_vs_xtalEn.pdf")

canvas.Clear()
second_vs_first.Scale(1./second_vs_first.Integral())
second_vs_first.Draw("COLZ")
second_vs_first.GetXaxis().SetTitle("leading neighbour xtal energy [GeV]")
second_vs_first.GetYaxis().SetTitle("subleading neighbour xtal energy [GeV]")
canvas.SaveAs("plots/second_vs_first.png")
canvas.SaveAs("plots/second_vs_first.pdf")

canvas.Clear()
third_vs_first.Scale(1./third_vs_first.Integral())
third_vs_first.Draw("COLZ")
third_vs_first.GetXaxis().SetTitle("leading crystal energy [GeV]")
third_vs_first.GetYaxis().SetTitle("third most energetic crystal energy [GeV]")
canvas.SaveAs("plots/third_vs_first.png")
canvas.SaveAs("plots/third_vs_first.pdf")

canvas.Clear()
first_vs_sum8.Scale(1./first_vs_sum8.Integral())
first_vs_sum8.Draw("COLZ")
first_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
first_vs_sum8.GetYaxis().SetTitle("leading crystal energy [GeV]")
canvas.SaveAs("plots/first_vs_sum8.png")
canvas.SaveAs("plots/first_vs_sum8.pdf")

canvas.Clear()
second_vs_sum8.Scale(1./second_vs_sum8.Integral())
second_vs_sum8.Draw("COLZ")
second_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
second_vs_sum8.GetYaxis().SetTitle("subleading crystal energy [GeV]")
canvas.SaveAs("plots/second_vs_sum8.png")
canvas.SaveAs("plots/second_vs_sum8.pdf")

canvas.Clear()
third_vs_sum8.Scale(1./third_vs_sum8.Integral())
third_vs_sum8.Draw("COLZ")
third_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
third_vs_sum8.GetYaxis().SetTitle("third most energetic crystal energy [GeV]")
canvas.SaveAs("plots/third_vs_sum8.png")
canvas.SaveAs("plots/third_vs_sum8.pdf")


canvas.Clear()
hE1oS_vs_sum8.Scale(1./hE1oS_vs_sum8.Integral())
hE1oS_vs_sum8.Draw("COLZ")
hE1oS_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
hE1oS_vs_sum8.GetYaxis().SetTitle("E1/sum8 [GeV]")
canvas.SaveAs("plots/hE1oS_vs_sum8.png")
canvas.SaveAs("plots/hE1oS_vs_sum8.pdf")

canvas.Clear()
hE2oS_vs_sum8.Scale(1./hE2oS_vs_sum8.Integral())
hE2oS_vs_sum8.Draw("COLZ")
hE2oS_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
hE2oS_vs_sum8.GetYaxis().SetTitle("E2/sum8 [GeV]")
canvas.SaveAs("plots/hE2oS_vs_sum8.png")
canvas.SaveAs("plots/hE2oS_vs_sum8.pdf")

canvas.Clear()
hE3oS_vs_sum8.Scale(1./hE3oS_vs_sum8.Integral())
hE3oS_vs_sum8.Draw("COLZ")
hE3oS_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
hE3oS_vs_sum8.GetYaxis().SetTitle("E3/sum8 [GeV]")
canvas.SaveAs("plots/hE3oS_vs_sum8.png")
canvas.SaveAs("plots/hE3oS_vs_sum8.pdf")

canvas.Clear()
hE4oS_vs_sum8.Scale(1./hE4oS_vs_sum8.Integral())
hE4oS_vs_sum8.Draw("COLZ")
hE4oS_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
hE4oS_vs_sum8.GetYaxis().SetTitle("E4/sum8 [GeV]")
canvas.SaveAs("plots/hE4oS_vs_sum8.png")
canvas.SaveAs("plots/hE4oS_vs_sum8.pdf")

canvas.Clear()
hE6oS_vs_sum8.Scale(1./hE6oS_vs_sum8.Integral())
hE6oS_vs_sum8.Draw("COLZ")
hE6oS_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
hE6oS_vs_sum8.GetYaxis().SetTitle("E6/sum8 [GeV]")
canvas.SaveAs("plots/hE6oS_vs_sum8.png")
canvas.SaveAs("plots/hE6oS_vs_sum8.pdf")

canvas.Clear()
hE7oS_vs_sum8.Scale(1./hE7oS_vs_sum8.Integral())
hE7oS_vs_sum8.Draw("COLZ")
hE7oS_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
hE7oS_vs_sum8.GetYaxis().SetTitle("E7/sum8 [GeV]")
canvas.SaveAs("plots/hE7oS_vs_sum8.png")
canvas.SaveAs("plots/hE7oS_vs_sum8.pdf")

canvas.Clear()
hE8oS_vs_sum8.Scale(1./hE8oS_vs_sum8.Integral())
hE8oS_vs_sum8.Draw("COLZ")
hE8oS_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
hE8oS_vs_sum8.GetYaxis().SetTitle("E8/sum8 [GeV]")
canvas.SaveAs("plots/hE8oS_vs_sum8.png")
canvas.SaveAs("plots/hE8oS_vs_sum8.pdf")

canvas.Clear()
hE9oS_vs_sum8.Scale(1./hE9oS_vs_sum8.Integral())
hE9oS_vs_sum8.Draw("COLZ")
hE9oS_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
hE9oS_vs_sum8.GetYaxis().SetTitle("E9/sum8 [GeV]")
canvas.SaveAs("plots/hE9oS_vs_sum8.png")
canvas.SaveAs("plots/hE9oS_vs_sum8.pdf")


canvas.Clear()
canvas.SetLogy(1)
cornerOverCross.Scale(1./cornerOverCross.Integral())
cornerOverCross.Draw("HIST")
cornerOverCross.GetXaxis().SetTitle("corner energy /cross energy")
canvas.SaveAs("plots/cornerOverCross.png")
canvas.SaveAs("plots/cornerOverCross.pdf")
canvas.SetLogy(0)
canvas.Clear()
cornerOverCross_vs_sum8.Scale(1./cornerOverCross_vs_sum8.Integral())
cornerOverCross_vs_sum8.Draw("COLZ")
cornerOverCross_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
cornerOverCross_vs_sum8.GetYaxis().SetTitle("corner energy /cross energy")
canvas.SaveAs("plots/cornerOverCross_vs_sum8.png")
canvas.SaveAs("plots/cornerOverCross_vs_sum8.pdf")

canvas.Clear()
minTwoMaxTwoRatio.Scale(1./minTwoMaxTwoRatio.Integral())
minTwoMaxTwoRatio.Draw("HIST")
minTwoMaxTwoRatio.GetXaxis().SetTitle("two min energies / two max energies")
canvas.SaveAs("plots/minTwoMaxTwoRatio.png")
canvas.SaveAs("plots/minTwoMaxTwoRatio.pdf")

canvas.Clear()
minTwoMaxTwoRatio_vs_sum8.Scale(1./minTwoMaxTwoRatio_vs_sum8.Integral())
minTwoMaxTwoRatio_vs_sum8.Draw("COLZ")
minTwoMaxTwoRatio_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
minTwoMaxTwoRatio_vs_sum8.GetYaxis().SetTitle("two min energies / two max energies")
canvas.SaveAs("plots/minTwoMaxTwoRatio_vs_sum8.png")
canvas.SaveAs("plots/minTwoMaxTwoRatio_vs_sum8.pdf")

canvas.Clear()
minTwoMaxTwoRatio_vs_invSum8.Scale(1./minTwoMaxTwoRatio_vs_invSum8.Integral())
minTwoMaxTwoRatio_vs_invSum8.Draw("COLZ")
minTwoMaxTwoRatio_vs_invSum8.GetXaxis().SetTitle("1/sum8 [1/GeV]")
minTwoMaxTwoRatio_vs_invSum8.GetYaxis().SetTitle("two min energies / two max energies")
canvas.SaveAs("plots/minTwoMaxTwoRatio_vs_invSum8.png")
canvas.SaveAs("plots/minTwoMaxTwoRatio_vs_invSum8.pdf")

canvas.Clear()
maxTwoMinTwoRatio_vs_sum8.Scale(1./maxTwoMinTwoRatio_vs_sum8.Integral())
maxTwoMinTwoRatio_vs_sum8.Draw("COLZ")
maxTwoMinTwoRatio_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
maxTwoMinTwoRatio_vs_sum8.GetYaxis().SetTitle("two max energies / two min energies")
canvas.SaveAs("plots/maxTwoMinTwoRatio_vs_sum8.png")
canvas.SaveAs("plots/maxTwoMinTwoRatio_vs_sum8.pdf")

canvas.Clear()
cornerOverCross_vs_invSum8.Scale(1./cornerOverCross_vs_invSum8.Integral())
cornerOverCross_vs_invSum8.Draw("COLZ")
cornerOverCross_vs_invSum8.GetXaxis().SetTitle("1/sum8 [1/GeV]")
cornerOverCross_vs_invSum8.GetYaxis().SetTitle("corner energy /cross energy")
canvas.SaveAs("plots/cornerOverCross_vs_invSum8.png")
canvas.SaveAs("plots/cornerOverCross_vs_invSum8.pdf")

canvas.Clear()
firstNeighIndex_vs_cornerOverCross.Scale(1./firstNeighIndex_vs_cornerOverCross.Integral())
firstNeighIndex_vs_cornerOverCross.Draw("COLZ")
firstNeighIndex_vs_cornerOverCross.GetXaxis().SetTitle("corner energy /cross energy")
firstNeighIndex_vs_cornerOverCross.GetYaxis().SetTitle("leading crystal index")
canvas.SaveAs("plots/firstNeighIndex_vs_cornerOverCross.png")
canvas.SaveAs("plots/firstNeighIndex_vs_cornerOverCross.pdf")

canvas.Clear()
secondNeighIndex_vs_cornerOverCross.Scale(1./secondNeighIndex_vs_cornerOverCross.Integral())
secondNeighIndex_vs_cornerOverCross.Draw("COLZ")
secondNeighIndex_vs_cornerOverCross.GetXaxis().SetTitle("corner energy /cross energy")
secondNeighIndex_vs_cornerOverCross.GetYaxis().SetTitle("subleading crystal index")
canvas.SaveAs("plots/secondNeighIndex_vs_cornerOverCross.png")
canvas.SaveAs("plots/secondNeighIndex_vs_cornerOverCross.pdf")

canvas.Clear()
leftOverSum8.Scale(1./leftOverSum8.Integral()) 
leftOverSum8.Draw("HIST")
leftOverSum8.GetXaxis().SetTitle("(E1+E2+E3)/sum8")
canvas.SaveAs("plots/leftOverSum8.png")
canvas.SaveAs("plots/leftOverSum8.pdf")

canvas.Clear()
rightOverSum8.Scale(1./rightOverSum8.Integral())
rightOverSum8.Draw("HIST")
rightOverSum8.GetXaxis().SetTitle("(E7+E8+E9)/sum8")
canvas.SaveAs("plots/rightOverSum8.png")
canvas.SaveAs("plots/rightOverSum8.pdf")

canvas.Clear()
leftRightOverSum8.Scale(1./leftRightOverSum8.Integral()) 
leftRightOverSum8.Draw("HIST")
leftRightOverSum8.GetXaxis().SetTitle("(E1+E2+E3+E7+E8+E9)/sum8")
canvas.SaveAs("plots/leftRightOverSum8.png")
canvas.SaveAs("plots/leftRightOverSum8.pdf")

canvas.Clear()
left_vs_sum8.Scale(1./left_vs_sum8.Integral())
left_vs_sum8.Draw("COLZ")
left_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
left_vs_sum8.GetYaxis().SetTitle("E1+E2+E3 [GeV]")
canvas.SaveAs("plots/left_vs_sum8.png")
canvas.SaveAs("plots/left_vs_sum8.pdf")

canvas.Clear()
right_vs_sum8.Scale(1./right_vs_sum8.Integral())
right_vs_sum8.Draw("COLZ")
right_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
right_vs_sum8.GetYaxis().SetTitle("E7+E8+E9 [GeV]")
canvas.SaveAs("plots/right_vs_sum8.png")
canvas.SaveAs("plots/right_vs_sum8.pdf")

canvas.Clear()
leftRight_vs_sum8.Scale(1./left_vs_sum8.Integral())
leftRight_vs_sum8.Draw("COLZ")
leftRight_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
leftRight_vs_sum8.GetYaxis().SetTitle("E1+E2+E3+E7+E8+E9 [GeV]")
canvas.SaveAs("plots/leftRight_vs_sum8.png")
canvas.SaveAs("plots/leftRight_vs_sum8.pdf")

canvas.Clear()
leftRightOverSum8_vs_sum8.Scale(1./leftRightOverSum8_vs_sum8.Integral())
leftRightOverSum8_vs_sum8.Draw("COLZ")
leftRightOverSum8_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
leftRightOverSum8_vs_sum8.GetYaxis().SetTitle("(E1+E2+E3+E7+E8+E9)/sum8")
canvas.SaveAs("plots/leftRightOverSum8_vs_sum8.png")
canvas.SaveAs("plots/leftRightOverSum8_vs_sum8.pdf")


canvas.Clear()
canvas.SetLogy(1)

firstEnFrac.Scale(1./firstEnFrac.Integral())
secondEnFrac.Scale(1./secondEnFrac.Integral())
thirdEnFrac.Scale(1./thirdEnFrac.Integral())
fourthEnFrac.Scale(1./fourthEnFrac.Integral())
fifthEnFrac.Scale(1./fifthEnFrac.Integral())
sixthEnFrac.Scale(1./sixthEnFrac.Integral())
seventhEnFrac.Scale(1./seventhEnFrac.Integral())
eighthEnFrac.Scale(1./eighthEnFrac.Integral())

firstEnFrac.Draw("HIST")
firstEnFrac.SetLineColor(1)
secondEnFrac.Draw("HISTSAME")
secondEnFrac.SetLineColor(2)
thirdEnFrac.Draw("HISTSAME")
thirdEnFrac.SetLineColor(3)
fourthEnFrac.Draw("HISTSAME")
fourthEnFrac.SetLineColor(4)
fifthEnFrac.Draw("HISTSAME")
fifthEnFrac.SetLineColor(5)
sixthEnFrac.Draw("HISTSAME")
sixthEnFrac.SetLineColor(6)
seventhEnFrac.Draw("HISTSAME")
seventhEnFrac.SetLineColor(7)
eighthEnFrac.Draw("HISTSAME")
eighthEnFrac.SetLineColor(9)

mymax = max([firstEnFrac.GetBinContent(firstEnFrac.GetMaximumBin()), secondEnFrac.GetBinContent(secondEnFrac.GetMaximumBin()),
             thirdEnFrac.GetBinContent(thirdEnFrac.GetMaximumBin()), fourthEnFrac.GetBinContent(fourthEnFrac.GetMaximumBin()),
             fifthEnFrac.GetBinContent(fifthEnFrac.GetMaximumBin()), sixthEnFrac.GetBinContent(sixthEnFrac.GetMaximumBin()),
             seventhEnFrac.GetBinContent(seventhEnFrac.GetMaximumBin()), eighthEnFrac.GetBinContent(eighthEnFrac.GetMaximumBin())
             ])
firstEnFrac.GetYaxis().SetRangeUser(0.001, 1.2*mymax)
firstEnFrac.GetXaxis().SetTitle("Fraction of sum8")


leg = ROOT.TLegend(0.55, 0.85, 0.8, 0.5)
leg.AddEntry(firstEnFrac, "most energetic", "l")
leg.AddEntry(secondEnFrac, "2nd", "l")
leg.AddEntry(thirdEnFrac, "3rd", "l")
leg.AddEntry(fourthEnFrac, "4th", "l")
leg.AddEntry(fifthEnFrac, "5th", "l")
leg.AddEntry(sixthEnFrac, "6th", "l")
leg.AddEntry(seventhEnFrac, "7th", "l")
leg.AddEntry(eighthEnFrac, "8th", "l")

leg.Draw()

canvas.SaveAs("plots/fractionOfSum8.png")
canvas.SaveAs("plots/fractionOfSum8.pdf")


canvas.SetLogy(0)
canvas.Clear()
secondfirstEnAsymmetry_vs_sum8.Scale(1./secondfirstEnAsymmetry_vs_sum8.Integral())
secondfirstEnAsymmetry_vs_sum8.Draw("COLZ")
secondfirstEnAsymmetry_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
secondfirstEnAsymmetry_vs_sum8.GetYaxis().SetTitle("2nd-1st Asymmetry")
text=ROOT.TText(150, 0.3, "2nd - 1st Asymmetry vs sum8")
text.SetTextFont(42)
text.Draw()
canvas.SaveAs("plots/secondfirstEnAsymmetry_vs_sum8.png")
canvas.SaveAs("plots/secondfirstEnAsymmetry_vs_sum8.pdf")

canvas.Clear()
thirdfirstEnAsymmetry_vs_sum8.Scale(1./thirdfirstEnAsymmetry_vs_sum8.Integral())
thirdfirstEnAsymmetry_vs_sum8.Draw("COLZ")
thirdfirstEnAsymmetry_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
thirdfirstEnAsymmetry_vs_sum8.GetYaxis().SetTitle("3rd-1st Asymmetry")
text.Clear()
text.DrawText(150, 0.3, "3rd - 1st Asymmetry vs sum8")
canvas.SaveAs("plots/thirdfirstEnAsymmetry_vs_sum8.png")
canvas.SaveAs("plots/thirdfirstEnAsymmetry_vs_sum8.pdf")


canvas.Clear()
fourthfirstEnAsymmetry_vs_sum8.Scale(1./fourthfirstEnAsymmetry_vs_sum8.Integral())
fourthfirstEnAsymmetry_vs_sum8.Draw("COLZ")
fourthfirstEnAsymmetry_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
fourthfirstEnAsymmetry_vs_sum8.GetYaxis().SetTitle("4th-1st Asymmetry")
text.Clear()
text.DrawText(150, 0.3, "4th - 1st Asymmetry vs sum8")
canvas.SaveAs("plots/fourthfirstEnAsymmetry_vs_sum8.png")
canvas.SaveAs("plots/fourthfirstEnAsymmetry_vs_sum8.pdf")

canvas.Clear()
fifthfirstEnAsymmetry_vs_sum8.Scale(1./fifthfirstEnAsymmetry_vs_sum8.Integral())
fifthfirstEnAsymmetry_vs_sum8.Draw("COLZ")
fifthfirstEnAsymmetry_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
fifthfirstEnAsymmetry_vs_sum8.GetYaxis().SetTitle("5th-1st Asymmetry")
text.Clear()
text.DrawText(150, 0.3, "5th - 1st Asymmetry vs sum8")
canvas.SaveAs("plots/fifthfirstEnAsymmetry_vs_sum8.png")
canvas.SaveAs("plots/fifthfirstEnAsymmetry_vs_sum8.pdf")

canvas.Clear()
sixthfirstEnAsymmetry_vs_sum8.Scale(1./sixthfirstEnAsymmetry_vs_sum8.Integral())
sixthfirstEnAsymmetry_vs_sum8.Draw("COLZ")
sixthfirstEnAsymmetry_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
sixthfirstEnAsymmetry_vs_sum8.GetYaxis().SetTitle("6th-1st Asymmetry")
text.Clear()
text.DrawText(150, 0.3, "6th - 1st Asymmetry vs sum8")
canvas.SaveAs("plots/sixthfirstEnAsymmetry_vs_sum8.png")
canvas.SaveAs("plots/sixthfirstEnAsymmetry_vs_sum8.pdf")

canvas.Clear()
seventhfirstEnAsymmetry_vs_sum8.Scale(1./seventhfirstEnAsymmetry_vs_sum8.Integral())
seventhfirstEnAsymmetry_vs_sum8.Draw("COLZ")
seventhfirstEnAsymmetry_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
seventhfirstEnAsymmetry_vs_sum8.GetYaxis().SetTitle("7th-1st Asymmetry")
text.Clear()
text.DrawText(150, 0.3, "7th - 1st Asymmetry vs sum8")
canvas.SaveAs("plots/seventhfirstEnAsymmetry_vs_sum8.png")
canvas.SaveAs("plots/seventhfirstEnAsymmetry_vs_sum8.pdf")

canvas.Clear()
eighthfirstEnAsymmetry_vs_sum8.Scale(1./eighthfirstEnAsymmetry_vs_sum8.Integral())
eighthfirstEnAsymmetry_vs_sum8.Draw("COLZ")
eighthfirstEnAsymmetry_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
eighthfirstEnAsymmetry_vs_sum8.GetYaxis().SetTitle("8th-1st Asymmetry")
text.Clear()
text.DrawText(150, 0.3, "8th - 1st Asymmetry vs sum8")
canvas.SaveAs("plots/eighthfirstEnAsymmetry_vs_sum8.png")
canvas.SaveAs("plots/eighthfirstEnAsymmetry_vs_sum8.pdf")

canvas.Clear()
secondEnFrac_vs_firstEnFrac.Scale(1./secondEnFrac_vs_firstEnFrac.Integral())
secondEnFrac_vs_firstEnFrac.Draw("COLZ")
secondEnFrac_vs_firstEnFrac.GetXaxis().SetTitle("1st xtal En fraction")
secondEnFrac_vs_firstEnFrac.GetYaxis().SetTitle("2nd xtal En fraction")
text.Clear()
text.DrawText(0.3, 0.7, "2nd vs 1st En fraction")
canvas.SaveAs("plots/secondEnFrac_vs_firstEnFrac.png")
canvas.SaveAs("plots/secondEnFrac_vs_firstEnFrac.pdf")

canvas.Clear()
thirdEnFrac_vs_firstEnFrac.Scale(1./thirdEnFrac_vs_firstEnFrac.Integral())
thirdEnFrac_vs_firstEnFrac.Draw("COLZ")
thirdEnFrac_vs_firstEnFrac.GetXaxis().SetTitle("1st xtal En fraction")
thirdEnFrac_vs_firstEnFrac.GetYaxis().SetTitle("3rd xtal En fraction")
text.Clear()
text.DrawText(0.3, 0.7, "3rd vs 1st En fraction")
canvas.SaveAs("plots/thirdEnFrac_vs_firstEnFrac.png")
canvas.SaveAs("plots/thirdEnFrac_vs_firstEnFrac.pdf")

canvas.Clear()
upOverSum8.Scale(1./upOverSum8.Integral()) 
upOverSum8.Draw("HIST")
upOverSum8.GetXaxis().SetTitle("(E1+E2+E3)/sum8")
canvas.SaveAs("plots/upOverSum8.png")
canvas.SaveAs("plots/upOverSum8.pdf")

canvas.Clear()
downOverSum8.Scale(1./downOverSum8.Integral())
downOverSum8.Draw("HIST")
downOverSum8.GetXaxis().SetTitle("(E7+E8+E9)/sum8")
canvas.SaveAs("plots/downOverSum8.png")
canvas.SaveAs("plots/downOverSum8.pdf")

canvas.Clear()
upDownOverSum8.Scale(1./upDownOverSum8.Integral()) 
upDownOverSum8.Draw("HIST")
upDownOverSum8.GetXaxis().SetTitle("(E1+E2+E3+E7+E8+E9)/sum8")
canvas.SaveAs("plots/upDownOverSum8.png")
canvas.SaveAs("plots/upDownOverSum8.pdf")

canvas.Clear()
up_vs_sum8.Scale(1./up_vs_sum8.Integral())
up_vs_sum8.Draw("COLZ")
up_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
up_vs_sum8.GetYaxis().SetTitle("E1+E2+E3 [GeV]")
canvas.SaveAs("plots/up_vs_sum8.png")
canvas.SaveAs("plots/up_vs_sum8.pdf")

canvas.Clear()
down_vs_sum8.Scale(1./down_vs_sum8.Integral())
down_vs_sum8.Draw("COLZ")
down_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
down_vs_sum8.GetYaxis().SetTitle("E7+E8+E9 [GeV]")
canvas.SaveAs("plots/down_vs_sum8.png")
canvas.SaveAs("plots/down_vs_sum8.pdf")

canvas.Clear()
upDown_vs_sum8.Scale(1./up_vs_sum8.Integral())
upDown_vs_sum8.Draw("COLZ")
upDown_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
upDown_vs_sum8.GetYaxis().SetTitle("E1+E2+E3+E7+E8+E9 [GeV]")
canvas.SaveAs("plots/upDown_vs_sum8.png")
canvas.SaveAs("plots/upDown_vs_sum8.pdf")

canvas.Clear()
upDownOverSum8_vs_sum8.Scale(1./upDownOverSum8_vs_sum8.Integral())
upDownOverSum8_vs_sum8.Draw("COLZ")
upDownOverSum8_vs_sum8.GetXaxis().SetTitle("sum8 [GeV]")
upDownOverSum8_vs_sum8.GetYaxis().SetTitle("(E1+E2+E3+E7+E8+E9)/sum8")
canvas.SaveAs("plots/upDownOverSum8_vs_sum8.png")
canvas.SaveAs("plots/upDownOverSum8_vs_sum8.pdf")



file0.Close()


