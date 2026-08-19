#This is the file that create the plots. The data has been collected and this file creates the plots. 



void mmanalys::plot_Raw()

{
    
   // 	Setstyle();
    	gStyle->SetOptStat(111111);
    
    
 //   gStyle->SetOptStat(111111);

    // ==========================================================
    // --- DRAW 2D RESIDUAL VS STRIP CANVASES ---
    // ==========================================================
    TCanvas *c_res_vs_stripX = new TCanvas("c_res_vs_stripX", "Residuals vs Strip Number (X Layers)", 1200, 800);
    c_res_vs_stripX->Divide(2, 2);
    
    c_res_vs_stripX->cd(1); residualstrips1X->Draw("colz");
    c_res_vs_stripX->cd(2); residualstrips2X->Draw("colz");
    c_res_vs_stripX->cd(3); residualstrips3X->Draw("colz");
    c_res_vs_stripX->cd(4); residualstrips4X->Draw("colz");

    TCanvas *c_X1_vs_X2 = new TCanvas("c_X1_vs_X2", "Hardware Diagnostic: Layer 1X vs Layer 2X", 800, 800);
    c_X1_vs_X2->cd();
    h_stripXY_X1X2->Draw("colz");
    TCanvas *c_X3_vs_X4 = new TCanvas("c_X3_vs_X4", "Hardware Diagnostic: Layer 3X vs Layer 4X", 800, 800);
    c_X3_vs_X4->cd();
    h_stripXY_X3X4->Draw("colz");

    // At the end, create canvases for the new histograms
    TCanvas *c_angleX_centroid = new TCanvas("c_angleX_centroid", "Track Angle X (centroid only)", 1200, 800);
    h_angleX_centroid->Draw();
    TCanvas *c_angleX_utpc = new TCanvas("c_angleX_utpc", "Track Angle X (uTPC only)", 1200, 800);
    h_angleX_utpc->Draw();
    TCanvas *c_angleY_centroid = new TCanvas("c_angleY_centroid", "Track Angle Y (centroid only)", 1200, 800);
    h_angleY_centroid->Draw();
    TCanvas *c_angleY_utpc = new TCanvas("c_angleY_utpc", "Track Angle Y (uTPC only)", 1200, 800);
    h_angleY_utpc->Draw();


    TCanvas *c_chi2_2D_X = new TCanvas("c_chi2_2D_X", "Chi2 vs Strip Heatmap (X Layers)", 1200, 800);
    c_chi2_2D_X->Divide(2,2);
    c_chi2_2D_X->cd(1); h_chi2_vs_strip_1X->Draw("colz");
    c_chi2_2D_X->cd(2); h_chi2_vs_strip_2X->Draw("colz");
    c_chi2_2D_X->cd(3); h_chi2_vs_strip_3X->Draw("colz");
    c_chi2_2D_X->cd(4); h_chi2_vs_strip_4X->Draw("colz");

    TCanvas *c_chi2_2D_Y = new TCanvas("c_chi2_2D_Y", "Chi2 vs Strip Heatmap (Y Layers)", 1200, 800);
    c_chi2_2D_Y->Divide(2,2);
    c_chi2_2D_Y->cd(1); h_chi2_vs_strip_1Y->Draw("colz");
    c_chi2_2D_Y->cd(2); h_chi2_vs_strip_2Y->Draw("colz");
    c_chi2_2D_Y->cd(3); h_chi2_vs_strip_3Y->Draw("colz");
    c_chi2_2D_Y->cd(4); h_chi2_vs_strip_4Y->Draw("colz");

    TCanvas *c55 = new TCanvas("c55","Single-strips-X", 150,0,1200,800);
        c55->Divide(2,2);
    c55->cd(1); h_singlestrip_1X->Draw();
    c55->cd(2); h_singlestrip_2X->Draw();
    c55->cd(3); h_singlestrip_3X->Draw();
    c55->cd(4); h_singlestrip_4X->Draw();

    TCanvas *c57 = new TCanvas("c57","Single-strips-Y", 150,0,1200,800);
        c57->Divide(2,2);
    c57->cd(1); h_singlestrip_1Y->Draw();
    c57->cd(2); h_singlestrip_2Y->Draw();
    c57->cd(3); h_singlestrip_3Y->Draw();
    c57->cd(4); h_singlestrip_4Y->Draw();
  
    TCanvas *c56 = new TCanvas("c56","Single-strips-clusters-X", 150,0,1200,800);
        c56->Divide(2,2);
    c56->cd(1); h_singlestripcluster_1X->Draw();
    c56->cd(2); h_singlestripcluster_2X->Draw();
    c56->cd(3); h_singlestripcluster_3X->Draw();
    c56->cd(4); h_singlestripcluster_4X->Draw();

    TCanvas *c58 = new TCanvas("c58","Strips-X-chi2", 150,0,1200,800);
        c58->Divide(2,2);
    c58->cd(1); h_stripschi2_1X->Draw();
    c58->cd(2); h_stripschi2_2X->Draw();
    c58->cd(3); h_stripschi2_3X->Draw();
    c58->cd(4); h_stripschi2_4X->Draw();

    TCanvas *c59 = new TCanvas("c59","Strips-Y-chi2", 150,0,1200,800);
        c59->Divide(2,2);
    c59->cd(1); h_stripschi2_1Y->Draw();
    c59->cd(2); h_stripschi2_2Y->Draw();
    c59->cd(3); h_stripschi2_3Y->Draw();
    c59->cd(4); h_stripschi2_4Y->Draw();

    TCanvas *c60 = new TCanvas("h_stripXY_X1Y1cl0", "Strip Correlation Xn vs Yn)", 1200, 900);
    c60->Divide(2, 2);
    c60->cd(1); h_stripXY_X1Y1cl0->Draw("colz");
    c60->cd(2); h_stripXY_X2Y2cl0->Draw("colz");
    c60->cd(3); h_stripXY_X3Y3cl0->Draw("colz");
    c60->cd(4); h_stripXY_X4Y4cl0->Draw("colz");
/*
    TCanvas *c61 = new TCanvas("h_stripXY_X1Y1cl1", "Strip Correlation Xn vs Yn", 1200, 900);
    c61->Divide(2, 2);
    c61->cd(1); h_stripXY_X1Y1cl1->Draw("colz");
    c61->cd(2); h_stripXY_X2Y2cl1->Draw("colz");
    c61->cd(3); h_stripXY_X3Y3cl1->Draw("colz");
    c61->cd(4); h_stripXY_X4Y4cl1->Draw("colz");

    TCanvas *c62 = new TCanvas("h_stripXY_X1Y1cl2", "Strip Correlation Xn vs Yn", 1200, 900);
    c62->Divide(2, 2);
    c62->cd(1); h_stripXY_X1Y1cl2->Draw("colz");
    c62->cd(2); h_stripXY_X2Y2cl2->Draw("colz");
    c62->cd(3); h_stripXY_X3Y3cl2->Draw("colz");
    c62->cd(4); h_stripXY_X4Y4cl2->Draw("colz");
*/

/*          -----ΠΡΕΠΕΙ ΝΑ ΤΑ ΔΩ ΜΕ ΤΟ CHI^2 ΑΛΛΑ ΜΕΤΑ-----
//_________ 0 LEADING CHARGE - STRIPS LAYERS ____
// Canvas για Leading Charge (X)
TCanvas *c0_charge_X = new TCanvas("c0_charge_X", "0 Leading Charge X Layers", 1200, 800);
c0_charge_X->Divide(2,2);
c0_charge_X->cd(1); h_lead0_charge_X1->Draw(); c0_charge_X->cd(2); h_lead0_charge_X2->Draw();
c0_charge_X->cd(3); h_lead0_charge_X3->Draw(); c0_charge_X->cd(4); h_lead0_charge_X4->Draw();

// Canvas για Leading Strips (X)
TCanvas *c0_strips_X = new TCanvas("c0_strips_X", "0 Leading Strips X Layers", 1200, 800);
c0_strips_X->Divide(2,2);
c0_strips_X->cd(1); h_lead0_strips_X1->Draw(); c0_strips_X->cd(2); h_lead0_strips_X2->Draw();
c0_strips_X->cd(3); h_lead0_strips_X3->Draw(); c0_strips_X->cd(4); h_lead0_strips_X4->Draw();

// Canvas για Leading Charge 1(X)
TCanvas *c1_charge_X = new TCanvas("c1_charge_X", "1 Leading Charge X Layers", 1200, 800);
c1_charge_X->Divide(2,2);
c1_charge_X->cd(1); h_lead1_charge_X1->Draw(); c1_charge_X->cd(2); h_lead1_charge_X2->Draw();
c1_charge_X->cd(3); h_lead1_charge_X3->Draw(); c1_charge_X->cd(4); h_lead1_charge_X4->Draw();

// Canvas για Leading Strips 1(X)
TCanvas *c1_strips_X = new TCanvas("c1_strips_X", "1 Leading Strips X Layers", 1200, 800);
c1_strips_X->Divide(2,2);
c1_strips_X->cd(1); h_lead1_strips_X1->Draw(); c1_strips_X->cd(2); h_lead1_strips_X2->Draw();
c1_strips_X->cd(3); h_lead1_strips_X3->Draw(); c1_strips_X->cd(4); h_lead1_strips_X4->Draw();

// Canvas για Leading Charge 2(X)
TCanvas *c2_charge_X = new TCanvas("c2_charge_X", "2 Leading Charge X Layers", 1200, 800);
c2_charge_X->Divide(2,2);
c2_charge_X->cd(1); h_lead2_charge_X1->Draw(); c2_charge_X->cd(2); h_lead2_charge_X2->Draw();
c2_charge_X->cd(3); h_lead2_charge_X3->Draw(); c2_charge_X->cd(4); h_lead2_charge_X4->Draw();

// Canvas για Leading Strips 2(X)
TCanvas *c2_strips_X = new TCanvas("c2_strips_X", "2 Leading Strips X Layers", 1200, 800);
c2_strips_X->Divide(2,2);
c2_strips_X->cd(1); h_lead2_strips_X1->Draw(); c2_strips_X->cd(2); h_lead2_strips_X2->Draw();
c2_strips_X->cd(3); h_lead2_strips_X3->Draw(); c2_strips_X->cd(4); h_lead2_strips_X4->Draw();

// Canvas για Leading Charge 0(Y)
TCanvas *c0_charge_Y = new TCanvas("c0_charge_Y", "0 Leading Charge Y Layers", 1200, 800);
c0_charge_Y->Divide(2,2);
c0_charge_Y->cd(1); h_lead0_charge_Y1->Draw(); c0_charge_Y->cd(2); h_lead0_charge_Y2->Draw();
c0_charge_Y->cd(3); h_lead0_charge_Y3->Draw(); c0_charge_Y->cd(4); h_lead0_charge_Y4->Draw();

// Canvas για Leading Strips 0(Y)
TCanvas *c0_strips_Y = new TCanvas("c0_strips_Y", "0 Leading Strips Y Layers", 1200, 800);
c0_strips_Y->Divide(2,2);
c0_strips_Y->cd(1); h_lead0_strips_Y1->Draw(); c0_strips_Y->cd(2); h_lead0_strips_Y2->Draw();
c0_strips_Y->cd(3); h_lead0_strips_Y3->Draw(); c0_strips_Y->cd(4); h_lead0_strips_Y4->Draw();

// Canvas για Leading Charge 1(Y)
TCanvas *c1_charge_Y = new TCanvas("c1_charge_Y", "1 Leading Charge Y Layers", 1200, 800);
c1_charge_Y->Divide(2,2);
c1_charge_Y->cd(1); h_lead1_charge_Y1->Draw(); c1_charge_Y->cd(2); h_lead1_charge_Y2->Draw();
c1_charge_Y->cd(3); h_lead1_charge_Y3->Draw(); c1_charge_Y->cd(4); h_lead1_charge_Y4->Draw();

// Canvas για Leading Strips 1(Y)
TCanvas *c1_strips_Y = new TCanvas("c1_strips_Y", "1 Leading Strips Y Layers", 1200, 800);
c1_strips_Y->Divide(2,2);
c1_strips_Y->cd(1); h_lead1_strips_Y1->Draw(); c1_strips_Y->cd(2); h_lead1_strips_Y2->Draw();
c1_strips_Y->cd(3); h_lead1_strips_Y3->Draw(); c1_strips_Y->cd(4); h_lead1_strips_Y4->Draw();

// Canvas για Leading Charge 2(Y)
TCanvas *c2_charge_Y = new TCanvas("c2_charge_Y", "2 Leading Charge Y Layers", 1200, 800);
c2_charge_Y->Divide(2,2);
c2_charge_Y->cd(1); h_lead2_charge_Y1->Draw(); c2_charge_Y->cd(2); h_lead2_charge_Y2->Draw();
c2_charge_Y->cd(3); h_lead2_charge_Y3->Draw(); c2_charge_Y->cd(4); h_lead2_charge_Y4->Draw();

// Canvas για Leading Strips 2(Y)
TCanvas *c2_strips_Y = new TCanvas("c2_strips_Y", "2 Leading Strips Y Layers", 1200, 800);
c2_strips_Y->Divide(2,2);
c2_strips_Y->cd(1); h_lead2_strips_Y1->Draw(); c2_strips_Y->cd(2); h_lead2_strips_Y2->Draw();
c2_strips_Y->cd(3); h_lead2_strips_Y3->Draw(); c2_strips_Y->cd(4); h_lead2_strips_Y4->Draw();
*/
// --- ΔΕΝ ΔΟΥΛΕΥΕΙ ΓΙΑ ΚΑΠΟΙΟ ΛΟΓΟ ΟΠΩΣ ΘΕΛΩ ΛΟΓΙΚΑ ΘΑ ΤΑ ΞΑΝΑ ΟΡΙΣΩ ΑΛΛΙΩΣ ---

/*
// --- Size vs Charge canvases (X1..X4 separated by angle region) ---
TCanvas *cSizeChargeX_60 = new TCanvas("cSizeChargeX_60", "Cluster Size vs Charge X (60-85)", 1200, 800);
cSizeChargeX_60->Divide(2,2);
cSizeChargeX_60->cd(1); h_cluSize_vs_cluCharge_60X1->Draw("colz");
cSizeChargeX_60->cd(2); h_cluSize_vs_cluCharge_60X2->Draw("colz");
cSizeChargeX_60->cd(3); h_cluSize_vs_cluCharge_60X3->Draw("colz");
cSizeChargeX_60->cd(4); h_cluSize_vs_cluCharge_60X4->Draw("colz");

TCanvas *cSizeChargeX_90 = new TCanvas("cSizeChargeX_90", "Cluster Size vs Charge X (85-95)", 1200, 800);
cSizeChargeX_90->Divide(2,2);
cSizeChargeX_90->cd(1); h_cluSize_vs_cluCharge_90X1->Draw("colz");
cSizeChargeX_90->cd(2); h_cluSize_vs_cluCharge_90X2->Draw("colz");
cSizeChargeX_90->cd(3); h_cluSize_vs_cluCharge_90X3->Draw("colz");
cSizeChargeX_90->cd(4); h_cluSize_vs_cluCharge_90X4->Draw("colz");

TCanvas *cSizeChargeX_120 = new TCanvas("cSizeChargeX_120", "Cluster Size vs Charge X (95-120)", 1200, 800);
cSizeChargeX_120->Divide(2,2);
cSizeChargeX_120->cd(1); h_cluSize_vs_cluCharge_120X1->Draw("colz");
cSizeChargeX_120->cd(2); h_cluSize_vs_cluCharge_120X2->Draw("colz");
cSizeChargeX_120->cd(3); h_cluSize_vs_cluCharge_120X3->Draw("colz");
cSizeChargeX_120->cd(4); h_cluSize_vs_cluCharge_120X4->Draw("colz");
*/
/*
// --- Size vs Charge canvases (Y1..Y4 separated by angle region) ---
TCanvas *cSizeChargeY_60 = new TCanvas("cSizeChargeY_60", "Cluster Size vs Charge Y (60-85)", 1200, 800);
cSizeChargeY_60->Divide(2,2);
cSizeChargeY_60->cd(1); h_cluSize_vs_cluCharge_60Y1->Draw("colz");
cSizeChargeY_60->cd(2); h_cluSize_vs_cluCharge_60Y2->Draw("colz");
cSizeChargeY_60->cd(3); h_cluSize_vs_cluCharge_60Y3->Draw("colz");
cSizeChargeY_60->cd(4); h_cluSize_vs_cluCharge_60Y4->Draw("colz");

TCanvas *cSizeChargeY_90 = new TCanvas("cSizeChargeY_90", "Cluster Size vs Charge Y (85-95)", 1200, 800);
cSizeChargeY_90->Divide(2,2);
cSizeChargeY_90->cd(1); h_cluSize_vs_cluCharge_90Y1->Draw("colz");
cSizeChargeY_90->cd(2); h_cluSize_vs_cluCharge_90Y2->Draw("colz");
cSizeChargeY_90->cd(3); h_cluSize_vs_cluCharge_90Y3->Draw("colz");
cSizeChargeY_90->cd(4); h_cluSize_vs_cluCharge_90Y4->Draw("colz");

TCanvas *cSizeChargeY_120 = new TCanvas("cSizeChargeY_120", "Cluster Size vs Charge Y (95-120)", 1200, 800);
cSizeChargeY_120->Divide(2,2);
cSizeChargeY_120->cd(1); h_cluSize_vs_cluCharge_120Y1->Draw("colz");
cSizeChargeY_120->cd(2); h_cluSize_vs_cluCharge_120Y2->Draw("colz");
cSizeChargeY_120->cd(3); h_cluSize_vs_cluCharge_120Y3->Draw("colz");
cSizeChargeY_120->cd(4); h_cluSize_vs_cluCharge_120Y4->Draw("colz");
*/

// --- Separate cluster size distributions per layer ---
TCanvas *cClusterSizeX = new TCanvas("cClusterSizeX", "Cluster Size X1-X4", 1200, 800);
cClusterSizeX->Divide(2,2);
cClusterSizeX->cd(1); cluster_length_Tmm1X->Draw();
cClusterSizeX->cd(2); cluster_length_Tmm2X->Draw();
cClusterSizeX->cd(3); cluster_length_Tmm3X->Draw();
cClusterSizeX->cd(4); cluster_length_Tmm4X->Draw();

TCanvas *cClusterSizeY = new TCanvas("cClusterSizeY", "Cluster Size Y1-Y4", 1200, 800);
cClusterSizeY->Divide(2,2);
cClusterSizeY->cd(1); cluster_length_Tmm1Y->Draw();
cClusterSizeY->cd(2); cluster_length_Tmm2Y->Draw();
cClusterSizeY->cd(3); cluster_length_Tmm3Y->Draw();
cClusterSizeY->cd(4); cluster_length_Tmm4Y->Draw();

// --- Separate cluster charge distributions per layer ---
TCanvas *cClusterChargeX = new TCanvas("cClusterChargeX", "Cluster Charge X1-X4", 1200, 800);
cClusterChargeX->Divide(2,2);
cClusterChargeX->cd(1); cl_Charge_Tmm1X->Draw();
cClusterChargeX->cd(2); cl_Charge_Tmm2X->Draw();
cClusterChargeX->cd(3); cl_Charge_Tmm3X->Draw();
cClusterChargeX->cd(4); cl_Charge_Tmm4X->Draw();

TCanvas *cClusterChargeY = new TCanvas("cClusterChargeY", "Cluster Charge Y1-Y4", 1200, 800);
cClusterChargeY->Divide(2,2);
cClusterChargeY->cd(1); cl_Charge_Tmm1Y->Draw();
cClusterChargeY->cd(2); cl_Charge_Tmm2Y->Draw();
cClusterChargeY->cd(3); cl_Charge_Tmm3Y->Draw();
cClusterChargeY->cd(4); cl_Charge_Tmm4Y->Draw();

/*
TCanvas *cX_charge_gt0 = new TCanvas("cX_charge_gt0", "Charge vs Angle X (nclus>0)", 1200, 800);
cX_charge_gt0->Divide(2,2);
cX_charge_gt0->cd(1); h_cluCharge_vs_AngleX1_gt0->Draw("colz");
cX_charge_gt0->cd(2); h_cluCharge_vs_AngleX2_gt0->Draw("colz");
cX_charge_gt0->cd(3); h_cluCharge_vs_AngleX3_gt0->Draw("colz");
cX_charge_gt0->cd(4); h_cluCharge_vs_AngleX4_gt0->Draw("colz");

TCanvas *cX_charge_gt1 = new TCanvas("cX_charge_gt1", "Charge vs Angle X (nclus>1)", 1200, 800);
cX_charge_gt1->Divide(2,2);
cX_charge_gt1->cd(1); h_cluCharge_vs_AngleX1_gt1->Draw("colz");
cX_charge_gt1->cd(2); h_cluCharge_vs_AngleX2_gt1->Draw("colz");
cX_charge_gt1->cd(3); h_cluCharge_vs_AngleX3_gt1->Draw("colz");
cX_charge_gt1->cd(4); h_cluCharge_vs_AngleX4_gt1->Draw("colz");

TCanvas *cX_charge_gt2 = new TCanvas("cX_charge_gt2", "Charge vs Angle X (nclus>2)", 1200, 800);
cX_charge_gt2->Divide(2,2);
cX_charge_gt2->cd(1); h_cluCharge_vs_AngleX1_gt2->Draw("colz");
cX_charge_gt2->cd(2); h_cluCharge_vs_AngleX2_gt2->Draw("colz");
cX_charge_gt2->cd(3); h_cluCharge_vs_AngleX3_gt2->Draw("colz");
cX_charge_gt2->cd(4); h_cluCharge_vs_AngleX4_gt2->Draw("colz");

TCanvas *cY_charge_gt0 = new TCanvas("cY_charge_gt0", "Charge vs Angle Y (nclus>0)", 1200, 800);
cY_charge_gt0->Divide(2,2);
cY_charge_gt0->cd(1); h_cluCharge_vs_AngleY1_gt0->Draw("colz");
cY_charge_gt0->cd(2); h_cluCharge_vs_AngleY2_gt0->Draw("colz");
cY_charge_gt0->cd(3); h_cluCharge_vs_AngleY3_gt0->Draw("colz");
cY_charge_gt0->cd(4); h_cluCharge_vs_AngleY4_gt0->Draw("colz");

TCanvas *cY_charge_gt1 = new TCanvas("cY_charge_gt1", "Charge vs Angle Y (nclus>1)", 1200, 800);
cY_charge_gt1->Divide(2,2);
cY_charge_gt1->cd(1); h_cluCharge_vs_AngleY1_gt1->Draw("colz");
cY_charge_gt1->cd(2); h_cluCharge_vs_AngleY2_gt1->Draw("colz");
cY_charge_gt1->cd(3); h_cluCharge_vs_AngleY3_gt1->Draw("colz");
cY_charge_gt1->cd(4); h_cluCharge_vs_AngleY4_gt1->Draw("colz");

TCanvas *cY_charge_gt2 = new TCanvas("cY_charge_gt2", "Charge vs Angle Y (nclus>2)", 1200, 800);
cY_charge_gt2->Divide(2,2);
cY_charge_gt2->cd(1); h_cluCharge_vs_AngleY1_gt2->Draw("colz");
cY_charge_gt2->cd(2); h_cluCharge_vs_AngleY2_gt2->Draw("colz");
cY_charge_gt2->cd(3); h_cluCharge_vs_AngleY3_gt2->Draw("colz");
cY_charge_gt2->cd(4); h_cluCharge_vs_AngleY4_gt2->Draw("colz");
*/
/*
// --- Καμβάδες για τον άξονα Χ ---
TCanvas *clszanglx0 = new TCanvas("clszanglx", "Cluster Size 0 vs Angle X", 800, 600);
clszanglx0->Divide(2, 2);
clszanglx0->cd(1); h_cluSize_vs_AngleX1_gt0->Draw("hist");
clszanglx0->cd(2); h_cluSize_vs_AngleX2_gt0->Draw("hist");
clszanglx0->cd(3); h_cluSize_vs_AngleX3_gt0->Draw("hist");
clszanglx0->cd(4); h_cluSize_vs_AngleX4_gt0->Draw("hist");

TCanvas *clszanglx1 = new TCanvas("clszanglx1", "Cluster Size 1 vs Angle X", 800, 600);
clszanglx1->Divide(2, 2);
clszanglx1->cd(1); h_cluSize_vs_AngleX1_gt1->Draw("hist");
clszanglx1->cd(2); h_cluSize_vs_AngleX2_gt1->Draw("hist");
clszanglx1->cd(3); h_cluSize_vs_AngleX3_gt1->Draw("hist");
clszanglx1->cd(4); h_cluSize_vs_AngleX4_gt1->Draw("hist");

TCanvas *clszanglx2 = new TCanvas("clszanglx2", "Cluster Size 2 vs Angle X", 800, 600);
clszanglx2->Divide(2, 2);
clszanglx2->cd(1); h_cluSize_vs_AngleX1_gt2->Draw("hist");
clszanglx2->cd(2); h_cluSize_vs_AngleX2_gt2->Draw("hist");
clszanglx2->cd(3); h_cluSize_vs_AngleX3_gt2->Draw("hist");
clszanglx2->cd(4); h_cluSize_vs_AngleX4_gt2->Draw("hist");

TCanvas *clszangly0 = new TCanvas("clszangly0", "Cluster Size 0 vs Angle Y", 800, 600);
clszangly0->Divide(2, 2);
clszangly0->cd(1); h_cluSize_vs_AngleY1_gt0->Draw("hist");
clszangly0->cd(2); h_cluSize_vs_AngleY2_gt0->Draw("hist");
clszangly0->cd(3); h_cluSize_vs_AngleY3_gt0->Draw("hist");
clszangly0->cd(4); h_cluSize_vs_AngleY4_gt0->Draw("hist");

TCanvas *clszangly1 = new TCanvas("clszangly1", "Cluster Size 1 vs Angle Y", 800, 600);
clszangly1->Divide(2, 2);
clszangly1->cd(1); h_cluSize_vs_AngleY1_gt1->Draw("hist");
clszangly1->cd(2); h_cluSize_vs_AngleY2_gt1->Draw("hist");
clszangly1->cd(3); h_cluSize_vs_AngleY3_gt1->Draw("hist");
clszangly1->cd(4); h_cluSize_vs_AngleY4_gt1->Draw("hist");

TCanvas *clszangly2 = new TCanvas("clszangly2", "Cluster Size 2 vs Angle Y", 800, 600);
clszangly2->Divide(2, 2);
clszangly2->cd(1); h_cluSize_vs_AngleY1_gt2->Draw("hist");
clszangly2->cd(2); h_cluSize_vs_AngleY2_gt2->Draw("hist");
clszangly2->cd(3); h_cluSize_vs_AngleY3_gt2->Draw("hist");
clszangly2->cd(4); h_cluSize_vs_AngleY4_gt2->Draw("hist");
*/

TCanvas *cX_ang = new TCanvas("cX_ang", "X-Axis: CLSZX Angle Analysis", 1000, 400);
cX_ang->Divide(3, 1);        // three histos now
cX_ang->cd(1); h_cluSize_vs_AngleX_gt0->Draw("hist");
cX_ang->cd(2); h_cluSize_vs_AngleX_gt1->Draw("hist");
cX_ang->cd(3); h_cluSize_vs_AngleX_gt2->Draw("hist");

TCanvas *cY_ang = new TCanvas("cY_ang", "Y-Axis: CLSZY Angle Analysis", 1000, 400);
cY_ang->Divide(3, 1);
cY_ang->cd(1); h_cluSize_vs_AngleY_gt0->Draw("hist");
cY_ang->cd(2); h_cluSize_vs_AngleY_gt1->Draw("hist");
cY_ang->cd(3); h_cluSize_vs_AngleY_gt2->Draw("hist");

TCanvas *cX_angch = new TCanvas("cX_angch", "X-Axis: Cluster Angle Analysis", 1000, 400);
cX_angch->Divide(3, 1);
cX_angch->cd(1); h_cluCharge_vs_AngleX_gt0->Draw("hist");
cX_angch->cd(2); h_cluCharge_vs_AngleX_gt1->Draw("hist");
cX_angch->cd(3); h_cluCharge_vs_AngleX_gt2->Draw("hist");

TCanvas *cY_angch = new TCanvas("cY_angch", "Y-Axis: Cluster Angle Analysis", 1000, 400);
cY_angch->Divide(3, 1);
cY_angch->cd(1); h_cluCharge_vs_AngleY_gt0->Draw("hist");
cY_angch->cd(2); h_cluCharge_vs_AngleY_gt1->Draw("hist");
cY_angch->cd(3); h_cluCharge_vs_AngleY_gt2->Draw("hist");


// --- αυτα τα εχω κανει. πρεπει να ειναι σωστα ---

TCanvas *cX_led = new TCanvas("cX_led", "X-Axis: Leading Charge per Strip", 800, 800);
cX_led->Divide(2, 2);
cX_led->cd(1); h_leadClu_Charge_vs_Strips1X->Draw("hist");
cX_led->cd(2); h_leadClu_Charge_vs_Strips2X->Draw("hist");
cX_led->cd(3); h_leadClu_Charge_vs_Strips3X->Draw("hist");
cX_led->cd(4); h_leadClu_Charge_vs_Strips4X->Draw("hist");

TCanvas *cY_led = new TCanvas("cY_led", "Y-Axis: Leading Charge per Strip", 800, 800);
cY_led->Divide(2, 2);
cY_led->cd(1); h_leadClu_Charge_vs_Strips1Y->Draw("hist");
cY_led->cd(2); h_leadClu_Charge_vs_Strips2Y->Draw("hist");
cY_led->cd(3); h_leadClu_Charge_vs_Strips3Y->Draw("hist");
cY_led->cd(4); h_leadClu_Charge_vs_Strips4Y->Draw("hist");
/*
TCanvas *cStripXY_60_80 = new TCanvas("cStripXY_60_80", "Strip Correlation Xn vs Yn (60-80 deg, nclus>0)", 1200, 900);
cStripXY_60_80->Divide(2, 2);
cStripXY_60_80->cd(1); h_stripXY_60_80_X1Y1->Draw("colz");
cStripXY_60_80->cd(2); h_stripXY_60_80_X2Y2->Draw("colz");
cStripXY_60_80->cd(3); h_stripXY_60_80_X3Y3->Draw("colz");
cStripXY_60_80->cd(4); h_stripXY_60_80_X4Y4->Draw("colz");

TCanvas *cStripXY_60_80_cl1 = new TCanvas("cStripXY_60_80_cl1", "Strip Correlation Xn vs Yn (60-80 deg, nclus>1)", 1200, 900);
cStripXY_60_80_cl1->Divide(2, 2);
cStripXY_60_80_cl1->cd(1); h_stripXY_60_80_X1Y1_cl1->Draw("colz");
cStripXY_60_80_cl1->cd(2); h_stripXY_60_80_X2Y2_cl1->Draw("colz");
cStripXY_60_80_cl1->cd(3); h_stripXY_60_80_X3Y3_cl1->Draw("colz");
cStripXY_60_80_cl1->cd(4); h_stripXY_60_80_X4Y4_cl1->Draw("colz");

TCanvas *cStripXY_60_80_cl2 = new TCanvas("cStripXY_60_80_cl2", "Strip Correlation Xn vs Yn (60-80 deg, nclus>2)", 1200, 900);
cStripXY_60_80_cl2->Divide(2, 2);
cStripXY_60_80_cl2->cd(1); h_stripXY_60_80_X1Y1_cl2->Draw("colz");
cStripXY_60_80_cl2->cd(2); h_stripXY_60_80_X2Y2_cl2->Draw("colz");
cStripXY_60_80_cl2->cd(3); h_stripXY_60_80_X3Y3_cl2->Draw("colz");
cStripXY_60_80_cl2->cd(4); h_stripXY_60_80_X4Y4_cl2->Draw("colz");

TCanvas *cStripXY_85_95 = new TCanvas("cStripXY_85_95", "Strip Correlation Xn vs Yn (85-95 deg, nclus>0)", 1200, 900);
cStripXY_85_95->Divide(2, 2);
cStripXY_85_95->cd(1); h_stripXY_85_95_X1Y1->Draw("colz");
cStripXY_85_95->cd(2); h_stripXY_85_95_X2Y2->Draw("colz");
cStripXY_85_95->cd(3); h_stripXY_85_95_X3Y3->Draw("colz");
cStripXY_85_95->cd(4); h_stripXY_85_95_X4Y4->Draw("colz");

TCanvas *cStripXY_85_95_cl1 = new TCanvas("cStripXY_85_95_cl1", "Strip Correlation Xn vs Yn (85-95 deg, nclus>1)", 1200, 900);
cStripXY_85_95_cl1->Divide(2, 2);
cStripXY_85_95_cl1->cd(1); h_stripXY_85_95_X1Y1_cl1->Draw("colz");
cStripXY_85_95_cl1->cd(2); h_stripXY_85_95_X2Y2_cl1->Draw("colz");
cStripXY_85_95_cl1->cd(3); h_stripXY_85_95_X3Y3_cl1->Draw("colz");
cStripXY_85_95_cl1->cd(4); h_stripXY_85_95_X4Y4_cl1->Draw("colz");

TCanvas *cStripXY_85_95_cl2 = new TCanvas("cStripXY_85_95_cl2", "Strip Correlation Xn vs Yn (85-95 deg, nclus>2)", 1200, 900);
cStripXY_85_95_cl2->Divide(2, 2);
cStripXY_85_95_cl2->cd(1); h_stripXY_85_95_X1Y1_cl2->Draw("colz");
cStripXY_85_95_cl2->cd(2); h_stripXY_85_95_X2Y2_cl2->Draw("colz");
cStripXY_85_95_cl2->cd(3); h_stripXY_85_95_X3Y3_cl2->Draw("colz");
cStripXY_85_95_cl2->cd(4); h_stripXY_85_95_X4Y4_cl2->Draw("colz");

TCanvas *cStripXY_100_120 = new TCanvas("cStripXY_100_120", "Strip Correlation Xn vs Yn (100-120 deg, nclus>0)", 1200, 900);
cStripXY_100_120->Divide(2, 2);
cStripXY_100_120->cd(1); h_stripXY_100_120_X1Y1->Draw("colz");
cStripXY_100_120->cd(2); h_stripXY_100_120_X2Y2->Draw("colz");
cStripXY_100_120->cd(3); h_stripXY_100_120_X3Y3->Draw("colz");
cStripXY_100_120->cd(4); h_stripXY_100_120_X4Y4->Draw("colz");

TCanvas *cStripXY_100_120_cl1 = new TCanvas("cStripXY_100_120_cl1", "Strip Correlation Xn vs Yn (100-120 deg, nclus>1)", 1200, 900);
cStripXY_100_120_cl1->Divide(2, 2);
cStripXY_100_120_cl1->cd(1); h_stripXY_100_120_X1Y1_cl1->Draw("colz");
cStripXY_100_120_cl1->cd(2); h_stripXY_100_120_X2Y2_cl1->Draw("colz");
cStripXY_100_120_cl1->cd(3); h_stripXY_100_120_X3Y3_cl1->Draw("colz");
cStripXY_100_120_cl1->cd(4); h_stripXY_100_120_X4Y4_cl1->Draw("colz");

TCanvas *cStripXY_100_120_cl2 = new TCanvas("cStripXY_100_120_cl2", "Strip Correlation Xn vs Yn (100-120 deg, nclus>2)", 1200, 900);
cStripXY_100_120_cl2->Divide(2, 2);
cStripXY_100_120_cl2->cd(1); h_stripXY_100_120_X1Y1_cl2->Draw("colz");
cStripXY_100_120_cl2->cd(2); h_stripXY_100_120_X2Y2_cl2->Draw("colz");
cStripXY_100_120_cl2->cd(3); h_stripXY_100_120_X3Y3_cl2->Draw("colz");
cStripXY_100_120_cl2->cd(4); h_stripXY_100_120_X4Y4_cl2->Draw("colz");
*/


TCanvas *c0_charge_2D = new TCanvas("c0_charge_2D", "0 Cluster Charge Distribution X", 1600, 800);
c0_charge_2D->Divide(4,2);
c0_charge_2D->cd(1); h_charge_strip_1X1->Draw("hist");
c0_charge_2D->cd(2); h_charge_strip_1X2->Draw("hist");
c0_charge_2D->cd(3); h_charge_strip_1X3->Draw("hist");
c0_charge_2D->cd(4); h_charge_strip_1X4->Draw("hist");
c0_charge_2D->cd(5); h_charge_strip_1Y1->Draw("hist");
c0_charge_2D->cd(6); h_charge_strip_1Y2->Draw("hist");
c0_charge_2D->cd(7); h_charge_strip_1Y3->Draw("hist");
c0_charge_2D->cd(8); h_charge_strip_1Y4->Draw("hist");
/*
TCanvas *c1_charge_2D = new TCanvas("c1_charge_2D", "1 Cluster Charge Distribution X", 1600, 800);
c1_charge_2D->Divide(4,2);
c1_charge_2D->cd(1); h_charge_strip_2X1->Draw("hist");
c1_charge_2D->cd(2); h_charge_strip_2X2->Draw("hist");
c1_charge_2D->cd(3); h_charge_strip_2X3->Draw("hist");
c1_charge_2D->cd(4); h_charge_strip_2X4->Draw("hist");
c1_charge_2D->cd(5); h_charge_strip_2Y1->Draw("hist");
c1_charge_2D->cd(6); h_charge_strip_2Y2->Draw("hist");
c1_charge_2D->cd(7); h_charge_strip_2Y3->Draw("hist");
c1_charge_2D->cd(8); h_charge_strip_2Y4->Draw("hist");

TCanvas *c3_charge_2D = new TCanvas("c3_charge_2D", "3 Cluster Charge Distribution X", 1600, 800);
c3_charge_2D->Divide(4,2);
c3_charge_2D->cd(1); h_charge_strip_3X1->Draw("hist");
c3_charge_2D->cd(2); h_charge_strip_3X2->Draw("hist");
c3_charge_2D->cd(3); h_charge_strip_3X3->Draw("hist");
c3_charge_2D->cd(4); h_charge_strip_3X4->Draw("hist");
c3_charge_2D->cd(5); h_charge_strip_3Y1->Draw("hist");
c3_charge_2D->cd(6); h_charge_strip_3Y2->Draw("hist");
c3_charge_2D->cd(7); h_charge_strip_3Y3->Draw("hist");
c3_charge_2D->cd(8); h_charge_strip_3Y4->Draw("hist");
*/

TCanvas *c_slopes = new TCanvas("c_slopes", "Slope Analysis", 1200, 800);
c_slopes->Divide(3, 2);
c_slopes->cd(1); h_slope_X->Draw();
c_slopes->cd(2); h_intercept_X->Draw();
c_slopes->cd(3); h_slope_vs_intercept_X->Draw("colz");
c_slopes->cd(4); h_slope_Y->Draw();
c_slopes->cd(5); h_intercept_Y->Draw();
c_slopes->cd(6); h_slope_vs_intercept_Y->Draw("colz");


// ==========================================
// --- 3. DRAWING CANVASES ---
// ==========================================

// === 6 2D GRAPHS: Angle vs nclus ===
TCanvas *cX_nclus_cl0 = new TCanvas("cX_nclus_cl0", "nclusX>0 (all layers > 0)", 900, 400);
cX_nclus_cl0->Divide(2, 1);
cX_nclus_cl0->cd(1); h_clus_angle_X_cl0->Draw("colz");
cX_nclus_cl0->cd(2); h_charge_angle_X_cl0->Draw("colz");

TCanvas *cX_nclus_cl1 = new TCanvas("cX_nclus_cl1", "nclusX>1 (all layers > 1)", 900, 400);
cX_nclus_cl1->Divide(2, 1);
cX_nclus_cl1->cd(1); h_clus_angle_X_cl1->Draw("colz");
cX_nclus_cl1->cd(2); h_charge_angle_X_cl1->Draw("colz");

TCanvas *cX_nclus_cl2 = new TCanvas("cX_nclus_cl2", "nclusX>2 (all layers > 2)", 900, 400);
cX_nclus_cl2->Divide(2, 1);
cX_nclus_cl2->cd(1); h_clus_angle_X_cl2->Draw("colz");
cX_nclus_cl2->cd(2); h_charge_angle_X_cl2->Draw("colz");

TCanvas *cY_nclus_cl0 = new TCanvas("cY_nclus_cl0", "nclusY>0 (all layers > 0)", 900, 400);
cY_nclus_cl0->Divide(2, 1);
cY_nclus_cl0->cd(1); h_clus_angle_Y_cl0->Draw("colz");
cY_nclus_cl0->cd(2); h_charge_angle_Y_cl0->Draw("colz");

TCanvas *cY_nclus_cl1 = new TCanvas("cY_nclus_cl1", "nclusY>1 (all layers > 1)", 900, 400);
cY_nclus_cl1->Divide(2, 1);
cY_nclus_cl1->cd(1); h_clus_angle_Y_cl1->Draw("colz");
cY_nclus_cl1->cd(2); h_charge_angle_Y_cl1->Draw("colz");

TCanvas *cY_nclus_cl2 = new TCanvas("cY_nclus_cl2", "nclusY>2 (all layers > 2)", 900, 400);
cY_nclus_cl2->Divide(2, 1);
cY_nclus_cl2->cd(1); h_clus_angle_Y_cl2->Draw("colz");
cY_nclus_cl2->cd(2); h_charge_angle_Y_cl2->Draw("colz");

TCanvas *cAngles = new TCanvas("cAngles", "Leading-Strip Angle Diagrams", 150, 0, 1200, 800);
cAngles->Divide(2,2);
cAngles->cd(1); slopeX->Draw();
cAngles->cd(2); slopeY->Draw();
cAngles->cd(3); slopeXY->Draw("colz");
cAngles->cd(4); h_map_XY->Draw("colz");

TCanvas *cAnglesChi2Cut = new TCanvas("cAnglesChi2Cut", "Angle Diagrams for #chi^{2}", 150, 0, 1200, 800);
cAnglesChi2Cut->Divide(2,2);
cAnglesChi2Cut->cd(1); h_slopeX_chi2lt->Draw("hist");
cAnglesChi2Cut->cd(2); h_slopeY_chi2lt->Draw("hist");
cAnglesChi2Cut->cd(3); h_slopeXY_chi2lt->Draw("colz");
cAnglesChi2Cut->cd(4); h_track_chi2->Draw("hist");

TCanvas *cTrackComboAngles = new TCanvas("cTrackComboAngles", "Track Angle Combinations", 150, 0, 1200, 800);
cTrackComboAngles->Divide(3,2);
cTrackComboAngles->cd(1); hslope4321->Draw("hist");
cTrackComboAngles->cd(2); hslope432->Draw("hist");
cTrackComboAngles->cd(3); hslope43->Draw("hist");
cTrackComboAngles->cd(4); hslope321->Draw("hist");
cTrackComboAngles->cd(5); hslope32->Draw("hist");
cTrackComboAngles->cd(6); hslope21->Draw("hist");

TCanvas *cTrackProjections = new TCanvas("cTrackProjections", "Track Projections", 150, 0, 1200, 800);
cTrackProjections->Divide(2,2);
cTrackProjections->cd(1); projection1->Draw("hist");
cTrackProjections->cd(2); projection2->Draw("hist");
cTrackProjections->cd(3); projection3->Draw("hist");
cTrackProjections->cd(4); projection4->Draw("hist");

TCanvas *cDeltaPos = new TCanvas("cDeltaPos", "Delta Position", 150, 0, 1200, 800);
cDeltaPos->Divide(2,2);
cDeltaPos->cd(1); Delta_pos1->Draw("hist");
cDeltaPos->cd(2); Delta_pos2->Draw("hist");
cDeltaPos->cd(3); Delta_pos3->Draw("hist");
cDeltaPos->cd(4); Delta_pos4->Draw("hist");

TCanvas *cDeltaAngle = new TCanvas("cDeltaAngle", "Delta Angle", 150, 0, 1200, 800);
cDeltaAngle->Divide(2,2);
cDeltaAngle->cd(1); delta_angle1->Draw("hist");
cDeltaAngle->cd(2); delta_angle2->Draw("hist");
cDeltaAngle->cd(3); delta_angle3->Draw("hist");
cDeltaAngle->cd(4); delta_angle4->Draw("hist");

TCanvas *cDeltaPosAngle = new TCanvas("cDeltaPosAngle", "Delta Position vs Angle", 150, 0, 1200, 800);
cDeltaPosAngle->Divide(2,2);
cDeltaPosAngle->cd(1); Delta_pos_angle1->Draw("colz");
cDeltaPosAngle->cd(2); Delta_pos_angle2->Draw("colz");
cDeltaPosAngle->cd(3); Delta_pos_angle3->Draw("colz");
cDeltaPosAngle->cd(4); Delta_pos_angle4->Draw("colz");

TCanvas *cSlopeTpcAngle = new TCanvas("cSlopeTpcAngle", "Track Slope vs TPC Angle", 150, 0, 1200, 800);
cSlopeTpcAngle->Divide(2,2);
cSlopeTpcAngle->cd(1); slope_tpc_ang1->Draw("colz");
cSlopeTpcAngle->cd(2); slope_tpc_ang2->Draw("colz");
cSlopeTpcAngle->cd(3); slope_tpc_ang3->Draw("colz");
cSlopeTpcAngle->cd(4); slope_tpc_ang4->Draw("colz");

TCanvas *cCentroidUtpc = new TCanvas("cCentroidUtpc", "Centroid vs uTPC Position", 150, 0, 1200, 800);
cCentroidUtpc->Divide(2,2);
cCentroidUtpc->cd(1); centroid_utpc1->Draw("colz");
cCentroidUtpc->cd(2); centroid_utpc2->Draw("colz");
cCentroidUtpc->cd(3); centroid_utpc3->Draw("colz");
cCentroidUtpc->cd(4); centroid_utpc4->Draw("colz");
/*
TCanvas *cChargePerStrip = new TCanvas("cChargePerStrip", "Charge per Strip", 150, 0, 1200, 800);
cChargePerStrip->Divide(2,2);
cChargePerStrip->cd(1); charge_per_strip_tmm1->Draw("colz");
cChargePerStrip->cd(2); charge_per_strip_tmm2->Draw("colz");
cChargePerStrip->cd(3); charge_per_strip_tmm3->Draw("colz");
cChargePerStrip->cd(4); charge_per_strip_tmm4->Draw("colz");
*/

TCanvas *cTimeTPC = new TCanvas("cTimeTPC", "TPC Time", 150, 0, 1200, 800);
cTimeTPC->Divide(2,2);
cTimeTPC->cd(1); timeT1->Draw("hist");
cTimeTPC->cd(2); timeT2->Draw("hist");
cTimeTPC->cd(3); timeT3->Draw("hist");
cTimeTPC->cd(4); timeT4->Draw("hist");

TCanvas *cUtpcSlope = new TCanvas("cUtpcSlope", "uTPC Slope", 150, 0, 1200, 800);
cUtpcSlope->Divide(2,2);
cUtpcSlope->cd(1); utpc_slop1->Draw("hist");
cUtpcSlope->cd(2); utpc_slop2->Draw("hist");
cUtpcSlope->cd(3); utpc_slop3->Draw("hist");
cUtpcSlope->cd(4); utpc_slop4->Draw("hist");

    // --- Canvas για την Κατανομή Landau ---
    TCanvas *cLandauX = new TCanvas("cLandauX", "Landau Distribution - X Layers", 1200, 800);
    cLandauX->Divide(2, 2);
    cLandauX->cd(1); h_landau_charge_X1->Draw("hist");
    cLandauX->cd(2); h_landau_charge_X2->Draw("hist");
    cLandauX->cd(3); h_landau_charge_X3->Draw("hist");
    cLandauX->cd(4); h_landau_charge_X4->Draw("hist");

    TCanvas *cLandauY = new TCanvas("cLandauY", "Landau Distribution - Y Layers", 1200, 800);
    cLandauY->Divide(2, 2);
    cLandauY->cd(1); h_landau_charge_Y1->Draw("hist");
    cLandauY->cd(2); h_landau_charge_Y2->Draw("hist");
    cLandauY->cd(3); h_landau_charge_Y3->Draw("hist");
    cLandauY->cd(4); h_landau_charge_Y4->Draw("hist");
    
TCanvas *c_delta_angle = new TCanvas("c_delta_angle", "Chi2 vs Deviation from 90", 600, 600);
c_delta_angle->Divide(1, 1);

TCanvas *c1_delta_angle = new TCanvas("c1_delta_angle", "Chi2 vs Deviation from 90", 600, 600);
c1_delta_angle->Divide(1, 1);
// Left side: The 3D Heatmap (Residuals as colors)
c1_delta_angle->cd(1);


// 1. Declare your 2D matrix ONCE
    double strip_noise_floor[8][640]; 

    // 2. Put your 8 TH2F histogram pointers into an array
    TH2F* chargeHists[8] = {
        h_charge_strip_1X1, h_charge_strip_1X2, h_charge_strip_1X3, h_charge_strip_1X4,
        h_charge_strip_1Y1, h_charge_strip_1Y2, h_charge_strip_1Y3, h_charge_strip_1Y4
    };

    // 3. Populate the matrix
    extractLowestChargeMatrix(chargeHists, strip_noise_floor);

    // 4. Print out or export the matrix to a text file
    std::ofstream outFile("strip_noise_thresholds.txt"); 
    outFile << "# Layer Strip LowestCharge_ADC" << std::endl;
    
    for (int lay = 0; lay < 8; lay++) {
        for (int strp = 0; strp < 640; strp++) {
            outFile << lay << " " << strp << " " << strip_noise_floor[lay][strp] << "\n";
        }
    }
    outFile.close(); 
    // =================================================================
/*
TCanvas *cAngleVsCharge234 = new TCanvas("cAngleVsCharge234", "Angle vs Cluster Charge 234", 150, 0, 1200, 800);
cAngleVsCharge234->Divide(2,2);
cAngleVsCharge234->cd(1); h1_angle_vs_q_234->Draw("colz");
cAngleVsCharge234->cd(2); h2_angle_vs_q_234->Draw("colz");
cAngleVsCharge234->cd(3); h3_angle_vs_q_234->Draw("colz");
cAngleVsCharge234->cd(4); h4_angle_vs_q_234->Draw("colz");

    TCanvas *c1 = new TCanvas("c1","strips-X", 50,0,1200,800);
        c1->Divide(2,2);
        c1->cd(1);h_strps_Tmm1X->Draw();
        c1->cd(2);h_strps_Tmm2X->Draw();
        c1->cd(3);h_strps_Tmm3X->Draw();
        c1->cd(4);h_strps_Tmm4X->Draw();

    gStyle->SetOptStat(1111);

    TCanvas *c2 = new TCanvas("c2","strips-Y", 150,0,1200,800);
        c2->Divide(2,2);
        c2->cd(1);h_strps_Tmm1Y->Draw();
        c2->cd(2);h_strps_Tmm2Y->Draw();
        c2->cd(3);h_strps_Tmm3Y->Draw();
        c2->cd(4);h_strps_Tmm4Y->Draw();

    TCanvas *c3 = new TCanvas("c3","EVENT TYPE ", 150,0,1600,800);
    c3->Divide(2,1);
    c3->cd(1); h_Event_typeX->Draw();
    c3->cd(2); h_Event_typeY->Draw();

    TCanvas *c100 = new TCanvas("c100","EVENT TYPE ", 150,0,1600,800);
   // c100->Divide(2,1);
    //c100->cd(1);
    h_Event_typeXY->Draw();
    //c100->cd(2); h_Event_typeY->Draw();

    
 //   gStyle->SetPadTickX(1);
 //   gStyle->SetPadTickY(1);
    gROOT->ForceStyle();
    
}




    TCanvas *c12 = new TCanvas("c12","Cluster/evnt X", 50,0,800,800);
        c12->Divide(2,2);
        c12->cd(1);clu_per_eventTmm1X->Draw();
        c12->cd(2);clu_per_eventTmm2X->Draw();
        c12->cd(3);clu_per_eventTmm3X->Draw();
        c12->cd(4);clu_per_eventTmm4X->Draw();
    TCanvas *c13 = new TCanvas("c13","Cluster/evnt Y", 900,0,800,800);
        c13->Divide(2,2);
        c13->cd(1);clu_per_eventTmm1Y->Draw();
        c13->cd(2);clu_per_eventTmm2Y->Draw();
        c13->cd(3);clu_per_eventTmm3Y->Draw();
        c13->cd(4);clu_per_eventTmm4Y->Draw();



    
    TCanvas *c91 = new TCanvas("c91","SLOPEs ", 100,0,1200, 800);
    c91->Divide(2,1);
    c91->cd(1);
    //hslope4321->SetLineWidth(3);
    slopeX->SetLineColor(2);         slopeX->Draw();
    c91->cd(2);
    slopeY->SetLineColor(4);         slopeY->Draw();

    TCanvas *c92 = new TCanvas("c92","SLOPEs ", 100,0, 800, 800);
    slopeXY->Draw();

    TCanvas *c200 = new TCanvas("c200"," MAP ", 150,0,800,800);
    h_map_XY->Draw();

  
    c10->Print("plots.pdf(","Title:page1");
    c11->Print("plots.pdf","Title:page2");
    c12->Print("plots.pdf","Title:page3");
    c13->Print("plots.pdf","Title:page4");
    c14->Print("plots.pdf","Title:page5");
    c15->Print("plots.pdf","Title:page6");
//    c100->Print("plots.pdf","Title:page7");
    c91->Print("plots.pdf)","Title:page8");
 //   c200->Print("plots.pdf)","Title:page9");

    TFile f1("outplots.root", "RECREATE"); //Open file to store histos ;
    h_map_XY->Write();                      // Write histo to file ;
    slopeX->Write();
    slopeY->Write();
    f1.Close();
    
       gROOT->ForceStyle();

    
}
void mmanalys::plot_Residuals(){

    gROOT->ForceStyle(kTRUE);
    gStyle->SetOptStat(111111);
    gStyle->SetHistLineWidth(5);

  
 
  /*
    TCanvas *c91 = new TCanvas("c91","SLOPEs ", 100,0,1200, 800);
 //   c91->Divide(3,2);
    c91->cd(1);
    //hslope4321->SetLineWidth(3);
    hslope4321->SetLineColor(2);    hslope4321->Draw();
*/
    // 1. Ορίζουμε τη συνάρτηση της Διπλής Gaussian μία φορά στην αρχή
    TF1 *dGaus = new TF1("dGaus", "gaus(0) + gaus(3)", -0.8, 0.8);
    dGaus->SetLineColor(2); // Κόκκινο χρώμα για το συνολικό fit
    dGaus->SetParNames("Const1", "Mean1", "Sigma1", "Const2", "Mean2", "Sigma2");

    // ---------------- RESIDUALS X ----------------
    TCanvas *c92 = new TCanvas("c92","Residuals ", 100,0,1200, 1200);
    c92->Divide(2,2);
    
    TH1* histX[4] = {Residual1, Residual2, Residual3, Residual4};
    
    for(int i = 0; i < 4; i++) {
        c92->cd(i+1);
        histX[i]->Draw();
        
        // Δυναμικά όρια: παίρνουμε τα όρια Χ από το ίδιο το ιστόγραμμα
        double minX = histX[i]->GetXaxis()->GetXmin();
        double maxX = histX[i]->GetXaxis()->GetXmax();
        
        TString fNameX = Form("fitX_%d", i);
        TF1 *dGausX = new TF1(fNameX, "gaus(0) + gaus(3)", minX, maxX);
        
        dGausX->SetLineColor(2); // Κόκκινο χρώμα
        dGausX->SetLineWidth(3); // Πιο χοντρή γραμμή για να φαίνεται ξεκάθαρα
        
        double maxVal  = histX[i]->GetMaximum();
        double meanVal = histX[i]->GetMean();
        double rmsVal  = histX[i]->GetRMS();
        
        // Προστασία: Αν το ιστόγραμμα είναι άδειο, να μην κρασάρει η διαίρεση
        if (rmsVal == 0) rmsVal = 0.1; 
        
        
        dGausX->SetParameters(maxVal, 0.0, 0.1, maxVal/5.0, 0.0, 0.5);
        
        // ΑΛΛΑΓΗ: Βγάλαμε το "Q" και βάλαμε το "+"
        histX[i]->Fit(dGausX, "+"); 
        
        // Αναγκαστική σχεδίαση της γραμμής πάνω στο ιστόγραμμα
        TF1 *g1 = new TF1(Form("g1X_%d", i), "gaus", minX, maxX);
        TF1 *g2 = new TF1(Form("g2X_%d", i), "gaus", minX, maxX);
        
        // Αντιγράφουμε τις παραμέτρους από το επιτυχημένο συνολικό Fit
        g1->SetParameters(dGausX->GetParameter(0), dGausX->GetParameter(1), dGausX->GetParameter(2));
        g2->SetParameters(dGausX->GetParameter(3), dGausX->GetParameter(4), dGausX->GetParameter(5));
        
        // Ζωγραφίζουμε τη Στενή (πυρήνας) - π.χ. Ματζέντα χρώμα με διακεκομμένη γραμμή
        g1->SetLineColor(6); 
        g1->SetLineStyle(2); 
        g1->Draw("SAME");
        
        // Ζωγραφίζουμε τη Φαρδιά (ουρές) - π.χ. Πράσινο χρώμα με διακεκομμένη γραμμή
        g2->SetLineColor(3); 
        g2->SetLineStyle(2); 
        g2->Draw("SAME");
    }

    // ---------------- RESIDUALS Y ----------------
    TCanvas *c94 = new TCanvas("c94","ResidualsY ", 100,0,1200, 1200);
    c94->Divide(2,2);
    
    TH1* histY[4] = {Residual1Y, Residual2Y, Residual3Y, Residual4Y};
    
    for(int i = 0; i < 4; i++) {
        c94->cd(i+1);
        histY[i]->Draw();
        
        double minX = histY[i]->GetXaxis()->GetXmin();
        double maxX = histY[i]->GetXaxis()->GetXmax();
        
        TString fNameY = Form("fitY_%d", i);
        TF1 *dGausY = new TF1(fNameY, "gaus(0) + gaus(3)", minX, maxX);
        
        dGausY->SetLineColor(4); // Μπλε χρώμα για το Y
        dGausY->SetLineWidth(3);
        
        double maxVal  = histY[i]->GetMaximum();
        double meanVal = histY[i]->GetMean();
        double rmsVal  = histY[i]->GetRMS();
        
        if (rmsVal == 0) rmsVal = 0.1;
        

        dGausY->SetParameters(maxVal, 0.0, 0.1, maxVal/5.0, 0.0, 0.5);
        
        histY[i]->Fit(dGausY, "+");
        
        TF1 *g1 = new TF1(Form("g1X_%d", i), "gaus", minX, maxX);
        TF1 *g2 = new TF1(Form("g2X_%d", i), "gaus", minX, maxX);
        
        // Αντιγράφουμε τις παραμέτρους από το επιτυχημένο συνολικό Fit
        g1->SetParameters(dGausY->GetParameter(0), dGausY->GetParameter(1), dGausY->GetParameter(2));
        g2->SetParameters(dGausY->GetParameter(3), dGausY->GetParameter(4), dGausY->GetParameter(5));
        
        // Ζωγραφίζουμε τη Στενή (πυρήνας) - π.χ. Ματζέντα χρώμα με διακεκομμένη γραμμή
        g1->SetLineColor(6); 
        g1->SetLineStyle(2); 
        g1->Draw("SAME");
        
        // Ζωγραφίζουμε τη Φαρδιά (ουρές) - π.χ. Πράσινο χρώμα με διακεκομμένη γραμμή
        g2->SetLineColor(3); 
        g2->SetLineStyle(2); 
        g2->Draw("SAME");    }
    
    TCanvas *c5LS = new TCanvas("c5LS","Residuals LS  ", 100,0,800,800);
    c5LS->Divide(2,2);
    gausf->SetParameters(480, 0.01, 0.090, 50, 0.01, 0.3);
    c5LS->cd(1);Residual2L->Draw();Residual2L->Fit("gaus");//Residual1->Fit("gausf", "rq");
    c5LS->cd(2);Residual3L->Draw();Residual3L->Fit("gaus");//Residual2->Fit("gausf", "rq");
    c5LS->cd(3);Residual2S->Draw();Residual2S->Fit("gaus");//Residual3->Fit("gausf", "rq");
    c5LS->cd(4);Residual3S->Draw();Residual3S->Fit("gaus");//Residual4->Fit("gausf", "rq");


    c92->Print("Residual.pdf(","Title:page1");
    c94->Print("ResidualY.pdf","Title:page12");
    c5LS->Print("Residual.pdf)","Title:page3");

    
  
       gROOT->ForceStyle();

 
    
}

void mmanalys::Setstyle()
{
    {// Add my own options here:
        
        gROOT->Reset();
        gStyle->SetOptStat(111111);

        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetPadGridX(1);
        gStyle->SetPadGridY(1);

        gStyle->SetLabelSize  (0.55,"X");
        gStyle->SetLabelSize  (0.055,"Y");
        gStyle->SetLabelSize  (0.055,"Z");
        
        // Various
          gStyle->SetTickLength (0.03 ,"XYZ");
          gStyle->SetTitleSize  (0.055,"XYZ");
          gStyle->SetTitleFont  (42   ,"XYZ");
          gStyle->SetLabelFont  (42   ,"XYZ");
          gStyle->SetNdivisions (405  ,"XYZ");
          
          gStyle->SetTitleOffset(1.15,"x");
          gStyle->SetTitleOffset(1.25,"Y");
          gStyle->SetTitleOffset(1.15,"Z");
          
          gStyle->SetLabelOffset(0.004,"x");
          gStyle->SetLabelOffset(0.015,"Y");
          gStyle->SetLabelOffset(0.015,"Z");
          
          gStyle->SetLabelSize  (0.055,"x");
          gStyle->SetLabelSize  (0.055,"Y");
          gStyle->SetLabelSize  (0.055,"Z");
          
          gStyle->SetStatFont   (42);
          gStyle->SetTitleFont  (42);
          
        
        
        cout << "Styles are Set!" << endl;
        gROOT->ForceStyle();

        return; // atlasStyle;

    }
}
  

void mmanalys::lsm( double *x, double *y, int n , double *a, double *b)
{
    double sumx=0,sumy=0,sumxy=0,sumx2=0;
    double a1,b1;
    
    for(int i=0;i<=n-1;i++)
    {
        sumx=sumx +x[i];
        sumx2=sumx2 +x[i]*x[i];
        sumy=sumy +y[i];
        sumxy=sumxy +x[i]*y[i];
        
    }
    a1=((sumx2*sumy -sumx*sumxy)*1.0/(n*sumx2-sumx*sumx)*1.0);  // intercept
    b1=((n*sumxy-sumx*sumy)*1.0/(n*sumx2-sumx*sumx)*1.0);  //  slop
//    printf("\n\nThe line is Y=%3.3f +%3.3f X",a1,b1);
    *a=a1; *b=b1;
}





/*
       // Set font
        Int_t font = 42 ; // Arial->43

        atlasStyle->SetTextFont(font);
        atlasStyle->SetTextSize(tsize);

        atlasStyle->SetLabelFont(font, "x");
        atlasStyle->SetTitleFont(font, "x");
        atlasStyle->SetLabelFont(font, "y");
        atlasStyle->SetTitleFont(font, "y");
        atlasStyle->SetLabelFont(font, "z");
        atlasStyle->SetTitleFont(font, "z");

        atlasStyle->SetLabelSize(tsize, "x");
        atlasStyle->SetTitleSize(tsize, "x");
        atlasStyle->SetLabelSize(tsize, "y");
        atlasStyle->SetTitleSize(tsize, "y");
        atlasStyle->SetLabelSize(tsize, "z");
        atlasStyle->SetTitleSize(tsize, "z");

        atlasStyle->SetPadTickX(1);
        atlasStyle->SetPadTickY(1);
        // set title offsets (for axis label)
        atlasStyle->SetTitleXOffset(1.4);
        atlasStyle->SetTitleYOffset(1.4);
        atlasStyle->SetStatFontSize(0.04);
    //    atlasStyle->SetStatFont();

        
        // Use bold lines and markers
        atlasStyle->SetMarkerStyle(25);
        atlasStyle->SetMarkerSize(2.9);
        atlasStyle->SetHistLineWidth(2);
        atlasStyle->SetLineStyleString(2, "[12 12]");// # postscript dashes

        // Get rid of error bar caps
        atlasStyle->SetEndErrorSize(0.);

        // Do not display any of the standard histogram decorations
        atlasStyle->SetOptTitle(0);
        atlasStyle->SetOptStat(111111);
        atlasStyle->SetOptFit(11111);

        // Put tick marks on top and RHS of plots
        atlasStyle->SetPadTickX(1);
        atlasStyle->SetPadTickY(1);

        // Remove legend borders and set font
        atlasStyle->SetLegendBorderSize(0);
        atlasStyle->SetLegendFillColor(0);
        atlasStyle->SetLegendFont(43);
        atlasStyle->SetLegendTextSize(20);

        gROOT->SetStyle("atlasStyle");

*/
  

/*

void mmanalys::remove_holes(vector<MMHit> &a1, map<int, vector<MMHit> > &m, int hole_size){
    int hole =hole_size;
    int start=0;
    int stop=0;
    //   int nclu=0;
    int clu_id=0;
    vector<MMHit> cluster;
    
    if (a1.size()==1) {
        cluster.push_back(a1[0]);
        m.insert(pair<int, vector<MMHit> >(clu_id, cluster));
    }
    else // more that 1 hit
    {
        for (unsigned int s=0; s<(a1.size()-1); s++) {
            if ((a1[s+1].getx()-a1[s].getx())>hole+1) {
                stop=s;
                for (int i=start; i<=stop; i++) {
                    cluster.push_back(a1[i]);
                }
                //    for (unsigned int k=0; k<cluster.size(); k++) {    (cluster.at(k)).show();}
                m.insert(pair<int, vector<MMHit> >(clu_id, cluster));
                clu_id++;
                start=s+1;
                cluster.clear();
            }
        }
        stop=(a1.size()-1);
        for (int i=start; i<=stop; i++) {
            cluster.push_back(a1[i]);
            //      cout<<"Perase!!!!!"<<endl;
        }
        m.insert(pair<int, vector<MMHit> >(clu_id, cluster));
    }
}


float mmanalys::get_cluster_charge(vector<MMHit> &z)
{
    float qtot=0;
    for(unsigned int i=0; i<z.size(); i++) qtot+= (z.at(i)).getq();
    //    *clq=qtot;
    return qtot;
}

float mmanalys::get_cluster_length(vector<MMHit> &z )
{
    float cl_length=0;
    int size_cl = z.size();
    cl_length=(z.at(size_cl-1)).getx()-(z.at(0)).getx()+1;
    
    return cl_length;
}

float mmanalys::get_cluster_maxq(vector<MMHit> &z)
{
    float maxq=0;
    for(unsigned int i=0; i<z.size(); i++)
    if(maxq <= (z.at(i)).getq() ) maxq = (z.at(i)).getq() ;
    
    return maxq;
}

void mmanalys::cluster_pos(vector<MMHit> &z1m, float *pos, int chamber , float &pitch )  // returns the clusterpos in mm
{
    float xq=0, sq=0;
    float ppp= pitch ; // dummy, for not getting warning, z1m is allready in mm, no need of pitch
    if (chamber==1){
        for (unsigned int i=0; i<z1m.size(); i++){
            if((z1m.at(i)).getq()>0){
                xq+=(z1m.at(i)).getx()*(z1m.at(i)).getq();
                sq+=(z1m.at(i)).getq();
            }
        }
        if(sq!=0)  *pos=xq/sq;
    }
    else {  // chamber ==-1 doesn't make any sense, identical withl if
        for (unsigned int i=0; i<z1m.size(); i++){
            if((z1m.at(i)).getq()>0){
                xq+=(z1m.at(i).getx())*(z1m.at(i)).getq();
                //                xq+=(z1m.at(i).getx()+1)*(z1m.at(i)).getq();
                sq+=(z1m.at(i)).getq();
                //   x[i]=256-x[i]+1;
            }
        }
        if(sq!=0) *pos=xq/sq;
    }
    if (sq==0)     *pos= -300;
    //       cout << "pos= "<<*pos << "  q= "<< sq << endl;
    
    
}




void mmanalys::transform(vector<MMHit> &hitlist,int chamber,float pitch){   // transforms the strips x in mm
    float x_assist;
    unsigned int size1=hitlist.size();
    if (chamber==1){
        for (unsigned int i=0; i<size1; i++) {
            x_assist=(hitlist[i]).getx()*pitch;
            (hitlist[i]).setx(x_assist);
        }
    }
    else { //if (chamber==-1)
        for (unsigned int i=0; i<size1; i++) {
            x_assist=(256 -(hitlist[i]).getx()+1)*pitch;
            (hitlist[i]).setx(x_assist);
        }
    }
}


void mmanalys::lsm( float *x, float *y, int n , float *a, float *b)
{
    float sumx=0,sumy=0,sumxy=0,sumx2=0;
    double a1,b1;
    
    for(int i=0;i<=n-1;i++)
    {
        sumx=sumx +x[i];
        sumx2=sumx2 +x[i]*x[i];
        sumy=sumy +y[i];
        sumxy=sumxy +x[i]*y[i];
        
    }
    a1=((sumx2*sumy -sumx*sumxy)*1.0/(n*sumx2-sumx*sumx)*1.0);
    b1=((n*sumxy-sumx*sumy)*1.0/(n*sumx2-sumx*sumx)*1.0);
    //    printf("\n\nThe line is Y=%3.3f +%3.3f X",a1,b1);
    *a=a1; *b=b1;
}

void mmanalys::cluster_pos_utpc(vector<MMHit> &z,float drift_gap,float drift_vel,float pitch,float *pos,float *a1gr,float *a0gr,float *siga1,float *siga0,float *chi2)
// returns the clusterpos in mm using microtpc  (non aligned)
{
    int zs=z.size();
    float xt[400],yt[400],qt[400],yer[400], ex[400], ey[400];
    float aa0=0, aa1=0, aa0er=0, aa1er=0,  x2=0, Px2=0;
    int ndf=0;
    
    TF1 *fline = new TF1 ("fline","[0]+x*[1]",0,120);
    
    if(z.size()<=4) {
        *pos= -400;
        aa1=-1000;    aa0=-1000;
        aa1er=-1000;    aa0er=-1000;
        x2=-1000;        *chi2=-1000;
        *a1gr=aa1;  *a0gr=aa0; *siga1=aa1er;    *siga0=aa0er;    *chi2=x2/(zs-4);
    }
    else
    {
        for (int i=1; i<zs-1 ; i++)
        {
            xt[i-1]=(z.at(i)).getx();
            yt[i-1]=(z.at(i)).gety();
            qt[i-1]=(z.at(i)).getq();
            yer[i-1]=(z.at(i)).getyer();
            //    ey[i-1]=sqrt(pow(yer[i-1]*drift_vel,2)+pow(0.04*yt[i-1],2));
            ///    ex[i-1]=sqrt(pow(pitch/sqrt(12),2));//+ pow(0.08*yt[i-1],2));
            ey[i-1]=yer[i-1]*drift_vel+ 600/qt[i-1];  // Kordas
            //            ey[i-1]=yer[i-1]*drift_vel*35;
            ex[i-1]=2.1*pitch/sqrt(12);
        }
        //      TGraph *gr = new TGraph(zs-2,xt,yt);
        TGraphErrors *gr = new TGraphErrors(zs-2,xt,yt,ex,ey);
        //        if(evntnum == 113) gr->Fit("pol1"); else gr->Fit("pol1","Q");//,"Q"
        fline->SetParameters(0, xt[0]);
        if(xt[0]>xt[1])    fline->SetParameters(1, -1);
        else                fline->SetParameters(1, 1);
        gr->Fit("fline","Q");
        
        aa1        = gr->GetFunction("fline")->GetParameter(1);
        aa1er    = gr->GetFunction("fline")->GetParError(1);
        aa0        = gr->GetFunction("fline")->GetParameter(0);
        aa0er    = gr->GetFunction("fline")->GetParError(0);
        x2        = gr->GetFunction("fline")->GetChisquare();
        Px2        = gr->GetFunction("fline")->GetProb();
        ndf        = gr->GetFunction("fline")->GetNDF();
        
        
        if(aa1!=0 && Px2>0.002){
            *pos=((drift_gap/2-aa0)/aa1);
            for (int i=0; i<zs-2 ; i++) {
                x2+= (yt[i]- (aa0+ aa1*xt[i]))* (yt[i]- (aa0+ aa1*xt[i]));
            }
        }
        else {
            *pos= -400;
            x2=1000;
        }
        *a1gr=aa1;  *a0gr=aa0; *siga1=aa1er;  *siga0=aa0er;  *chi2=Px2;
        // actually is not x2 but Prob(x2)
    }   // end of else
    
}

void mmanalys::find_longest(map<int, vector<MMHit> > &m, vector<MMHit> &hitlist ){
    map<int, vector<MMHit> >::iterator it_fl;
    int max_length=0;
    int max_id=0;
    for(it_fl=m.begin(); it_fl!=m.end();++it_fl ){
        int size_cl = (it_fl->second).size();
        int length1=((it_fl->second).at(size_cl-1)).getx()-((it_fl->second).at(0)).getx()+1;
        
        if (length1>max_length) {
            max_length=length1;
            max_id=(it_fl->first);
        }
        
    }
    hitlist=(m.find(max_id)->second);
}

void mmanalys::find_heaviest(map<int, vector<MMHit> > &m, vector<MMHit> &hitlist ){
    map<int, vector<MMHit> >::iterator it_fl;
    float max_q=0;
    int max_id=0;
    float qsum=0;
    for(it_fl=m.begin(); it_fl!=m.end();++it_fl ){
        int size_cl = (it_fl->second).size();
        for (int j=0; j<size_cl; j++) {
            qsum=qsum+((it_fl->second).at(j)).getq();
        }
        if (qsum>max_q) {
            //           cout<<"q sum "<< qsum << "max_q "<< max_q<<endl;
            max_q=qsum;
            max_id=(it_fl->first);
        }
        qsum=0;
    }
    hitlist=(m.find(max_id)->second);
    //    cout<<"Heaviest"<<endl;
}



#include "TROOT.h" #include "TStyle.h" #include "TH1.h" #include "TCanvas.h"
void test ()
{
    gROOT->ForceStyle(kTRUE);
    gStyle->SetHistLineWidth(6);
    TH1D *h = new TH1D("h", "", 100, 0, 100);
    for (int i = 0; i < 100; i++) { h->SetBinContent(i+1, i); }
    TCanvas *c = new TCanvas("c");
    h->Draw();
    c->Modified();
    c->Update(); }
 */


 
