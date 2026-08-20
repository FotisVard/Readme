
void mmanalys::plot_Raw()

{
    
   // 	Setstyle();
    	gStyle->SetOptStat(111111);
    
    
 //   gStyle->SetOptStat(111111);

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



void mmanalys::plot_Clusters(){

    gROOT->ForceStyle(kTRUE);
    gStyle->SetOptStat(111111);
    gStyle->SetHistLineWidth(5);

    TCanvas *c10 = new TCanvas("c10","Cluster CHARGE X", 50,0,800,800);
        c10->Divide(2,2);
        c10->cd(1);cl_Charge_Tmm1X->Draw();
        c10->cd(2);cl_Charge_Tmm2X->Draw();
        c10->cd(3);cl_Charge_Tmm3X->Draw();
        c10->cd(4);cl_Charge_Tmm4X->Draw();
    TCanvas *c11 = new TCanvas("c11","Cluster CHARGE Y", 800,0,800,800);
        c11->Divide(2,2);
        c11->cd(1);cl_Charge_Tmm1Y->Draw();
        c11->cd(2);cl_Charge_Tmm2Y->Draw();
        c11->cd(3);cl_Charge_Tmm3Y->Draw();
        c11->cd(4);cl_Charge_Tmm4Y->Draw();

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

    TCanvas *c14 = new TCanvas("c14","Cluster width X ", 50,0,800,800);
        c14->Divide(2,2);
        c14->cd(1);cluster_length_Tmm1X->Draw();
        c14->cd(2);cluster_length_Tmm2X->Draw();
        c14->cd(3);cluster_length_Tmm3X->Draw();
        c14->cd(4);cluster_length_Tmm4X->Draw();
    TCanvas *c15 = new TCanvas("c15","Cluster width Y ", 900,0,800,800);
        c15->Divide(2,2);
        c15->cd(1);cluster_length_Tmm1Y->Draw();
        c15->cd(2);cluster_length_Tmm2Y->Draw();
        c15->cd(3);cluster_length_Tmm3Y->Draw();
        c15->cd(4);cluster_length_Tmm4Y->Draw();


    
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
