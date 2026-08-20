// Best 3-layer angleX, angleY, and angleXY histograms
TH1F *h_best3LayerAngleX = new TH1F("h_best3LayerAngleX", "Best 3-layer angleX;angleX (deg);events", 180, 40, 140);
TH1F *h_best3LayerAngleY = new TH1F("h_best3LayerAngleY", "Best 3-layer angleY;angleY (deg);events", 180, 40, 140);
TH2F *h_best3LayerAngleXY = new TH2F("h_best3LayerAngleXY", "Best 3-layer angleXY;angleX (deg);angleY (deg)", 180, 40, 140, 180, 40, 140);
TH1F *h_strps_Tmm1X= new TH1F("h_strps_Tmm1X", ";strip#",641,0,641);
TH1F *h_strps_Tmm2X= new TH1F("h_strps_Tmm2X", ";strip#",641,0,641);
TH1F *h_strps_Tmm3X= new TH1F("h_strps_Tmm3X", ";strip#",641,0,641);
TH1F *h_strps_Tmm4X= new TH1F("h_strps_Tmm4X", ";strip#",641,0,641);

TH1F *h_strps_Tmm1Y= new TH1F("h_strps_Tmm1Y", ";strip#",641,0,641);
TH1F *h_strps_Tmm2Y= new TH1F("h_strps_Tmm2Y", ";strip#",641,0,641);
TH1F *h_strps_Tmm3Y= new TH1F("h_strps_Tmm3Y", ";strip#",641,0,641);
TH1F *h_strps_Tmm4Y= new TH1F("h_strps_Tmm4Y", ";strip#",641,0,641);

TH1F *clu_per_eventTmm1X= new TH1F("clu_per_eventTmm1X", ";num_of_clu/ev", 10,0, 10);
TH1F *clu_per_eventTmm2X= new TH1F("clu_per_eventTmm2X", ";num_of_clu/ev", 10,0, 10);
TH1F *clu_per_eventTmm3X= new TH1F("clu_per_eventTmm3X", ";num_of_clu/ev", 10,0, 10);
TH1F *clu_per_eventTmm4X= new TH1F("clu_per_eventTmm4X", ";num_of_clu/ev", 10,0, 10);
TH1F *clu_per_eventTmm1Y= new TH1F("clu_per_eventTmm1Y", ";num_of_clu/ev", 10,0, 10);
TH1F *clu_per_eventTmm2Y= new TH1F("clu_per_eventTmm2Y", ";num_of_clu/ev", 10,0, 10);
TH1F *clu_per_eventTmm3Y= new TH1F("clu_per_eventTmm3Y", ";num_of_clu/ev", 10,0, 10);
TH1F *clu_per_eventTmm4Y= new TH1F("clu_per_eventTmm4Y", ";num_of_clu/ev", 10,0, 10);

TH1F *cl_Charge_Tmm1X= new TH1F("cl_Charge_Tmm1X", ";charge",300,0,5000);
TH1F *cl_Charge_Tmm2X= new TH1F("cl_Charge_Tmm2X", ";charge",300,0,5000);
TH1F *cl_Charge_Tmm3X= new TH1F("cl_Charge_Tmm3X", ";charge",300,0,5000);
TH1F *cl_Charge_Tmm4X= new TH1F("cl_Charge_Tmm4X", ";charge",300,0,5000);
TH1F *cl_Charge_Tmm1Y= new TH1F("cl_Charge_Tmm1Y", ";charge",300,0,9000);
TH1F *cl_Charge_Tmm2Y= new TH1F("cl_Charge_Tmm2Y", ";charge",300,0,9000);
TH1F *cl_Charge_Tmm3Y= new TH1F("cl_Charge_Tmm3Y", ";charge",300,0,9000);
TH1F *cl_Charge_Tmm4Y= new TH1F("cl_Charge_Tmm4Y", ";charge",300,0,9000);


TH1F *cluster_length_Tmm1X=new TH1F("cluster_length_Tmm1X", ";cluster length (strips)", 30,0,30);
TH1F *cluster_length_Tmm2X=new TH1F("cluster_length_Tmm2X", ";cluster length (strips)", 30,0,30);
TH1F *cluster_length_Tmm3X=new TH1F("cluster_length_Tmm3X", ";cluster length (strips)", 30,0,30);
TH1F *cluster_length_Tmm4X=new TH1F("cluster_length_Tmm4X", ";cluster length (strips)", 30,0,30);
TH1F *cluster_length_Tmm1Y=new TH1F("cluster_length_Tmm1Y", ";cluster length (strips)", 30,0,30);
TH1F *cluster_length_Tmm2Y=new TH1F("cluster_length_Tmm2Y", ";cluster length (strips)", 30,0,30);
TH1F *cluster_length_Tmm3Y=new TH1F("cluster_length_Tmm3Y", ";cluster length (strips)", 30,0,30);
TH1F *cluster_length_Tmm4Y=new TH1F("cluster_length_Tmm4Y", ";cluster length (strips)", 30,0,30);



TH1F *cl_withmaxq_Tmm1= new TH1F("cl_withmaxq_Tmm1", ";Cluster Charge",50,0,3000);
TH1F *cl_withmaxq_Tmm2= new TH1F("cl_withmaxq_Tmm2", ";Cluster Charge",50,0,3000);
TH1F *cl_withmaxq_Tmm3= new TH1F("cl_withmaxq_Tmm3", ";Cluster Charge",50,0,3000);
TH1F *cl_withmaxq_Tmm4= new TH1F("cl_withmaxq_Tmm4", ";Cluster Charge",50,0,3000);

TH1F *vertical_clq_X1= new TH1F("vertical_clq_X1", ";Cluster Charge",50,0,3000);
TH1F *vertical_clq_X2= new TH1F("vertical_clq_X2", ";Cluster Charge",50,0,3000);
TH1F *vertical_clq_X3= new TH1F("vertical_clq_X3", ";Cluster Charge",50,0,3000);
TH1F *vertical_clq_X4= new TH1F("vertical_clq_X4", ";Cluster Charge",50,0,3000);

TH1F *chargetmm1= new TH1F("chargetmm1",";",500,0,2500);
TH1F *chargetmm2= new TH1F("chargetmm2",";",500,0,2500);
TH1F *chargetmm3= new TH1F("chargetmm3",";",500,0,2500);
TH1F *chargetmm4= new TH1F("chargetmm4",";",500,0,2500);

TH1F *clustersize1= new TH1F("clustersize1",";",100,0,50);//lily
TH1F *clustersize2= new TH1F("clustersize2",";",100,0,35);//lily
TH1F *clustersize3= new TH1F("clustersize3",";",100,0,30);//lily
TH1F *clustersize4= new TH1F("clustersize4",";",100,0,30);//lily


TH1F *h_cl_evnt_Tmm1X= new TH1F("h_cl_evnt_Tmm1X", ";clusters/event",10,0,10);
TH1F *h_cl_evnt_Tmm2X= new TH1F("h_cl_evnt_Tmm2X", ";clusters/event",10,0,10);
TH1F *h_cl_evnt_Tmm3X= new TH1F("h_cl_evnt_Tmm3X", ";clusters/event",10,0,10);
TH1F *h_cl_evnt_Tmm4X= new TH1F("h_cl_evnt_Tmm4X", ";clusters/event",10,0,10);

TH1F *h_Event_typeX= new TH1F("h_Event_typeX", ";one cluster in LAYER ",20,0,20);
TH1F *h_Event_typeY= new TH1F("h_Event_typeY", ";one cluster in LAYER ",20,0,20);

TH2F *h_Event_typeXY= new TH2F("h_Event_typeXY", ";one cluster in X LAYER;one cluster in Y LAYER ",20,0,20,20,0,20);



TH1F *pos_Tmm1X= new TH1F("pos_Tmm1X", "; X (mm);",100, 0,100);
TH1F *pos_Tmm2X= new TH1F("pos_Tmm2X", "; X (mm);",100, 0,100);
TH1F *pos_Tmm3X= new TH1F("pos_Tmm3X", "; X (mm);",100, 0,100);
TH1F *pos_Tmm4X= new TH1F("pos_Tmm4X", "; X (mm);",100, 0,100);

TH1F *pos_Tmm1X_utpc= new TH1F("pos_Tmm1X_utpc", "; X (mm);",100, -50,150);
TH1F *pos_Tmm2X_utpc= new TH1F("pos_Tmm2X_utpc", "; X (mm);",100, -50,150);
TH1F *pos_Tmm3X_utpc= new TH1F("pos_Tmm3X_utpc", "; X (mm);",100, -50,150);
TH1F *pos_Tmm4X_utpc= new TH1F("pos_Tmm4X_utpc", "; X (mm);",100, -50,150);

TH1F *Residual1= new TH1F("Residual1", "; residual_1 (mm);",100, -12, 12); //-0.5,0.5);
TH1F *Residual2= new TH1F("Residual2", "; residual_2 (mm);",100, -12, 12); //-0.5,0.5);
TH1F *Residual3= new TH1F("Residual3", "; residual_3 (mm);",100, -12, 12); //-0.5,0.5);
TH1F *Residual4= new TH1F("Residual4", "; residual_4 (mm);",100, -12, 12); //-0.5,0.5);

TH1F *Residual1Y= new TH1F("Residual1Y", "; residual_1 (mm);",100, -12, 12); //-0.5,0.5);
TH1F *Residual2Y= new TH1F("Residual2Y", "; residual_2 (mm);",100, -12, 12); //-0.5,0.5);
TH1F *Residual3Y= new TH1F("Residual3Y", "; residual_3 (mm);",100, -12, 12); //-0.5,0.5);
TH1F *Residual4Y= new TH1F("Residual4Y", "; residual_4 (mm);",100, -12, 12); //-0.5,0.5);


TH1F *Residual1L= new TH1F("Residual1L", "; residual_1 (mm);",100, -1,1);
TH1F *Residual2L= new TH1F("Residual2L", "; residual_2 (mm);",100, -1,1);
TH1F *Residual3L= new TH1F("Residual3L", "; residual_3 (mm);",100, -1,1);
TH1F *Residual4L= new TH1F("Residual4L", "; residual_4 (mm);",100, -1,1);

TH1F *Residual1S= new TH1F("Residual1S", "; residual_1 (mm);",100, -1,1);
TH1F *Residual2S= new TH1F("Residual2S", "; residual_2 (mm);",100, -1,1);
TH1F *Residual3S= new TH1F("Residual3S", "; residual_3 (mm);",100, -1,1);
TH1F *Residual4S= new TH1F("Residual4S", "; residual_4 (mm);",100, -1,1);

TH1F *Resid_excl_1= new TH1F("Resid_excl_1", "; residuals1 excluding_1 (mm);",100, -1,1);
TH1F *Resid_excl_2= new TH1F("Resid_excl_2", "; residuals1 excluding_2 (mm);",100, -1,1);
TH1F *Resid_excl_3= new TH1F("Resid_excl_3", "; residuals1 excluding_3 (mm);",100, -1,1);
TH1F *Resid_excl_4= new TH1F("Resid_excl_4", "; residuals1 excluding_4 (mm);",100, -1,1);

TH1F *h_Rp123= new TH1F("h_Rp123", "; residual (mm);",100, -0.5,0.5);
TH1F *h_Rp124= new TH1F("h_Rp124", "; residual (mm);",100, -0.5,0.5);
TH1F *h_Rp234= new TH1F("h_Rp234", "; residual (mm);",100, -0.5,0.5);
TH1F *h_Rp134= new TH1F("h_Rp134", "; residual (mm);",100, -0.5,0.5);


TH1F *Align_1= new TH1F("Align_1", "; residuals 1 (mm);",100, -1,1);
TH1F *Align_2= new TH1F("Align_2", "; residuals 2 (mm);",100, -1,1);
TH1F *Align_3= new TH1F("Align_3", "; residuals 3 (mm);",100, -1,1);
TH1F *Align_4= new TH1F("Align_4", "; residuals 4 (mm);",100, -1,1);

TH1F *track_chi2= new TH1F("track_chi2",  ";  ;",10000, 0,20000);

TH1F *slopeX= new TH1F("slopeX", "; angle (deg);",180, 40,140);
TH1F *slopeY= new TH1F("slopeY", "; angle (deg);",180, 40,140);

TH2F *slopeXY= new TH2F("slopeXY", "; Xangle (deg);Yangle (deg)",180, 40,140, 180, 40,140);


TH1F *hslope4321= new TH1F("hslope4321", "; angle (deg);",180, 40,140);
TH1F *hslope432= new TH1F("hslope432", "; angle (deg);",180, 40,140);
TH1F *hslope43= new TH1F("hslope43", "; angle (deg);",180, 40,140);
TH1F *hslope321= new TH1F("hslope321", "; angle (deg);",180, 40,140);
TH1F *hslope32= new TH1F("hslope32", "; angle (deg);",180, 40,140);
TH1F *hslope21= new TH1F("hslope21", "; angle (deg);",180, 40,140);

TH1F *hslope4321Y= new TH1F("hslope4321Y", "; angle (deg);",180, 40,140);
TH1F *hslope432Y= new TH1F("hslope432Y", "; angle (deg);",180, 40,140);
TH1F *hslope43Y= new TH1F("hslope43Y", "; angle (deg);",180, 40,140);
TH1F *hslope321Y= new TH1F("hslope321Y", "; angle (deg);",180, 40,140);
TH1F *hslope32Y= new TH1F("hslope32Y", "; angle (deg);",180, 40,140);
TH1F *hslope21Y= new TH1F("hslope21Y", "; angle (deg);",180, 40,140);

TH1F *projection1 = new TH1F("projection1",    "; X (mm);",400, -400,400);
TH1F *projection2 = new TH1F("projection2",    "; X (mm);",400, -400,400);
TH1F *projection3 = new TH1F("projection3",    "; X (mm);",400, -400,400);
TH1F *projection4 = new TH1F("projection4",    "; X (mm);",400, -400,400);

TH1F *Delta_pos1 = new TH1F("Delta_pos1", ";Dx(mm)",500,-500, 500);
TH1F *Delta_pos2 = new TH1F("Delta_pos2", ";Dx(mm)",500,-500, 500);
TH1F *Delta_pos3 = new TH1F("Delta_pos3", ";Dx(mm)",500,-500, 500);
TH1F *Delta_pos4 = new TH1F("Delta_pos4", ";Dx(mm)",500,-500, 500);

TH1F *delta_angle1 = new TH1F("delta_angle1", ";Da(degrees)",100,-60, 60);
TH1F *delta_angle2 = new TH1F("delta_angle2", ";Da(degrees)",100,-60, 60);
TH1F *delta_angle3 = new TH1F("delta_angle3", ";Da(degrees)",100,-60, 60);
TH1F *delta_angle4 = new TH1F("delta_angle4", ";Da(degrees)",100,-60, 60);


TH2F *Delta_pos_angle1 = new TH2F("Delta_pos_angle1", ";Dx(mm);angle", 500,-500, 500, 100, 0,180);
TH2F *Delta_pos_angle2 = new TH2F("Delta_pos_angle2", ";Dx(mm);angle", 500,-500, 500, 100, 0,180);
TH2F *Delta_pos_angle3 = new TH2F("Delta_pos_angle3", ";Dx(mm);angle", 500,-500, 500, 100, 0,180);
TH2F *Delta_pos_angle4 = new TH2F("Delta_pos_angle4", ";Dx(mm);angle", 500,-500, 500, 100, 0,180);

TH2F *slope_tpc_ang1 = new TH2F("slope_tpc_ang1", ";Track Slope(mm);tpc_angle", 100,0, 180, 100, 0,180);
TH2F *slope_tpc_ang2 = new TH2F("slope_tpc_ang2", ";Track Slope(mm);tpc_angle", 100,0, 180, 100, 0,180);
TH2F *slope_tpc_ang3 = new TH2F("slope_tpc_ang3", ";Track Slope(mm);tpc_angle", 100,0, 180, 100, 0,180);
TH2F *slope_tpc_ang4 = new TH2F("slope_tpc_ang4", ";Track Slope(mm);tpc_angle", 100,0, 180, 100, 0,180);

TH2F *centroid_utpc1 = new TH2F("centroid_utpc1", ";Centroid-pos(mm);tpc-pos(mm)", 100,0, 100, 100, 0,100);
TH2F *centroid_utpc2 = new TH2F("centroid_utpc2", ";Centroid-pos(mm);tpc-pos(mm)", 100,0, 100, 100, 0,100);
TH2F *centroid_utpc3 = new TH2F("centroid_utpc3", ";Centroid-pos(mm);tpc-pos(mm)", 100,0, 100, 100, 0,100);
TH2F *centroid_utpc4 = new TH2F("centroid_utpc4", ";Centroid-pos(mm);tpc-pos(mm)", 100,0, 100, 100, 0,100);


TH2F *charge_per_strip_tmm1 = new TH2F("charge_per_strip_tmm1",";strip;charge",640,0,640,2500,-10,2500 );//lily
TH2F *charge_per_strip_tmm2 = new TH2F("charge_per_strip_tmm2",";strip;charge",640,0,640,2500,-10,2500 );//lily
TH2F *charge_per_strip_tmm3 = new TH2F("charge_per_strip_tmm3",";strip;charge",640,0,640,2500,-10,2500 );//lily
TH2F *charge_per_strip_tmm4 = new TH2F("charge_per_strip_tmm4",";strip;charge",640,0,640,2500,-10,2500 );//lily

TH1F *timeT1 = new TH1F("timeT1", ";time(ns)",160,1, 500);
TH1F *timeT2 = new TH1F("timeT2", ";time(ns)",160,1, 500);
TH1F *timeT3 = new TH1F("timeT3", ";time(ns)",160,1, 500);
TH1F *timeT4 = new TH1F("timeT4", ";time(ns)",160,1, 500);

TF1 *trackfit = new TF1("trackfit", "pol1", -100, 1100);    // One way of fitting
TF1 *gausfit = new TF1("gausfit", "gaus", -0.1, 0.1);    // One way of fitting
TF1 *gausfit2 = new TF1("gausfit2", "gaus", -0.25, 0.25);    // One way of fitting
TF1 *gausf = new TF1("gausf", "gaus+gaus", -1.9, 1.9);    // One way of fitting
TF1 *g1 = new TF1("g1", "gaus", 0, 180);    // for draw
TF1 *g2 = new TF1("g2", "gaus", 0, 180);    // for draw

TF1 *gaus2 = new TF1("gaus2", "gaus(0)+gaus(3)", 0, 180 );    // One way of fitting

TH1F *utpc_slop1 = new TH1F("utpc_slop1", ";slope(rad)",200,0,180);
TH1F *utpc_slop2 = new TH1F("utpc_slop2", ";slope(rad)",200,0,180);
TH1F *utpc_slop3 = new TH1F("utpc_slop3", ";slope(rad)",200,0,180);
TH1F *utpc_slop4 = new TH1F("utpc_slop4", ";slope(rad)",200,0,180);

TH2F *h1_angle_vs_q_234 = new TH2F("h1_angle_vs_q_234", "Cluster Charge; angle(deg)", 180, -90,90, 2200, 0, 2200 );
TH2F *h2_angle_vs_q_234 = new TH2F("h2_angle_vs_q_234", "Cluster Charge; angle(deg)", 180, -90,90, 2200, 0, 2200 );
TH2F *h3_angle_vs_q_234 = new TH2F("h3_angle_vs_q_234", "Cluster Charge; angle(deg)", 180, -90,90, 2200, 0, 2200 );
TH2F *h4_angle_vs_q_234 = new TH2F("h4_angle_vs_q_234", "Cluster Charge; angle(deg)", 180, -90,90, 2200, 0, 2200 );

//F1 *g1 = new TF1("g1","gaus",0,1023);
//TF1 *g2 = new TF1("g2","gaus",0,1023);
//TF1 *g3 = new TF1("g3","gaus",0,1023);
// total functon for fitting
//TF1 *fitfun = new TF1("fitfun","pol1(0)+gaus(2)+gaus(5)+gaus(8)",0,1023);

TH2F *h_map_Tmm1_XY= new TH2F("h_map_Tmm1_XY", ";mm;mm",50,0,100, 50,0,100);
TH2F *h_map_Tmm2_XY= new TH2F("h_map_Tmm2_XY", ";mm;mm",50,0,100, 50,0,100);
TH2F *h_map_Tmm3_XY= new TH2F("h_map_Tmm3_XY", ";mm;mm",50,0,100, 50,0,100);
TH2F *h_map_Tmm4_XY= new TH2F("h_map_Tmm4_XY", ";mm;mm",50,0,100, 50,0,100);

TH2F *h_map_XY= new TH2F("h_map_XY", ";mm;mm",400,-100,700, 400,-100,700);


// Add near other histogram declarations:
TH2F *h_clus_angle_X_cl0 = new TH2F("h_clus_angle_X_cl0", "Cluster Map vs Angle X (nclus>0);Angle (deg);Total Clusters", 100, 60, 120, 40, 0, 40);
TH2F *h_clus_angle_X_cl1 = new TH2F("h_clus_angle_X_cl1", "Cluster Map vs Angle X (nclus>1);Angle (deg);Total Clusters", 100, 60, 120, 40, 0, 40);
TH2F *h_clus_angle_X_cl2 = new TH2F("h_clus_angle_X_cl2", "Cluster Map vs Angle X (nclus>2);Angle (deg);Total Clusters", 100, 60, 120, 40, 0, 40);

TH2F *h_clus_angle_Y_cl0 = new TH2F("h_clus_angle_Y_cl0", "Cluster Map vs Angle Y (nclus>0);Angle (deg);Total Clusters", 100, 60, 120, 40, 0, 40);
TH2F *h_clus_angle_Y_cl1 = new TH2F("h_clus_angle_Y_cl1", "Cluster Map vs Angle Y (nclus>1);Angle (deg);Total Clusters", 100, 60, 120, 40, 0, 40);
TH2F *h_clus_angle_Y_cl2 = new TH2F("h_clus_angle_Y_cl2", "Cluster Map vs Angle Y (nclus>2);Angle (deg);Total Clusters", 100, 60, 120, 40, 0, 40);

TH2F *h_charge_angle_X_cl0 = new TH2F("h_charge_angle_X_cl0", "Charge Map vs Angle X (nclus>0);Angle (deg);Total Charge", 100, 60, 120, 200, 0, 20000);
TH2F *h_charge_angle_X_cl1 = new TH2F("h_charge_angle_X_cl1", "Charge Map vs Angle X (nclus>1);Angle (deg);Total Charge", 100, 60, 120, 200, 0, 20000);
TH2F *h_charge_angle_X_cl2 = new TH2F("h_charge_angle_X_cl2", "Charge Map vs Angle X (nclus>2);Angle (deg);Total Charge", 100, 60, 120, 200, 0, 20000);

TH2F *h_charge_angle_Y_cl0 = new TH2F("h_charge_angle_Y_cl0", "Charge Map vs Angle Y (nclus>0);Angle (deg);Total Charge", 100, 60, 120, 200, 0, 20000);
TH2F *h_charge_angle_Y_cl1 = new TH2F("h_charge_angle_Y_cl1", "Charge Map vs Angle Y (nclus>1);Angle (deg);Total Charge", 100, 60, 120, 200, 0, 20000);
TH2F *h_charge_angle_Y_cl2 = new TH2F("h_charge_angle_Y_cl2", "Charge Map vs Angle Y (nclus>2);Angle (deg);Total Charge", 100, 60, 120, 200, 0, 20000);




// --- Leading Cluster Charge (Ανά Επίπεδο) ---
TH1F *h_lead0_charge_X1 = new TH1F("h_lead0_charge_X1", "Leading Cluster Charge X1;Charge;Counts", 200, 0, 5000);
TH1F *h_lead0_charge_X2 = new TH1F("h_lead0_charge_X2", "Leading Cluster Charge X2;Charge;Counts", 200, 0, 5000);
TH1F *h_lead0_charge_X3 = new TH1F("h_lead0_charge_X3", "Leading Cluster Charge X3;Charge;Counts", 200, 0, 5000);
TH1F *h_lead0_charge_X4 = new TH1F("h_lead0_charge_X4", "Leading Cluster Charge X4;Charge;Counts", 200, 0, 5000);

TH1F *h_lead0_charge_Y1 = new TH1F("h_lead0_charge_Y1", "Leading Cluster Charge Y1;Charge;Counts", 200, 0, 5000);
TH1F *h_lead0_charge_Y2 = new TH1F("h_lead0_charge_Y2", "Leading Cluster Charge Y2;Charge;Counts", 200, 0, 5000);
TH1F *h_lead0_charge_Y3 = new TH1F("h_lead0_charge_Y3", "Leading Cluster Charge Y3;Charge;Counts", 200, 0, 5000);
TH1F *h_lead0_charge_Y4 = new TH1F("h_lead0_charge_Y4", "Leading Cluster Charge Y4;Charge;Counts", 200, 0, 5000);

// --- Leading Cluster Strips (Ανά Επίπεδο - 640 Strips) ---
TH1F *h_lead0_strips_X1 = new TH1F("h_lead0_strips_X1", "0 Leading Cluster Strips X1;Strip#;Counts", 640, 0, 640);
TH1F *h_lead0_strips_X2 = new TH1F("h_lead0_strips_X2", "0 Leading Cluster Strips X2;Strip#;Counts", 640, 0, 640);
TH1F *h_lead0_strips_X3 = new TH1F("h_lead0_strips_X3", "0 Leading Cluster Strips X3;Strip#;Counts", 640, 0, 640);
TH1F *h_lead0_strips_X4 = new TH1F("h_lead0_strips_X4", "0 Leading Cluster Strips X4;Strip#;Counts", 640, 0, 640);

TH1F *h_lead0_strips_Y1 = new TH1F("h_lead0_strips_Y1", "0 Leading Cluster Strips Y1;Strip#;Counts", 640, 0, 640);
TH1F *h_lead0_strips_Y2 = new TH1F("h_lead0_strips_Y2", "0 Leading Cluster Strips Y2;Strip#;Counts", 640, 0, 640);
TH1F *h_lead0_strips_Y3 = new TH1F("h_lead0_strips_Y3", "0 Leading Cluster Strips Y3;Strip#;Counts", 640, 0, 640);
TH1F *h_lead0_strips_Y4 = new TH1F("h_lead0_strips_Y4", "0 Leading Cluster Strips Y4;Strip#;Counts", 640, 0, 640);

// --- Leading Cluster Strips (Ανά Επίπεδο - 640 Strips) ---
TH1F *h_lead1_strips_X1 = new TH1F("h_lead1_strips_X1", "1 Leading Cluster Strips X1;Strip#;Counts", 640, 0, 640);
TH1F *h_lead1_strips_X2 = new TH1F("h_lead1_strips_X2", "1 Leading Cluster Strips X2;Strip#;Counts", 640, 0, 640);
TH1F *h_lead1_strips_X3 = new TH1F("h_lead1_strips_X3", "1 Leading Cluster Strips X3;Strip#;Counts", 640, 0, 640);
TH1F *h_lead1_strips_X4 = new TH1F("h_lead1_strips_X4", "1 Leading Cluster Strips X4;Strip#;Counts", 640, 0, 640);

TH1F *h_lead1_strips_Y1 = new TH1F("h_lead1_strips_Y1", "1 Leading Cluster Strips Y1;Strip#;Counts", 640, 0, 640);
TH1F *h_lead1_strips_Y2 = new TH1F("h_lead1_strips_Y2", "1 Leading Cluster Strips Y2;Strip#;Counts", 640, 0, 640);
TH1F *h_lead1_strips_Y3 = new TH1F("h_lead1_strips_Y3", "1 Leading Cluster Strips Y3;Strip#;Counts", 640, 0, 640);
TH1F *h_lead1_strips_Y4 = new TH1F("h_lead1_strips_Y4", "1 Leading Cluster Strips Y4;Strip#;Counts", 640, 0, 640);

// --- Leading Cluster Charge (Ανά Επίπεδο) ---
TH1F *h_lead1_charge_X1 = new TH1F("h_lead1_charge_X1", "Leading Cluster Charge X1;Charge;Counts", 200, 0, 5000);
TH1F *h_lead1_charge_X2 = new TH1F("h_lead1_charge_X2", "Leading Cluster Charge X2;Charge;Counts", 200, 0, 5000);
TH1F *h_lead1_charge_X3 = new TH1F("h_lead1_charge_X3", "Leading Cluster Charge X3;Charge;Counts", 200, 0, 5000);
TH1F *h_lead1_charge_X4 = new TH1F("h_lead1_charge_X4", "Leading Cluster Charge X4;Charge;Counts", 200, 0, 5000);

TH1F *h_lead1_charge_Y1 = new TH1F("h_lead1_charge_Y1", "Leading Cluster Charge Y1;Charge;Counts", 200, 0, 5000);
TH1F *h_lead1_charge_Y2 = new TH1F("h_lead1_charge_Y2", "Leading Cluster Charge Y2;Charge;Counts", 200, 0, 5000);
TH1F *h_lead1_charge_Y3 = new TH1F("h_lead1_charge_Y3", "Leading Cluster Charge Y3;Charge;Counts", 200, 0, 5000);
TH1F *h_lead1_charge_Y4 = new TH1F("h_lead1_charge_Y4", "Leading Cluster Charge Y4;Charge;Counts", 200, 0, 5000);

// --- Leading Cluster Strips (Ανά Επίπεδο - 640 Strips) ---
TH1F *h_lead2_strips_X1 = new TH1F("h_lead2_strips_X1", "2 Leading Cluster Strips X1;Strip#;Counts", 640, 0, 640);
TH1F *h_lead2_strips_X2 = new TH1F("h_lead2_strips_X2", "2 Leading Cluster Strips X2;Strip#;Counts", 640, 0, 640);
TH1F *h_lead2_strips_X3 = new TH1F("h_lead2_strips_X3", "2 Leading Cluster Strips X3;Strip#;Counts", 640, 0, 640);
TH1F *h_lead2_strips_X4 = new TH1F("h_lead2_strips_X4", "2 Leading Cluster Strips X4;Strip#;Counts", 640, 0, 640);

TH1F *h_lead2_strips_Y1 = new TH1F("h_lead2_strips_Y1", "2 Leading Cluster Strips Y1;Strip#;Counts", 640, 0, 640);
TH1F *h_lead2_strips_Y2 = new TH1F("h_lead2_strips_Y2", "2 Leading Cluster Strips Y2;Strip#;Counts", 640, 0, 640);
TH1F *h_lead2_strips_Y3 = new TH1F("h_lead2_strips_Y3", "2 Leading Cluster Strips Y3;Strip#;Counts", 640, 0, 640);
TH1F *h_lead2_strips_Y4 = new TH1F("h_lead2_strips_Y4", "2 Leading Cluster Strips Y4;Strip#;Counts", 640, 0, 640);

// --- Leading Cluster Charge (Ανά Επίπεδο) ---
TH1F *h_lead2_charge_X1 = new TH1F("h_lead2_charge_X1", "2 Leading Cluster Charge X1;Charge;Counts", 200, 0, 5000);
TH1F *h_lead2_charge_X2 = new TH1F("h_lead2_charge_X2", "2 Leading Cluster Charge X2;Charge;Counts", 200, 0, 5000);
TH1F *h_lead2_charge_X3 = new TH1F("h_lead2_charge_X3", "2 Leading Cluster Charge X3;Charge;Counts", 200, 0, 5000);
TH1F *h_lead2_charge_X4 = new TH1F("h_lead2_charge_X4", "2 Leading Cluster Charge X4;Charge;Counts", 200, 0, 5000);

TH1F *h_lead2_charge_Y1 = new TH1F("h_lead2_charge_Y1", "2 Leading Cluster Charge Y1;Charge;Counts", 200, 0, 5000);
TH1F *h_lead2_charge_Y2 = new TH1F("h_lead2_charge_Y2", "2 Leading Cluster Charge Y2;Charge;Counts", 200, 0, 5000);
TH1F *h_lead2_charge_Y3 = new TH1F("h_lead2_charge_Y3", "2 Leading Cluster Charge Y3;Charge;Counts", 200, 0, 5000);
TH1F *h_lead2_charge_Y4 = new TH1F("h_lead2_charge_Y4", "2 Leading Cluster Charge Y4;Charge;Counts", 200, 0, 5000);



// --- Διαγράμματα Strip Entries ανά Γωνία (80-88, 88-92, 92-100) ---
TH1F *h_strips_ang_80_88  = new TH1F("h_strips_ang_80_88",  "Strip Entries (80-88 deg);Strip#;Entries", 640, 0, 640);
TH1F *h_strips_ang_88_92  = new TH1F("h_strips_ang_88_92",  "Strip Entries (88-92 deg);Strip#;Entries", 640, 0, 640);
TH1F *h_strips_ang_92_100 = new TH1F("h_strips_ang_92_100", "Strip Entries (92-100 deg);Strip#;Entries", 640, 0, 640);



// --- 2D strip correlations Xn vs Yn in selected angle windows ---
TH2F *h_stripXY_60_80_X1Y1 = new TH2F("h_stripXY_60_80_X1Y1", "X1 vs Y1 (60-80 deg);X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_60_80_X2Y2 = new TH2F("h_stripXY_60_80_X2Y2", "X2 vs Y2 (60-80 deg);X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_60_80_X3Y3 = new TH2F("h_stripXY_60_80_X3Y3", "X3 vs Y3 (60-80 deg);X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_60_80_X4Y4 = new TH2F("h_stripXY_60_80_X4Y4", "X4 vs Y4 (60-80 deg);X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);

TH2F *h_stripXY_60_80_X1Y1_cl1 = new TH2F("h_stripXY_60_80_X1Y1_cl1", "X1 vs Y1 (60-80 deg, nclus>1);X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_60_80_X2Y2_cl1 = new TH2F("h_stripXY_60_80_X2Y2_cl1", "X2 vs Y2 (60-80 deg, nclus>1);X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_60_80_X3Y3_cl1 = new TH2F("h_stripXY_60_80_X3Y3_cl1", "X3 vs Y3 (60-80 deg, nclus>1);X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_60_80_X4Y4_cl1 = new TH2F("h_stripXY_60_80_X4Y4_cl1", "X4 vs Y4 (60-80 deg, nclus>1);X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);

TH2F *h_stripXY_60_80_X1Y1_cl2 = new TH2F("h_stripXY_60_80_X1Y1_cl2", "X1 vs Y1 (60-80 deg, nclus>2);X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_60_80_X2Y2_cl2 = new TH2F("h_stripXY_60_80_X2Y2_cl2", "X2 vs Y2 (60-80 deg, nclus>2);X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_60_80_X3Y3_cl2 = new TH2F("h_stripXY_60_80_X3Y3_cl2", "X3 vs Y3 (60-80 deg, nclus>2);X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_60_80_X4Y4_cl2 = new TH2F("h_stripXY_60_80_X4Y4_cl2", "X4 vs Y4 (60-80 deg, nclus>2);X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);

TH2F *h_stripXY_85_95_X1Y1 = new TH2F("h_stripXY_85_95_X1Y1", "X1 vs Y1 (85-95 deg);X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_85_95_X2Y2 = new TH2F("h_stripXY_85_95_X2Y2", "X2 vs Y2 (85-95 deg);X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_85_95_X3Y3 = new TH2F("h_stripXY_85_95_X3Y3", "X3 vs Y3 (85-95 deg);X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_85_95_X4Y4 = new TH2F("h_stripXY_85_95_X4Y4", "X4 vs Y4 (85-95 deg);X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);

TH2F *h_stripXY_85_95_X1Y1_cl1 = new TH2F("h_stripXY_85_95_X1Y1_cl1", "X1 vs Y1 (85-95 deg, nclus>1);X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_85_95_X2Y2_cl1 = new TH2F("h_stripXY_85_95_X2Y2_cl1", "X2 vs Y2 (85-95 deg, nclus>1);X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_85_95_X3Y3_cl1 = new TH2F("h_stripXY_85_95_X3Y3_cl1", "X3 vs Y3 (85-95 deg, nclus>1);X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_85_95_X4Y4_cl1 = new TH2F("h_stripXY_85_95_X4Y4_cl1", "X4 vs Y4 (85-95 deg, nclus>1);X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);

TH2F *h_stripXY_85_95_X1Y1_cl2 = new TH2F("h_stripXY_85_95_X1Y1_cl2", "X1 vs Y1 (85-95 deg, nclus>2);X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_85_95_X2Y2_cl2 = new TH2F("h_stripXY_85_95_X2Y2_cl2", "X2 vs Y2 (85-95 deg, nclus>2);X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_85_95_X3Y3_cl2 = new TH2F("h_stripXY_85_95_X3Y3_cl2", "X3 vs Y3 (85-95 deg, nclus>2);X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_85_95_X4Y4_cl2 = new TH2F("h_stripXY_85_95_X4Y4_cl2", "X4 vs Y4 (85-95 deg, nclus>2);X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);

TH2F *h_stripXY_100_120_X1Y1 = new TH2F("h_stripXY_100_120_X1Y1", "X1 vs Y1 (100-120 deg);X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_100_120_X2Y2 = new TH2F("h_stripXY_100_120_X2Y2", "X2 vs Y2 (100-120 deg);X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_100_120_X3Y3 = new TH2F("h_stripXY_100_120_X3Y3", "X3 vs Y3 (100-120 deg);X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_100_120_X4Y4 = new TH2F("h_stripXY_100_120_X4Y4", "X4 vs Y4 (100-120 deg);X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);

TH2F *h_stripXY_100_120_X1Y1_cl1 = new TH2F("h_stripXY_100_120_X1Y1_cl1", "X1 vs Y1 (100-120 deg, nclus>1);X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_100_120_X2Y2_cl1 = new TH2F("h_stripXY_100_120_X2Y2_cl1", "X2 vs Y2 (100-120 deg, nclus>1);X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_100_120_X3Y3_cl1 = new TH2F("h_stripXY_100_120_X3Y3_cl1", "X3 vs Y3 (100-120 deg, nclus>1);X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_100_120_X4Y4_cl1 = new TH2F("h_stripXY_100_120_X4Y4_cl1", "X4 vs Y4 (100-120 deg, nclus>1);X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);

TH2F *h_stripXY_100_120_X1Y1_cl2 = new TH2F("h_stripXY_100_120_X1Y1_cl2", "X1 vs Y1 (100-120 deg, nclus>2);X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_100_120_X2Y2_cl2 = new TH2F("h_stripXY_100_120_X2Y2_cl2", "X2 vs Y2 (100-120 deg, nclus>2);X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_100_120_X3Y3_cl2 = new TH2F("h_stripXY_100_120_X3Y3_cl2", "X3 vs Y3 (100-120 deg, nclus>2);X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_100_120_X4Y4_cl2 = new TH2F("h_stripXY_100_120_X4Y4_cl2", "X4 vs Y4 (100-120 deg, nclus>2);X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);


// --- Leading Clusters from 0-3 ---
// --- Leading Clusters 1X ---
TH2F *h_charge_strip_1X1 = new TH2F("h_charge_strip_1X1", "1 Cluster X1: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_1X2 = new TH2F("h_charge_strip_1X2", "1 Cluster X2: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_1X3 = new TH2F("h_charge_strip_1X3", "1 Cluster X3: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_1X4 = new TH2F("h_charge_strip_1X4", "1 Cluster X4: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);


// --- Leading Clusters 1Y ---
TH2F *h_charge_strip_1Y1 = new TH2F("h_charge_strip_1Y1", "1 Cluster Y1: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_1Y2 = new TH2F("h_charge_strip_1Y2", "1 Cluster Y2: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_1Y3 = new TH2F("h_charge_strip_1Y3", "1 Cluster Y3: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_1Y4 = new TH2F("h_charge_strip_1Y4", "1 Cluster Y4: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);   

// --- Leading Clusters 2X ---
TH2F *h_charge_strip_2X1 = new TH2F("h_charge_strip_2X1", "2 Cluster X1: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_2X2 = new TH2F("h_charge_strip_2X2", "2 Cluster X2: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_2X3 = new TH2F("h_charge_strip_2X3", "2 Cluster X3: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_2X4 = new TH2F("h_charge_strip_2X4", "2 Cluster X4: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);

// --- Leading Clusters 2Y ---
TH2F *h_charge_strip_2Y1 = new TH2F("h_charge_strip_2Y1", "2 Cluster Y1: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_2Y2 = new TH2F("h_charge_strip_2Y2", "2 Cluster Y2: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_2Y3 = new TH2F("h_charge_strip_2Y3", "2 Cluster Y3: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_2Y4 = new TH2F("h_charge_strip_2Y4", "2 Cluster Y4: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);  

// --- Leading Clusters 3X ---
TH2F *h_charge_strip_3X1 = new TH2F("h_charge_strip_3X1", "3 Cluster X1: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_3X2 = new TH2F("h_charge_strip_3X2", "3 Cluster X2: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_3X3 = new TH2F("h_charge_strip_3X3", "3 Cluster X3: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_3X4 = new TH2F("h_charge_strip_3X4", "3 Cluster X4: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);

// --- Leading Clusters 3Y ---
TH2F *h_charge_strip_3Y1 = new TH2F("h_charge_strip_3Y1", "3 Cluster Y1: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_3Y2 = new TH2F("h_charge_strip_3Y2", "3 Cluster Y2: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_3Y3 = new TH2F("h_charge_strip_3Y3", "3 Cluster Y3: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);
TH2F *h_charge_strip_3Y4 = new TH2F("h_charge_strip_3Y4", "3 Cluster Y4: Charge per Strip;Strip#;Total Charge", 640, 0, 640 , 500, 0, 2000);  
//--- End of Leading Clusters---


// --- Γραφήματα Γωνίας X (TH1F με βάρος το Size/Charge) ---

TH2F *h_cluSize_vs_AngleX_gt0 = new TH2F("h_cluSize_vs_AngleX_gt0","Size vs Angle X (nclus>0);Angle (deg);Total Size",100, 60, 120 , 100 , 0, 100);
TH2F *h_cluSize_vs_AngleX_gt1 = new TH2F("h_cluSize_vs_AngleX_gt1","Size vs Angle X (nclus>1);Angle (deg);Total Size",100, 60, 120 , 100 , 0, 100);
TH2F *h_cluSize_vs_AngleX_gt2 = new TH2F("h_cluSize_vs_AngleX_gt2","Size vs Angle X (nclus>2);Angle (deg);Total Size",100, 60, 120 , 100, 0, 100);
 
// per-layer size-vs-angle histograms (X1..X4)
TH2F *h_cluSize_vs_AngleX1_gt0 = new TH2F("h_cluSize_vs_AngleX1_gt0","Size 1 vs Angle X1 Angle (deg);Cluster 1 Count (X1)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleX2_gt0 = new TH2F("h_cluSize_vs_AngleX2_gt0","Size 1 vs Angle X2 Angle (deg);Cluster 1 Count (X2)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleX3_gt0 = new TH2F("h_cluSize_vs_AngleX3_gt0","Size 1 vs Angle X3 Angle (deg);Cluster 1 Count (X3)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleX4_gt0 = new TH2F("h_cluSize_vs_AngleX4_gt0","Size 1 vs Angle X4 Angle (deg);Cluster 1 Count (X4)",100,60,120,20,0,20);

// per-layer size-vs-angle histograms (X1..X4)
TH2F *h_cluSize_vs_AngleX1_gt1 = new TH2F("h_cluSize_vs_AngleX1_gt1","Size 2 vs Angle X1 Angle (deg);Cluster 2 Count (X1)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleX2_gt1 = new TH2F("h_cluSize_vs_AngleX2_gt1","Size 2 vs Angle X2 Angle (deg);Cluster 2 Count (X2)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleX3_gt1 = new TH2F("h_cluSize_vs_AngleX3_gt1","Size 2 vs Angle X3 Angle (deg);Cluster 2 Count (X3)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleX4_gt1 = new TH2F("h_cluSize_vs_AngleX4_gt1","Size 2 vs Angle X4 Angle (deg);Cluster 2 Count (X4)",100,60,120,20,0,20);

// per-layer size-vs-angle histograms (X1..X4)
TH2F *h_cluSize_vs_AngleX1_gt2 = new TH2F("h_cluSize_vs_AngleX1_gt2","Size 3 vs Angle X1 Angle (deg);Cluster 3 Count (X1)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleX2_gt2 = new TH2F("h_cluSize_vs_AngleX2_gt2","Size 3 vs Angle X2 Angle (deg);Cluster 3 Count (X2)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleX3_gt2 = new TH2F("h_cluSize_vs_AngleX3_gt2","Size 3 vs Angle X3 Angle (deg);Cluster 3 Count (X3)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleX4_gt2 = new TH2F("h_cluSize_vs_AngleX4_gt2","Size 3 vs Angle X4 Angle (deg);Cluster 3 Count (X4)",100,60,120,20,0,20);


TH2F *h_cluCharge_vs_AngleX_gt0 = new TH2F("h_cluCharge_vs_AngleX_gt0","Charge vs Angle X (nclus>0);Angle (deg);Total Charge",100,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleX_gt1 = new TH2F("h_cluCharge_vs_AngleX_gt1","Charge vs Angle X (nclus>4);Angle (deg);Total Charge",100,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleX_gt2 = new TH2F("h_cluCharge_vs_AngleX_gt2","Charge vs Angle X (nclus>8);Angle (deg);Total Charge",100,60,120,600,0,6000);

// --- Γραφήματα Γωνίας Y (TH1F με βάρος το Size/Charge) ---
TH2F *h_cluSize_vs_AngleY_gt0 = new TH2F("h_cluSize_vs_AngleY_gt0","Size vs Angle Y (nclus>0);Angle (deg);Total Size",100, 60, 120, 100 , 0, 100);
TH2F *h_cluSize_vs_AngleY_gt1 = new TH2F("h_cluSize_vs_AngleY_gt1","Size vs Angle Y (nclus>1);Angle (deg);Total Size",100, 60, 120, 100 , 0, 100);
TH2F *h_cluSize_vs_AngleY_gt2 = new TH2F("h_cluSize_vs_AngleY_gt2","Size vs Angle Y (nclus>2);Angle (deg);Total Size",100, 60, 120, 100 , 0, 100);

// per-layer size-vs-angle histograms (Y1..Y4)
TH2F *h_cluSize_vs_AngleY1_gt0 = new TH2F("h_cluSize_vs_AngleY1_gt0","Size 1 vs Angle Y1 Angle (deg);Cluster 1 Count (Y1)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleY2_gt0 = new TH2F("h_cluSize_vs_AngleY2_gt0","Size 1 vs Angle Y2 Angle (deg);Cluster 1 Count (Y2)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleY3_gt0 = new TH2F("h_cluSize_vs_AngleY3_gt0","Size 1 vs Angle Y3 Angle (deg);Cluster 1 Count (Y3)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleY4_gt0 = new TH2F("h_cluSize_vs_AngleY4_gt0","Size 1 vs Angle Y4 Angle (deg);Cluster 1 Count (Y4)",100,60,120,20,0,20);

// per-layer size-vs-angle histograms (Y1..Y4)
TH2F *h_cluSize_vs_AngleY1_gt1 = new TH2F("h_cluSize_vs_AngleY1_gt1","Size 2 vs Angle Y1 Angle (deg);Cluster 2 Count (Y1)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleY2_gt1 = new TH2F("h_cluSize_vs_AngleY2_gt1","Size 2 vs Angle Y2 Angle (deg);Cluster 2 Count (Y2)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleY3_gt1 = new TH2F("h_cluSize_vs_AngleY3_gt1","Size 2 vs Angle Y3 Angle (deg);Cluster 2 Count (Y3)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleY4_gt1 = new TH2F("h_cluSize_vs_AngleY4_gt1","Size 2 vs Angle Y4 Angle (deg);Cluster 2 Count (Y4)",100,60,120,20,0,20);

// per-layer size-vs-angle histograms (Y1..Y4)
TH2F *h_cluSize_vs_AngleY1_gt2 = new TH2F("h_cluSize_vs_AngleY1_gt2","Size 3 vs Angle Y1 Angle (deg);Cluster 3 Count (Y1)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleY2_gt2 = new TH2F("h_cluSize_vs_AngleY2_gt2","Size 3 vs Angle Y2 Angle (deg);Cluster 3 Count (Y2)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleY3_gt2 = new TH2F("h_cluSize_vs_AngleY3_gt2","Size 3 vs Angle Y3 Angle (deg);Cluster 3 Count (Y3)",100,60,120,20,0,20);
TH2F *h_cluSize_vs_AngleY4_gt2 = new TH2F("h_cluSize_vs_AngleY4_gt2","Size 3 vs Angle Y4 Angle (deg);Cluster 3 Count (Y4)",100,60,120,20,0,20);

TH2F *h_cluCharge_vs_AngleY_gt0 = new TH2F("h_cluCharge_vs_AngleY_gt0","Charge vs Angle Y (nclus>0);Angle (deg);Total Charge",100,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleY_gt1 = new TH2F("h_cluCharge_vs_AngleY_gt1","Charge vs Angle Y (nclus>4);Angle (deg);Total Charge",100,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleY_gt2 = new TH2F("h_cluCharge_vs_AngleY_gt2","Charge vs Angle Y (nclus>8);Angle (deg);Total Charge",100,60,120,600,0,6000);

// --- Size vs Charge (separate per X layer and angle region)
TH2F *h_cluSize_vs_cluCharge_60X1 = new TH2F("h_cluSize_vs_cluCharge_60X1","Size vs cluCharge X1 (60-85 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_60X2 = new TH2F("h_cluSize_vs_cluCharge_60X2","Size vs cluCharge X2 (60-85 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_60X3 = new TH2F("h_cluSize_vs_cluCharge_60X3","Size vs cluCharge X3 (60-85 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_60X4 = new TH2F("h_cluSize_vs_cluCharge_60X4","Size vs cluCharge X4 (60-85 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);

TH2F *h_cluSize_vs_cluCharge_90X1 = new TH2F("h_cluSize_vs_cluCharge_90X1","Size vs cluCharge X1 (85-95 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_90X2 = new TH2F("h_cluSize_vs_cluCharge_90X2","Size vs cluCharge X2 (85-95 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_90X3 = new TH2F("h_cluSize_vs_cluCharge_90X3","Size vs cluCharge X3 (85-95 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_90X4 = new TH2F("h_cluSize_vs_cluCharge_90X4","Size vs cluCharge X4 (85-95 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);

TH2F *h_cluSize_vs_cluCharge_120X1 = new TH2F("h_cluSize_vs_cluCharge_120X1","Size vs cluCharge X1 (95-120 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_120X2 = new TH2F("h_cluSize_vs_cluCharge_120X2","Size vs cluCharge X2 (95-120 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_120X3 = new TH2F("h_cluSize_vs_cluCharge_120X3","Size vs cluCharge X3 (95-120 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_120X4 = new TH2F("h_cluSize_vs_cluCharge_120X4","Size vs cluCharge X4 (95-120 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);

// --- Size vs Charge (separate per Y layer and angle region)
TH2F *h_cluSize_vs_cluCharge_60Y1 = new TH2F("h_cluSize_vs_cluCharge_60Y1","Size vs cluCharge Y1 (60-85 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_60Y2 = new TH2F("h_cluSize_vs_cluCharge_60Y2","Size vs cluCharge Y2 (60-85 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_60Y3 = new TH2F("h_cluSize_vs_cluCharge_60Y3","Size vs cluCharge Y3 (60-85 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_60Y4 = new TH2F("h_cluSize_vs_cluCharge_60Y4","Size vs cluCharge Y4 (60-85 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);

TH2F *h_cluSize_vs_cluCharge_90Y1 = new TH2F("h_cluSize_vs_cluCharge_90Y1","Size vs cluCharge Y1 (85-95 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_90Y2 = new TH2F("h_cluSize_vs_cluCharge_90Y2","Size vs cluCharge Y2 (85-95 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_90Y3 = new TH2F("h_cluSize_vs_cluCharge_90Y3","Size vs cluCharge Y3 (85-95 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_90Y4 = new TH2F("h_cluSize_vs_cluCharge_90Y4","Size vs cluCharge Y4 (85-95 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);

TH2F *h_cluSize_vs_cluCharge_120Y1 = new TH2F("h_cluSize_vs_cluCharge_120Y1","Size vs cluCharge Y1 (95-120 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_120Y2 = new TH2F("h_cluSize_vs_cluCharge_120Y2","Size vs cluCharge Y2 (95-120 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_120Y3 = new TH2F("h_cluSize_vs_cluCharge_120Y3","Size vs cluCharge Y3 (95-120 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);
TH2F *h_cluSize_vs_cluCharge_120Y4 = new TH2F("h_cluSize_vs_cluCharge_120Y4","Size vs cluCharge Y4 (95-120 deg);Cluster size;Cluster charge",20,0,20,600,0,6000);

// --- Γραφήματα Leading Cluster Charge vs Strips X (1X έως 4X) ---
TH1F *h_leadClu_Charge_vs_Strips1X = new TH1F("h_leadClu_Charge_vs_Strips1X", "Leading Charge vs Strips 1X;Strip#;Charge", 640, 0, 640);
TH1F *h_leadClu_Charge_vs_Strips2X = new TH1F("h_leadClu_Charge_vs_Strips2X", "Leading Charge vs Strips 2X;Strip#;Charge", 640, 0, 640);
TH1F *h_leadClu_Charge_vs_Strips3X = new TH1F("h_leadClu_Charge_vs_Strips3X", "Leading Charge vs Strips 3X;Strip#;Charge", 640, 0, 640);
TH1F *h_leadClu_Charge_vs_Strips4X = new TH1F("h_leadClu_Charge_vs_Strips4X", "Leading Charge vs Strips 4X;Strip#;Charge", 640, 0, 640);

// --- Γραφήματα Leading Cluster Charge vs Strips Y (1Y έως 4Y) ---
TH1F *h_leadClu_Charge_vs_Strips1Y = new TH1F("h_leadClu_Charge_vs_Strips1Y", "Leading Charge vs Strips 1Y;Strip#;Charge", 640, 0, 640);
TH1F *h_leadClu_Charge_vs_Strips2Y = new TH1F("h_leadClu_Charge_vs_Strips2Y", "Leading Charge vs Strips 2Y;Strip#;Charge", 640, 0, 640);
TH1F *h_leadClu_Charge_vs_Strips3Y = new TH1F("h_leadClu_Charge_vs_Strips3Y", "Leading Charge vs Strips 3Y;Strip#;Charge", 640, 0, 640);
TH1F *h_leadClu_Charge_vs_Strips4Y = new TH1F("h_leadClu_Charge_vs_Strips4Y", "Leading Charge vs Strips 4Y;Strip#;Charge", 640, 0, 640);

// --- Per-layer Charge vs Angle (nclus>0, >1, >2) ---
TH2F *h_cluCharge_vs_AngleX1_gt0 = new TH2F("h_cluCharge_vs_AngleX1_gt0","Charge vs Angle X1 (nclus>0);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleX2_gt0 = new TH2F("h_cluCharge_vs_AngleX2_gt0","Charge vs Angle X2 (nclus>0);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleX3_gt0 = new TH2F("h_cluCharge_vs_AngleX3_gt0","Charge vs Angle X3 (nclus>0);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleX4_gt0 = new TH2F("h_cluCharge_vs_AngleX4_gt0","Charge vs Angle X4 (nclus>0);Angle (deg);Charge",60,60,120,600,0,6000);

TH2F *h_cluCharge_vs_AngleX1_gt1 = new TH2F("h_cluCharge_vs_AngleX1_gt1","Charge vs Angle X1 (nclus>1);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleX2_gt1 = new TH2F("h_cluCharge_vs_AngleX2_gt1","Charge vs Angle X2 (nclus>1);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleX3_gt1 = new TH2F("h_cluCharge_vs_AngleX3_gt1","Charge vs Angle X3 (nclus>1);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleX4_gt1 = new TH2F("h_cluCharge_vs_AngleX4_gt1","Charge vs Angle X4 (nclus>1);Angle (deg);Charge",60,60,120,600,0,6000);

TH2F *h_cluCharge_vs_AngleX1_gt2 = new TH2F("h_cluCharge_vs_AngleX1_gt2","Charge vs Angle X1 (nclus>2);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleX2_gt2 = new TH2F("h_cluCharge_vs_AngleX2_gt2","Charge vs Angle X2 (nclus>2);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleX3_gt2 = new TH2F("h_cluCharge_vs_AngleX3_gt2","Charge vs Angle X3 (nclus>2);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleX4_gt2 = new TH2F("h_cluCharge_vs_AngleX4_gt2","Charge vs Angle X4 (nclus>2);Angle (deg);Charge",60,60,120,600,0,6000);

TH2F *h_cluCharge_vs_AngleY1_gt0 = new TH2F("h_cluCharge_vs_AngleY1_gt0","Charge vs Angle Y1 (nclus>0);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleY2_gt0 = new TH2F("h_cluCharge_vs_AngleY2_gt0","Charge vs Angle Y2 (nclus>0);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleY3_gt0 = new TH2F("h_cluCharge_vs_AngleY3_gt0","Charge vs Angle Y3 (nclus>0);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleY4_gt0 = new TH2F("h_cluCharge_vs_AngleY4_gt0","Charge vs Angle Y4 (nclus>0);Angle (deg);Charge",60,60,120,600,0,6000);

TH2F *h_cluCharge_vs_AngleY1_gt1 = new TH2F("h_cluCharge_vs_AngleY1_gt1","Charge vs Angle Y1 (nclus>1);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleY2_gt1 = new TH2F("h_cluCharge_vs_AngleY2_gt1","Charge vs Angle Y2 (nclus>1);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleY3_gt1 = new TH2F("h_cluCharge_vs_AngleY3_gt1","Charge vs Angle Y3 (nclus>1);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleY4_gt1 = new TH2F("h_cluCharge_vs_AngleY4_gt1","Charge vs Angle Y4 (nclus>1);Angle (deg);Charge",60,60,120,600,0,6000);

TH2F *h_cluCharge_vs_AngleY1_gt2 = new TH2F("h_cluCharge_vs_AngleY1_gt2","Charge vs Angle Y1 (nclus>2);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleY2_gt2 = new TH2F("h_cluCharge_vs_AngleY2_gt2","Charge vs Angle Y2 (nclus>2);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleY3_gt2 = new TH2F("h_cluCharge_vs_AngleY3_gt2","Charge vs Angle Y3 (nclus>2);Angle (deg);Charge",60,60,120,600,0,6000);
TH2F *h_cluCharge_vs_AngleY4_gt2 = new TH2F("h_cluCharge_vs_AngleY4_gt2","Charge vs Angle Y4 (nclus>2);Angle (deg);Charge",60,60,120,600,0,6000);
//--- End of Charge vs Angle per layer ---

// --- Slope Analysis Histograms ---
TH1F *h_slope_X = new TH1F("h_slope_X", "Track Slope (X direction);Slope (a);Counts", 151, -75, 75);
TH1F *h_slope_Y = new TH1F("h_slope_Y", "Track Slope (Y direction);Slope (a);Counts", 151, -75, 75);
TH1F *h_intercept_X = new TH1F("h_intercept_X", "Track Intercept (X direction);Intercept (b) [mm];Counts", 9500, -4500, 5500);
TH1F *h_intercept_Y = new TH1F("h_intercept_Y", "Track Intercept (Y direction);Intercept (b) [mm];Counts", 9500, -4500, 5500);
TH2F *h_slope_vs_intercept_X = new TH2F("h_slope_vs_intercept_X", "Slope vs Intercept (X);Slope (a);Intercept (b)", 151, -75, 75, 9500, -4500, 5500);
TH2F *h_slope_vs_intercept_Y = new TH2F("h_slope_vs_intercept_Y", "Slope vs Intercept (Y);Slope (a);Intercept (b)", 151, -75, 75, 9500, -4500, 5500);
TH2F *h_angleXY_all = new TH2F("h_angleXY_all", "AngleY vs AngleX (all layers);angleX (deg);angleY (deg)", 180, 40, 140, 180, 40, 140);
TH1F *h_angleX_all = new TH1F("h_angleX_all", "AngleX (all layers);angleX (deg);events", 180, 40, 140);
TH1F *h_angleY_all = new TH1F("h_angleY_all", "AngleY (all layers);angleY (deg);events", 180, 40, 140);
TH1F *h_best3LayerAngleX_R2 = new TH1F("h_best3LayerAngleX_R2", "Best 3-layer angle by R^{2} (X);angleX (deg);events", 180, 40, 140);
TH1F *h_best3LayerAngleY_R2 = new TH1F("h_best3LayerAngleY_R2", "Best 3-layer angle by R^{2} (Y);angleY (deg);events", 180, 40, 140);
TH1F *h_best3LayerR2 = new TH1F("h_best3LayerR2", "Best 3-layer fit quality;R^{2};entries", 100, 0.0, 1.0);
TH1F *h_4LayerAngleX_R2 = new TH1F("h_4LayerAngleX_R2", "4-layer angle by R^{2} (X);angleX (deg);events", 180, 40, 140);
TH1F *h_4LayerAngleY_R2 = new TH1F("h_4LayerAngleY_R2", "4-layer angle by R^{2} (Y);angleY (deg);events", 180, 40, 140);
TH1F *h_4LayerR2 = new TH1F("h_4LayerR2", "4-layer fit quality;R^{2};entries", 100, 0.0, 1.0);

// New angle histograms with reduced chi2 selection (independent from existing histograms).
TH1F *h_angleX_new_chi2 = new TH1F("h_angleX_new_chi2", "AngleX (new, reduced #chi^{2} cut);angleX (deg);events", 180, 40, 140);
TH1F *h_angleY_new_chi2 = new TH1F("h_angleY_new_chi2", "AngleY (new, reduced #chi^{2} cut);angleY (deg);events", 180, 40, 140);
TH2F *h_angleXY_new_chi2 = new TH2F("h_angleXY_new_chi2", "AngleY vs AngleX (new, reduced #chi^{2} cut);angleX (deg);angleY (deg)", 180, 40, 140, 180, 40, 140);
TH1F *h_redChi2_new_X = new TH1F("h_redChi2_new_X", "Reduced #chi^{2} X (new selection);#chi^{2}/NDF;events", 150, 0.0, 15.0);
TH1F *h_redChi2_new_Y = new TH1F("h_redChi2_new_Y", "Reduced #chi^{2} Y (new selection);#chi^{2}/NDF;events", 150, 0.0, 15.0);

