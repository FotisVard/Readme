#define mmanalys_cxx

#include "mmanalys.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <map>
#include "TMath.h"
#include "TStyle.h"
#include <stdio.h>
#include <string>
#include "TLegend.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "histogram_booking.h"
#include "plot.h"
#include "MMHit.h"
#include "MMCluster.h"
#include "MMPacmanAlgo.h"


using namespace std;
using std::ofstream;

long evntnum=0;

static double noise_matrix[8][640] ={{0.0}};
static bool noise_loaded = false;

// Empirical resolution (sigma) per cluster size, indexed by size (0 unused)
static double empSigmaX[11] = {0.0};   // size 1..10
static double empSigmaY[11] = {0.0};
static bool useEmpiricalSigma = false;


// =========================================================================
// FUNCTION: Extract minimum non-zero charge threshold for every strip
// =========================================================================
void extractLowestChargeMatrix(TH2F* histList[8], double minChargeMatrix[8][640]) {
    
    // Names of the 8 detector planes
    const char* layerNames[8] = {
        "1X", "2X", "3X", "4X", 
        "1Y", "2Y", "3Y", "4Y"
    };

    for (int layer = 0; layer < 8; layer++) {
        TH2F* h = histList[layer];
        if (!h) continue;

        int nBinsX = h->GetNbinsX(); // 640 strips
        int nBinsY = h->GetNbinsY(); // Charge bins

        for (int stripBin = 1; stripBin <= nBinsX && stripBin <= 640; stripBin++) {
            int stripNumber = (int)h->GetXaxis()->GetBinCenter(stripBin);
            double lowestCharge = -1.0; // Default if strip has no hits

            // Scan from bottom to top along the Y-axis (Charge)
            for (int chargeBin = 1; chargeBin <= nBinsY; chargeBin++) {
                double binContent = h->GetBinContent(stripBin, chargeBin);
                
                // First bin with at least 1 count is the minimum recorded charge
                if (binContent > 0) {
                    lowestCharge = h->GetYaxis()->GetBinLowEdge(chargeBin);
                    break;
                }
            }

            // Store in the matrix: row = layer index (0-7), col = strip index (0-639)
            minChargeMatrix[layer][stripBin - 1] = lowestCharge;
        }
        
        std::cout << "Extracted minimum charges for Layer " << layerNames[layer] << std::endl;
    }
}

// =========================================================================
// FUNCTION: Extract empirical sigma per cluster size from residual histograms
// =========================================================================
void UpdateEmpiricalSigmas() {
    // Get the sigma histograms created by FitSlicesY
    TH1D* sigmaX = (TH1D*)gDirectory->Get("h_res_vs_sizeX_2");
    TH1D* sigmaY = (TH1D*)gDirectory->Get("h_res_vs_sizeY_2");
    if (!sigmaX || !sigmaY) {
        std::cerr << "Error: could not retrieve sigma histograms.\n";
        return;
    }
    for (int i = 1; i <= 10; ++i) {
        empSigmaX[i] = sigmaX->GetBinContent(i);
        empSigmaY[i] = sigmaY->GetBinContent(i);
    }
    // Print them so you can see the values
    std::cout << "\n=== EMPIRICAL SIGMAS ===\n";
    for (int i = 1; i <= 10; ++i)
        std::cout << "Size " << i << " : X = " << empSigmaX[i] 
                  << " mm, Y = " << empSigmaY[i] << " mm\n";
    std::cout << "=========================\n\n";
}

// --- GLOBAL MAPPING VARIABLES ---
static int stripMap[800];
static bool mapInitialized = false;

void initializeMap() {
    if (mapInitialized) return; 
    
    // Default identity mapping
    for(int i = 0; i < 800; i++) {
        stripMap[i] = i;
    }

    // =======================================================
    // --- EXPLICIT NEWMAP2 CALIBRATION (768 STRIPS) ---
    // =======================================================
    
    // APV 1
    stripMap[32] = 96; stripMap[33] = 97; stripMap[31] = 95; stripMap[34] = 98;
    stripMap[30] = 94; stripMap[35] = 99; stripMap[29] = 93; stripMap[36] = 100;
    stripMap[28] = 92; stripMap[37] = 101; stripMap[27] = 91; stripMap[38] = 102;
    stripMap[26] = 90; stripMap[39] = 103; stripMap[25] = 89; stripMap[40] = 104;
    stripMap[24] = 88; stripMap[41] = 105; stripMap[23] = 87; stripMap[42] = 106;
    stripMap[22] = 86; stripMap[43] = 107; stripMap[21] = 85; stripMap[44] = 108;
    stripMap[20] = 84; stripMap[45] = 109; stripMap[19] = 83; stripMap[46] = 110;
    stripMap[18] = 82; stripMap[47] = 111; stripMap[17] = 81; stripMap[48] = 112;
    stripMap[16] = 80; stripMap[49] = 113; stripMap[15] = 79; stripMap[50] = 114;
    stripMap[14] = 78; stripMap[51] = 115; stripMap[13] = 77; stripMap[52] = 116;
    stripMap[12] = 76; stripMap[53] = 117; stripMap[11] = 75; stripMap[54] = 118;
    stripMap[10] = 74; stripMap[55] = 119; stripMap[9] = 73; stripMap[56] = 120;
    stripMap[8] = 72; stripMap[57] = 121; stripMap[7] = 71; stripMap[58] = 122;
    stripMap[6] = 70; stripMap[59] = 123; stripMap[5] = 69; stripMap[60] = 124;
    stripMap[4] = 68; stripMap[61] = 125; stripMap[3] = 67; stripMap[62] = 126;
    stripMap[2] = 66; stripMap[63] = 127; stripMap[1] = 65; stripMap[64] = 128;
    stripMap[128] = 64; stripMap[65] = 1; stripMap[127] = 63; stripMap[66] = 2;
    stripMap[126] = 62; stripMap[67] = 3; stripMap[125] = 61; stripMap[68] = 4;
    stripMap[124] = 60; stripMap[69] = 5; stripMap[123] = 59; stripMap[70] = 6;
    stripMap[122] = 58; stripMap[71] = 7; stripMap[121] = 57; stripMap[72] = 8;
    stripMap[120] = 56; stripMap[73] = 9; stripMap[119] = 55; stripMap[74] = 10;
    stripMap[118] = 54; stripMap[75] = 11; stripMap[117] = 53; stripMap[76] = 12;
    stripMap[116] = 52; stripMap[77] = 13; stripMap[115] = 51; stripMap[78] = 14;
    stripMap[114] = 50; stripMap[79] = 15; stripMap[113] = 49; stripMap[80] = 16;
    stripMap[112] = 48; stripMap[81] = 17; stripMap[111] = 47; stripMap[82] = 18;
    stripMap[110] = 46; stripMap[83] = 19; stripMap[109] = 45; stripMap[84] = 20;
    stripMap[108] = 44; stripMap[85] = 21; stripMap[107] = 43; stripMap[86] = 22;
    stripMap[106] = 42; stripMap[87] = 23; stripMap[105] = 41; stripMap[88] = 24;
    stripMap[104] = 40; stripMap[89] = 25; stripMap[103] = 39; stripMap[90] = 26;
    stripMap[102] = 38; stripMap[91] = 27; stripMap[101] = 37; stripMap[92] = 28;
    stripMap[100] = 36; stripMap[93] = 29; stripMap[99] = 35; stripMap[94] = 30;
    stripMap[98] = 34; stripMap[95] = 31; stripMap[97] = 33; stripMap[96] = 32;

    // APV 2
    stripMap[160] = 224; stripMap[161] = 225; stripMap[159] = 223; stripMap[162] = 226;
    stripMap[158] = 222; stripMap[163] = 227; stripMap[157] = 221; stripMap[164] = 228;
    stripMap[156] = 220; stripMap[165] = 229; stripMap[155] = 219; stripMap[166] = 230;
    stripMap[154] = 218; stripMap[167] = 231; stripMap[153] = 217; stripMap[168] = 232;
    stripMap[152] = 216; stripMap[169] = 233; stripMap[151] = 215; stripMap[170] = 234;
    stripMap[150] = 214; stripMap[171] = 235; stripMap[149] = 213; stripMap[172] = 236;
    stripMap[148] = 212; stripMap[173] = 237; stripMap[147] = 211; stripMap[174] = 238;
    stripMap[146] = 210; stripMap[175] = 239; stripMap[145] = 209; stripMap[176] = 240;
    stripMap[144] = 208; stripMap[177] = 241; stripMap[143] = 207; stripMap[178] = 242;
    stripMap[142] = 206; stripMap[179] = 243; stripMap[141] = 205; stripMap[180] = 244;
    stripMap[140] = 204; stripMap[181] = 245; stripMap[139] = 203; stripMap[182] = 246;
    stripMap[138] = 202; stripMap[183] = 247; stripMap[137] = 201; stripMap[184] = 248;
    stripMap[136] = 200; stripMap[185] = 249; stripMap[135] = 199; stripMap[186] = 250;
    stripMap[134] = 198; stripMap[187] = 251; stripMap[133] = 197; stripMap[188] = 252;
    stripMap[132] = 196; stripMap[189] = 253; stripMap[131] = 195; stripMap[190] = 254;
    stripMap[130] = 194; stripMap[191] = 255; stripMap[129] = 193; stripMap[192] = 256;
    stripMap[256] = 192; stripMap[193] = 129; stripMap[255] = 191; stripMap[194] = 130;
    stripMap[254] = 190; stripMap[195] = 131; stripMap[253] = 189; stripMap[196] = 132;
    stripMap[252] = 188; stripMap[197] = 133; stripMap[251] = 187; stripMap[198] = 134;
    stripMap[250] = 186; stripMap[199] = 135; stripMap[249] = 185; stripMap[200] = 136;
    stripMap[248] = 184; stripMap[201] = 137; stripMap[247] = 183; stripMap[202] = 138;
    stripMap[246] = 182; stripMap[203] = 139; stripMap[245] = 181; stripMap[204] = 140;
    stripMap[244] = 180; stripMap[205] = 141; stripMap[243] = 179; stripMap[206] = 142;
    stripMap[242] = 178; stripMap[207] = 143; stripMap[241] = 177; stripMap[208] = 144;
    stripMap[240] = 176; stripMap[209] = 145; stripMap[239] = 175; stripMap[210] = 146;
    stripMap[238] = 174; stripMap[211] = 147; stripMap[237] = 173; stripMap[212] = 148;
    stripMap[236] = 172; stripMap[213] = 149; stripMap[235] = 171; stripMap[214] = 150;
    stripMap[234] = 170; stripMap[215] = 151; stripMap[233] = 169; stripMap[216] = 152;
    stripMap[232] = 168; stripMap[217] = 153; stripMap[231] = 167; stripMap[218] = 154;
    stripMap[230] = 166; stripMap[219] = 155; stripMap[229] = 165; stripMap[220] = 156;
    stripMap[228] = 164; stripMap[221] = 157; stripMap[227] = 163; stripMap[222] = 158;
    stripMap[226] = 162; stripMap[223] = 159; stripMap[225] = 161; stripMap[224] = 160;

    // APV 3
    stripMap[288] = 352; stripMap[289] = 353; stripMap[287] = 351; stripMap[290] = 354;
    stripMap[286] = 350; stripMap[291] = 355; stripMap[285] = 349; stripMap[292] = 356;
    stripMap[284] = 348; stripMap[293] = 357; stripMap[283] = 347; stripMap[294] = 358;
    stripMap[282] = 346; stripMap[295] = 359; stripMap[281] = 345; stripMap[296] = 360;
    stripMap[280] = 344; stripMap[297] = 361; stripMap[279] = 343; stripMap[298] = 362;
    stripMap[278] = 342; stripMap[299] = 363; stripMap[277] = 341; stripMap[300] = 364;
    stripMap[276] = 340; stripMap[301] = 365; stripMap[275] = 339; stripMap[302] = 366;
    stripMap[274] = 338; stripMap[303] = 367; stripMap[273] = 337; stripMap[304] = 368;
    stripMap[272] = 336; stripMap[305] = 369; stripMap[271] = 335; stripMap[306] = 370;
    stripMap[270] = 334; stripMap[307] = 371; stripMap[269] = 333; stripMap[308] = 372;
    stripMap[268] = 332; stripMap[309] = 373; stripMap[267] = 331; stripMap[310] = 374;
    stripMap[266] = 330; stripMap[311] = 375; stripMap[265] = 329; stripMap[312] = 376;
    stripMap[264] = 328; stripMap[313] = 377; stripMap[263] = 327; stripMap[314] = 378;
    stripMap[262] = 326; stripMap[315] = 379; stripMap[261] = 325; stripMap[316] = 380;
    stripMap[260] = 324; stripMap[317] = 381; stripMap[259] = 323; stripMap[318] = 382;
    stripMap[258] = 322; stripMap[319] = 383; stripMap[257] = 321; stripMap[320] = 384;
    stripMap[384] = 320; stripMap[321] = 257; stripMap[383] = 319; stripMap[322] = 258;
    stripMap[382] = 318; stripMap[323] = 259; stripMap[381] = 317; stripMap[324] = 260;
    stripMap[380] = 316; stripMap[325] = 261; stripMap[379] = 315; stripMap[326] = 262;
    stripMap[378] = 314; stripMap[327] = 263; stripMap[377] = 313; stripMap[328] = 264;
    stripMap[376] = 312; stripMap[329] = 265; stripMap[375] = 311; stripMap[330] = 266;
    stripMap[374] = 310; stripMap[331] = 267; stripMap[373] = 309; stripMap[332] = 268;
    stripMap[372] = 308; stripMap[333] = 269; stripMap[371] = 307; stripMap[334] = 270;
    stripMap[370] = 306; stripMap[335] = 271; stripMap[369] = 305; stripMap[336] = 272;
    stripMap[368] = 304; stripMap[337] = 273; stripMap[367] = 303; stripMap[338] = 274;
    stripMap[366] = 302; stripMap[339] = 275; stripMap[365] = 301; stripMap[340] = 276;
    stripMap[364] = 300; stripMap[341] = 277; stripMap[363] = 299; stripMap[342] = 278;
    stripMap[362] = 298; stripMap[343] = 279; stripMap[361] = 297; stripMap[344] = 280;
    stripMap[360] = 296; stripMap[345] = 281; stripMap[359] = 295; stripMap[346] = 282;
    stripMap[358] = 294; stripMap[347] = 283; stripMap[357] = 293; stripMap[348] = 284;
    stripMap[356] = 292; stripMap[349] = 285; stripMap[355] = 291; stripMap[350] = 286;
    stripMap[354] = 290; stripMap[351] = 287; stripMap[353] = 289; stripMap[352] = 288;

    // APV 4
    stripMap[416] = 480; stripMap[417] = 481; stripMap[415] = 479; stripMap[418] = 482;
    stripMap[414] = 478; stripMap[419] = 483; stripMap[413] = 477; stripMap[420] = 484;
    stripMap[412] = 476; stripMap[421] = 485; stripMap[411] = 475; stripMap[422] = 486;
    stripMap[410] = 474; stripMap[423] = 487; stripMap[409] = 473; stripMap[424] = 488;
    stripMap[408] = 472; stripMap[425] = 489; stripMap[407] = 471; stripMap[426] = 490;
    stripMap[406] = 470; stripMap[427] = 491; stripMap[405] = 469; stripMap[428] = 492;
    stripMap[404] = 468; stripMap[429] = 493; stripMap[403] = 467; stripMap[430] = 494;
    stripMap[402] = 466; stripMap[431] = 495; stripMap[401] = 465; stripMap[432] = 496;
    stripMap[400] = 464; stripMap[433] = 497; stripMap[399] = 463; stripMap[434] = 498;
    stripMap[398] = 462; stripMap[435] = 499; stripMap[397] = 461; stripMap[436] = 500;
    stripMap[396] = 460; stripMap[437] = 501; stripMap[395] = 459; stripMap[438] = 502;
    stripMap[394] = 458; stripMap[439] = 503; stripMap[393] = 457; stripMap[440] = 504;
    stripMap[392] = 456; stripMap[441] = 505; stripMap[391] = 455; stripMap[442] = 506;
    stripMap[390] = 454; stripMap[443] = 507; stripMap[389] = 453; stripMap[444] = 508;
    stripMap[388] = 452; stripMap[445] = 509; stripMap[387] = 451; stripMap[446] = 510;
    stripMap[386] = 450; stripMap[447] = 511; stripMap[385] = 449; stripMap[448] = 512;
    stripMap[512] = 448; stripMap[449] = 385; stripMap[511] = 447; stripMap[450] = 386;
    stripMap[510] = 446; stripMap[451] = 387; stripMap[509] = 445; stripMap[452] = 388;
    stripMap[508] = 444; stripMap[453] = 389; stripMap[507] = 443; stripMap[454] = 390;
    stripMap[506] = 442; stripMap[455] = 391; stripMap[505] = 441; stripMap[456] = 392;
    stripMap[504] = 440; stripMap[457] = 393; stripMap[503] = 439; stripMap[458] = 394;
    stripMap[502] = 438; stripMap[459] = 395; stripMap[501] = 437; stripMap[460] = 396;
    stripMap[500] = 436; stripMap[461] = 397; stripMap[499] = 435; stripMap[462] = 398;
    stripMap[498] = 434; stripMap[463] = 399; stripMap[497] = 433; stripMap[464] = 400;
    stripMap[496] = 432; stripMap[465] = 401; stripMap[495] = 431; stripMap[466] = 402;
    stripMap[494] = 430; stripMap[467] = 403; stripMap[493] = 429; stripMap[468] = 404;
    stripMap[492] = 428; stripMap[469] = 405; stripMap[491] = 427; stripMap[470] = 406;
    stripMap[490] = 426; stripMap[471] = 407; stripMap[489] = 425; stripMap[472] = 408;
    stripMap[488] = 424; stripMap[473] = 409; stripMap[487] = 423; stripMap[474] = 410;
    stripMap[486] = 422; stripMap[475] = 411; stripMap[485] = 421; stripMap[476] = 412;
    stripMap[484] = 420; stripMap[477] = 413; stripMap[483] = 419; stripMap[478] = 414;
    stripMap[482] = 418; stripMap[479] = 415; stripMap[481] = 417; stripMap[480] = 416;

    // APV 5
    stripMap[544] = 608; stripMap[545] = 609; stripMap[543] = 607; stripMap[546] = 610;
    stripMap[542] = 606; stripMap[547] = 611; stripMap[541] = 605; stripMap[548] = 612;
    stripMap[540] = 604; stripMap[549] = 613; stripMap[539] = 603; stripMap[550] = 614;
    stripMap[538] = 602; stripMap[551] = 615; stripMap[537] = 601; stripMap[552] = 616;
    stripMap[536] = 600; stripMap[553] = 617; stripMap[535] = 599; stripMap[554] = 618;
    stripMap[534] = 598; stripMap[555] = 619; stripMap[533] = 597; stripMap[556] = 620;
    stripMap[532] = 596; stripMap[557] = 621; stripMap[531] = 595; stripMap[558] = 622;
    stripMap[530] = 594; stripMap[559] = 623; stripMap[529] = 593; stripMap[560] = 624;
    stripMap[528] = 592; stripMap[561] = 625; stripMap[527] = 591; stripMap[562] = 626;
    stripMap[526] = 590; stripMap[563] = 627; stripMap[525] = 589; stripMap[564] = 628;
    stripMap[524] = 588; stripMap[565] = 629; stripMap[523] = 587; stripMap[566] = 630;
    stripMap[522] = 586; stripMap[567] = 631; stripMap[521] = 585; stripMap[568] = 632;
    stripMap[520] = 584; stripMap[569] = 633; stripMap[519] = 583; stripMap[570] = 634;
    stripMap[518] = 582; stripMap[571] = 635; stripMap[517] = 581; stripMap[572] = 636;
    stripMap[516] = 580; stripMap[573] = 637; stripMap[515] = 579; stripMap[574] = 638;
    stripMap[514] = 578; stripMap[575] = 639; stripMap[513] = 577; stripMap[576] = 640;
    stripMap[640] = 576; stripMap[577] = 513; stripMap[639] = 575; stripMap[578] = 514;
    stripMap[638] = 574; stripMap[579] = 515; stripMap[637] = 573; stripMap[580] = 516;
    stripMap[636] = 572; stripMap[581] = 517; stripMap[635] = 571; stripMap[582] = 518;
    stripMap[634] = 570; stripMap[583] = 519; stripMap[633] = 569; stripMap[584] = 520;
    stripMap[632] = 568; stripMap[585] = 521; stripMap[631] = 567; stripMap[586] = 522;
    stripMap[630] = 566; stripMap[587] = 523; stripMap[629] = 565; stripMap[588] = 524;
    stripMap[628] = 564; stripMap[589] = 525; stripMap[627] = 563; stripMap[590] = 526;
    stripMap[626] = 562; stripMap[591] = 527; stripMap[625] = 561; stripMap[592] = 528;
    stripMap[624] = 560; stripMap[593] = 529; stripMap[623] = 559; stripMap[594] = 530;
    stripMap[622] = 558; stripMap[595] = 531; stripMap[621] = 557; stripMap[596] = 532;
    stripMap[620] = 556; stripMap[597] = 533; stripMap[619] = 555; stripMap[598] = 534;
    stripMap[618] = 554; stripMap[599] = 535; stripMap[617] = 553; stripMap[600] = 536;
    stripMap[616] = 552; stripMap[601] = 537; stripMap[615] = 551; stripMap[602] = 538;
    stripMap[614] = 550; stripMap[603] = 539; stripMap[613] = 549; stripMap[604] = 540;
    stripMap[612] = 548; stripMap[605] = 541; stripMap[611] = 547; stripMap[606] = 542;
    stripMap[610] = 546; stripMap[607] = 543; stripMap[609] = 545; stripMap[608] = 544;

    // APV 6 (Mirrored backwards)
    stripMap[673] = 736; stripMap[672] = 737; stripMap[674] = 735; stripMap[671] = 738;
    stripMap[675] = 734; stripMap[670] = 739; stripMap[676] = 733; stripMap[669] = 740;
    stripMap[677] = 732; stripMap[668] = 741; stripMap[678] = 731; stripMap[667] = 742;
    stripMap[679] = 730; stripMap[666] = 743; stripMap[680] = 729; stripMap[665] = 744;
    stripMap[681] = 728; stripMap[664] = 745; stripMap[682] = 727; stripMap[663] = 746;
    stripMap[683] = 726; stripMap[662] = 747; stripMap[684] = 725; stripMap[661] = 748;
    stripMap[685] = 724; stripMap[660] = 749; stripMap[686] = 723; stripMap[659] = 750;
    stripMap[687] = 722; stripMap[658] = 751; stripMap[688] = 721; stripMap[657] = 752;
    stripMap[689] = 720; stripMap[656] = 753; stripMap[690] = 719; stripMap[655] = 754;
    stripMap[691] = 718; stripMap[654] = 755; stripMap[692] = 717; stripMap[653] = 756;
    stripMap[693] = 716; stripMap[652] = 757; stripMap[694] = 715; stripMap[651] = 758;
    stripMap[695] = 714; stripMap[650] = 759; stripMap[696] = 713; stripMap[649] = 760;
    stripMap[697] = 712; stripMap[648] = 761; stripMap[698] = 711; stripMap[647] = 762;
    stripMap[699] = 710; stripMap[646] = 763; stripMap[700] = 709; stripMap[645] = 764;
    stripMap[701] = 708; stripMap[644] = 765; stripMap[702] = 707; stripMap[643] = 766;
    stripMap[703] = 706; stripMap[642] = 767; stripMap[704] = 705; stripMap[641] = 768;
    stripMap[705] = 704; stripMap[768] = 641; stripMap[706] = 703; stripMap[767] = 642;
    stripMap[707] = 702; stripMap[766] = 643; stripMap[708] = 701; stripMap[765] = 644;
    stripMap[709] = 700; stripMap[764] = 645; stripMap[710] = 699; stripMap[763] = 646;
    stripMap[711] = 698; stripMap[762] = 647; stripMap[712] = 697; stripMap[761] = 648;
    stripMap[713] = 696; stripMap[760] = 649; stripMap[714] = 695; stripMap[759] = 650;
    stripMap[715] = 694; stripMap[758] = 651; stripMap[716] = 693; stripMap[757] = 652;
    stripMap[717] = 692; stripMap[756] = 653; stripMap[718] = 691; stripMap[755] = 654;
    stripMap[719] = 690; stripMap[754] = 655; stripMap[720] = 689; stripMap[753] = 656;
    stripMap[721] = 688; stripMap[752] = 657; stripMap[722] = 687; stripMap[751] = 658;
    stripMap[723] = 686; stripMap[750] = 659; stripMap[724] = 685; stripMap[749] = 660;
    stripMap[725] = 684; stripMap[748] = 661; stripMap[726] = 683; stripMap[747] = 662;
    stripMap[727] = 682; stripMap[746] = 663; stripMap[728] = 681; stripMap[745] = 664;
    stripMap[729] = 680; stripMap[744] = 665; stripMap[730] = 679; stripMap[743] = 666;
    stripMap[731] = 678; stripMap[742] = 667; stripMap[732] = 677; stripMap[741] = 668;
    stripMap[733] = 676; stripMap[740] = 669; stripMap[734] = 675; stripMap[739] = 670;
    stripMap[735] = 674; stripMap[738] = 671; stripMap[736] = 673; stripMap[737] = 672;

    mapInitialized = true;
}
/*
void initializeMap() {
    if (mapInitialized) return;

    // Default identity mapping
    for(int i = 0; i < 800; i++) stripMap[i] = i;

        stripMap[0]=575; stripMap[1]=574; stripMap[2]=573; stripMap[3]=572; stripMap[4]=571; stripMap[5]=570; stripMap[6]=569; stripMap[7]=568; stripMap[8]=567; stripMap[9]=566;
        stripMap[10]=565; stripMap[11]=564; stripMap[12]=563; stripMap[13]=562; stripMap[14]=561; stripMap[15]=560; stripMap[16]=559; stripMap[17]=558; stripMap[18]=557; stripMap[19]=556;
        stripMap[20]=555; stripMap[21]=554; stripMap[22]=553; stripMap[23]=552; stripMap[24]=551; stripMap[25]=550; stripMap[26]=549; stripMap[27]=548; stripMap[28]=547; stripMap[29]=546;
        stripMap[30]=545; stripMap[31]=544; stripMap[32]=543; stripMap[33]=542; stripMap[34]=541; stripMap[35]=540; stripMap[36]=539; stripMap[37]=538; stripMap[38]=537; stripMap[39]=536;
        stripMap[40]=535; stripMap[41]=534; stripMap[42]=533; stripMap[43]=532; stripMap[44]=531; stripMap[45]=530; stripMap[46]=529; stripMap[47]=528; stripMap[48]=527; stripMap[49]=526;
        stripMap[50]=525; stripMap[51]=524; stripMap[52]=523; stripMap[53]=522; stripMap[54]=521; stripMap[55]=520; stripMap[56]=519; stripMap[57]=518; stripMap[58]=517; stripMap[59]=516;
        stripMap[60]=515; stripMap[61]=514; stripMap[62]=513; stripMap[63]=512;
        
        stripMap[64]=639; stripMap[65]=638; stripMap[66]=637; stripMap[67]=636; stripMap[68]=635; stripMap[69]=634; stripMap[70]=633; stripMap[71]=632; stripMap[72]=631; stripMap[73]=630;
        stripMap[74]=629; stripMap[75]=628; stripMap[76]=627; stripMap[77]=626; stripMap[78]=625; stripMap[79]=624; stripMap[80]=623; stripMap[81]=622; stripMap[82]=621; stripMap[83]=620;
        stripMap[84]=619; stripMap[85]=618; stripMap[86]=617; stripMap[87]=616; stripMap[88]=615; stripMap[89]=614; stripMap[90]=613; stripMap[91]=612; stripMap[92]=611; stripMap[93]=610;
        stripMap[94]=609; stripMap[95]=608; stripMap[96]=607; stripMap[97]=606; stripMap[98]=605; stripMap[99]=604; stripMap[100]=603; stripMap[101]=602; stripMap[102]=601; stripMap[103]=600;
        stripMap[104]=599; stripMap[105]=598; stripMap[106]=597; stripMap[107]=596; stripMap[108]=595; stripMap[109]=594; stripMap[110]=593; stripMap[111]=592; stripMap[112]=591; stripMap[113]=590;
        stripMap[114]=589; stripMap[115]=588; stripMap[116]=587; stripMap[117]=586; stripMap[118]=585; stripMap[119]=584; stripMap[120]=583; stripMap[121]=582; stripMap[122]=581; stripMap[123]=580;
        stripMap[124]=579; stripMap[125]=578; stripMap[126]=577; stripMap[127]=576;

        stripMap[128]=384; stripMap[129]=385; stripMap[130]=386; stripMap[131]=387; stripMap[132]=388; stripMap[133]=389; stripMap[134]=390; stripMap[135]=391; stripMap[136]=392; stripMap[137]=393;
        stripMap[138]=394; stripMap[139]=395; stripMap[140]=396; stripMap[141]=397; stripMap[142]=398; stripMap[143]=399; stripMap[144]=400; stripMap[145]=401; stripMap[146]=402; stripMap[147]=403;
        stripMap[148]=404; stripMap[149]=405; stripMap[150]=406; stripMap[151]=407; stripMap[152]=408; stripMap[153]=409; stripMap[154]=410; stripMap[155]=411; stripMap[156]=412; stripMap[157]=413;
        stripMap[158]=414; stripMap[159]=415; stripMap[160]=416; stripMap[161]=417; stripMap[162]=418; stripMap[163]=419; stripMap[164]=420; stripMap[165]=421; stripMap[166]=422; stripMap[167]=423;
        stripMap[168]=424; stripMap[169]=425; stripMap[170]=426; stripMap[171]=427; stripMap[172]=428; stripMap[173]=429; stripMap[174]=430; stripMap[175]=431; stripMap[176]=432; stripMap[177]=433;
        stripMap[178]=434; stripMap[179]=435; stripMap[180]=436; stripMap[181]=437; stripMap[182]=438; stripMap[183]=439; stripMap[184]=440; stripMap[185]=441; stripMap[186]=442; stripMap[187]=443;
        stripMap[188]=444; stripMap[189]=445; stripMap[190]=446; stripMap[191]=447; 

        stripMap[192]=511; stripMap[193]=510; stripMap[194]=509; stripMap[195]=508; stripMap[196]=507; stripMap[197]=506; stripMap[198]=505; stripMap[199]=504; stripMap[200]=503; stripMap[201]=502;
        stripMap[202]=501; stripMap[203]=500; stripMap[204]=499; stripMap[205]=498; stripMap[206]=497; stripMap[207]=496; stripMap[208]=495; stripMap[209]=494; stripMap[210]=493; stripMap[211]=492;
        stripMap[212]=491; stripMap[213]=490; stripMap[214]=489; stripMap[215]=488; stripMap[216]=487; stripMap[217]=486; stripMap[218]=485; stripMap[219]=484; stripMap[220]=483; stripMap[221]=482;
        stripMap[222]=481; stripMap[223]=480; stripMap[224]=479; stripMap[225]=478; stripMap[226]=477; stripMap[227]=476; stripMap[228]=475; stripMap[229]=474; stripMap[230]=473; stripMap[231]=472;
        stripMap[232]=471; stripMap[233]=470; stripMap[234]=469; stripMap[235]=468; stripMap[236]=467; stripMap[237]=466; stripMap[238]=465; stripMap[239]=464; stripMap[240]=463; stripMap[241]=462;
        stripMap[242]=461; stripMap[243]=460; stripMap[244]=459; stripMap[245]=458; stripMap[246]=457; stripMap[247]=456; stripMap[248]=455; stripMap[249]=454; stripMap[250]=453; stripMap[251]=452;
        stripMap[252]=451; stripMap[253]=450; stripMap[254]=449; stripMap[255]=448;

        stripMap[256]=319; stripMap[257]=318; stripMap[258]=317; stripMap[259]=316; stripMap[260]=315; stripMap[261]=314; stripMap[262]=313; stripMap[263]=312; stripMap[264]=311; stripMap[265]=310;
        stripMap[266]=309; stripMap[267]=308; stripMap[268]=307; stripMap[269]=306; stripMap[270]=305; stripMap[271]=304; stripMap[272]=303; stripMap[273]=302; stripMap[274]=301; stripMap[275]=300;
        stripMap[276]=299; stripMap[277]=298; stripMap[278]=297; stripMap[279]=296; stripMap[280]=295; stripMap[281]=294; stripMap[282]=293; stripMap[283]=292; stripMap[284]=291; stripMap[285]=290;
        stripMap[286]=289; stripMap[287]=288; stripMap[288]=287; stripMap[289]=286; stripMap[290]=285; stripMap[291]=284; stripMap[292]=283; stripMap[293]=282; stripMap[294]=281; stripMap[295]=280;
        stripMap[296]=279; stripMap[297]=278; stripMap[298]=277; stripMap[299]=276; stripMap[300]=275; stripMap[301]=274; stripMap[302]=273; stripMap[303]=272; stripMap[304]=271; stripMap[305]=270;
        stripMap[306]=269; stripMap[307]=268; stripMap[308]=267; stripMap[309]=266; stripMap[310]=265; stripMap[311]=264; stripMap[312]=263; stripMap[313]=262; stripMap[314]=261; stripMap[315]=260;
        stripMap[316]=259; stripMap[317]=258; stripMap[318]=257; stripMap[319]=256;
        
        stripMap[320]=383; stripMap[321]=382; stripMap[322]=381; stripMap[323]=380; stripMap[324]=379; stripMap[325]=378; stripMap[326]=377; stripMap[327]=376; stripMap[328]=375; stripMap[329]=374;
        stripMap[330]=373; stripMap[331]=372; stripMap[332]=371; stripMap[333]=370; stripMap[334]=369; stripMap[335]=368; stripMap[336]=367; stripMap[337]=366; stripMap[338]=365; stripMap[339]=364;
        stripMap[340]=363; stripMap[341]=362; stripMap[342]=361; stripMap[343]=360; stripMap[344]=359; stripMap[345]=358; stripMap[346]=357; stripMap[347]=356; stripMap[348]=355; stripMap[349]=354;
        stripMap[350]=353; stripMap[351]=352; stripMap[352]=351; stripMap[353]=350; stripMap[354]=349; stripMap[355]=348; stripMap[356]=347; stripMap[357]=346; stripMap[358]=345; stripMap[359]=344; 
        stripMap[360]=343; stripMap[361]=342; stripMap[362]=341; stripMap[363]=340; stripMap[364]=339; stripMap[365]=338; stripMap[366]=337; stripMap[367]=336; stripMap[368]=335; stripMap[369]=334; 
        stripMap[370]=333; stripMap[371]=332; stripMap[372]=331; stripMap[373]=330; stripMap[374]=329; stripMap[375]=328; stripMap[376]=327; stripMap[377]=326; stripMap[378]=325; stripMap[379]=324; 
        stripMap[380]=323; stripMap[381]=322; stripMap[382]=321; stripMap[383]=320;
        
        stripMap[384]=192; stripMap[385]=191; stripMap[386]=190; stripMap[387]=189; stripMap[388]=188; stripMap[389]=187; stripMap[390]=186; stripMap[391]=185; stripMap[392]=184; stripMap[393]=183;
        stripMap[394]=182; stripMap[395]=181; stripMap[396]=180; stripMap[397]=179; stripMap[398]=178; stripMap[399]=177; stripMap[400]=176; stripMap[401]=175; stripMap[402]=174; stripMap[403]=173;
        stripMap[404]=172; stripMap[405]=171; stripMap[406]=170; stripMap[407]=169; stripMap[408]=168; stripMap[409]=167; stripMap[410]=166; stripMap[411]=165; stripMap[412]=164; stripMap[413]=163;
        stripMap[414]=162; stripMap[415]=161; stripMap[416]=160; stripMap[417]=159; stripMap[418]=158; stripMap[419]=157; stripMap[420]=156; stripMap[421]=155; stripMap[422]=154; stripMap[423]=153;
        stripMap[424]=152; stripMap[425]=151; stripMap[426]=150; stripMap[427]=149; stripMap[428]=148; stripMap[429]=147; stripMap[430]=146; stripMap[431]=145; stripMap[432]=144; stripMap[433]=143;
        stripMap[434]=142; stripMap[435]=141; stripMap[436]=140; stripMap[437]=139; stripMap[438]=138; stripMap[439]=137; stripMap[440]=136; stripMap[441]=135; stripMap[442]=134; stripMap[443]=133;
        stripMap[444]=132; stripMap[445]=131; stripMap[446]=130; stripMap[447]=129; stripMap[448]=128;

        stripMap[449]=255; stripMap[450]=254; stripMap[451]=253; stripMap[452]=252; stripMap[453]=251; stripMap[454]=250; stripMap[455]=249; stripMap[456]=248; stripMap[457]=247; stripMap[458]=246;
        stripMap[459]=245; stripMap[460]=244; stripMap[461]=243; stripMap[462]=242; stripMap[463]=241; stripMap[464]=240; stripMap[465]=239; stripMap[466]=238; stripMap[467]=237; stripMap[468]=236;
        stripMap[469]=235; stripMap[470]=234; stripMap[471]=233; stripMap[472]=232; stripMap[473]=231; stripMap[474]=230; stripMap[475]=229; stripMap[476]=228; stripMap[477]=227; stripMap[478]=226;
        stripMap[479]=225; stripMap[480]=224; stripMap[481]=223; stripMap[482]=222; stripMap[483]=221; stripMap[484]=220; stripMap[485]=219; stripMap[486]=218; stripMap[487]=217; stripMap[488]=216;
        stripMap[489]=215; stripMap[490]=214; stripMap[491]=213; stripMap[492]=212; stripMap[493]=211; stripMap[494]=210; stripMap[495]=209; stripMap[496]=208; stripMap[497]=207; stripMap[498]=206;
        stripMap[499]=205; stripMap[500]=204; stripMap[501]=203; stripMap[502]=202; stripMap[503]=201; stripMap[504]=200; stripMap[505]=199; stripMap[506]=198; stripMap[507]=197; stripMap[508]=196;
        stripMap[509]=195; stripMap[510]=194; stripMap[511]=193;

        stripMap[512]=63; stripMap[513]=62; stripMap[514]=61; stripMap[515]=60; stripMap[516]=59; stripMap[517]=58; stripMap[518]=57; stripMap[519]=56; stripMap[520]=55; stripMap[521]=54;
        stripMap[522]=53; stripMap[523]=52; stripMap[524]=51; stripMap[525]=50; stripMap[526]=49; stripMap[527]=48; stripMap[528]=47; stripMap[529]=46; stripMap[530]=45; stripMap[531]=44;
        stripMap[532]=43; stripMap[533]=42; stripMap[534]=41; stripMap[535]=40; stripMap[536]=39; stripMap[537]=38; stripMap[538]=37; stripMap[539]=36; stripMap[540]=35; stripMap[541]=34;
        stripMap[542]=33; stripMap[543]=32; stripMap[544]=31; stripMap[545]=30; stripMap[546]=29; stripMap[547]=28; stripMap[548]=27; stripMap[549]=26; stripMap[550]=25; stripMap[551]=24;
        stripMap[552]=23; stripMap[553]=22; stripMap[554]=21; stripMap[555]=20; stripMap[556]=19; stripMap[557]=18; stripMap[558]=17; stripMap[559]=16; stripMap[560]=15; stripMap[561]=14;
        stripMap[562]=13; stripMap[563]=12; stripMap[564]=11; stripMap[565]=10; stripMap[566]=9;  stripMap[567]=8;  stripMap[568]=7;  stripMap[569]=6;  stripMap[570]=5;  stripMap[571]=4;
        stripMap[572]=3;  stripMap[573]=2;  stripMap[574]=1;  stripMap[575]=0;
        
        stripMap[576]=127; stripMap[577]=126; stripMap[578]=125; stripMap[579]=124; stripMap[580]=123; stripMap[581]=122; stripMap[582]=121; stripMap[583]=120; stripMap[584]=119; stripMap[585]=118;
        stripMap[586]=117; stripMap[587]=116; stripMap[588]=115; stripMap[589]=114; stripMap[590]=113; stripMap[591]=112; stripMap[592]=111; stripMap[593]=110; stripMap[594]=109; stripMap[595]=108;
        stripMap[596]=107; stripMap[597]=106; stripMap[598]=105; stripMap[599]=104; stripMap[600]=103; stripMap[601]=102; stripMap[602]=101; stripMap[603]=100; stripMap[604]=99;  stripMap[605]=98;
        stripMap[606]=97;  stripMap[607]=96;  stripMap[608]=95;  stripMap[609]=94;  stripMap[610]=93;  stripMap[611]=92;  stripMap[612]=91;  stripMap[613]=90;  stripMap[614]=89;  stripMap[615]=88;
        stripMap[616]=87;  stripMap[617]=86;  stripMap[618]=85;  stripMap[619]=84;  stripMap[620]=83;  stripMap[621]=82;  stripMap[622]=81;  stripMap[623]=80;  stripMap[624]=79;  stripMap[625]=78;
        stripMap[626]=77;  stripMap[627]=76;  stripMap[628]=75;  stripMap[629]=74;  stripMap[630]=73;  stripMap[631]=72;  stripMap[632]=71;  stripMap[633]=70;  stripMap[634]=69;  stripMap[635]=68;
        stripMap[636]=67;  stripMap[637]=66;  stripMap[638]=65;  stripMap[639]=64;
        mapInitialized = true;
    }
*/
/*
// --- INITIALIZE MAP FUNCTION ---
void initializeMap() {
    if (mapInitialized) return; 
    
    // Default identity mapping
    for(int i = 0; i < 800; i++) stripMap[i] = i;

    // Apply the specific swaps from your NEWMAP.txt calibration
    stripMap[32] = 33; stripMap[33] = 32; stripMap[31] = 34; stripMap[34] = 31;
    stripMap[30] = 35; stripMap[35] = 30; stripMap[29] = 36; stripMap[36] = 29;
    stripMap[28] = 37; stripMap[37] = 28; stripMap[27] = 38; stripMap[38] = 27;
    stripMap[26] = 39; stripMap[39] = 26; stripMap[25] = 40; stripMap[40] = 25;
    stripMap[24] = 41; stripMap[41] = 24; stripMap[23] = 42; stripMap[42] = 23;
    stripMap[22] = 43; stripMap[43] = 22; stripMap[21] = 44; stripMap[44] = 21;
    stripMap[20] = 45; stripMap[45] = 20; stripMap[19] = 46; stripMap[46] = 19;
    stripMap[18] = 47; stripMap[47] = 18; stripMap[17] = 48; stripMap[48] = 17;
    stripMap[16] = 49; stripMap[49] = 16; stripMap[15] = 50; stripMap[50] = 15;
    stripMap[14] = 51; stripMap[51] = 14; stripMap[13] = 52; stripMap[52] = 13;
    stripMap[12] = 53; stripMap[53] = 12; stripMap[11] = 54; stripMap[54] = 11;
    stripMap[10] = 55; stripMap[55] = 10; stripMap[9] = 56; stripMap[56] = 9;
    stripMap[8] = 57; stripMap[57] = 8; stripMap[7] = 58; stripMap[58] = 7;
    stripMap[6] = 59; stripMap[59] = 6; stripMap[5] = 60; stripMap[60] = 5;
    stripMap[4] = 61; stripMap[61] = 4; stripMap[3] = 62; stripMap[62] = 3;
    stripMap[2] = 63; stripMap[63] = 2; stripMap[1] = 64; stripMap[64] = 1;
    
    stripMap[128] = 65; stripMap[65] = 128; stripMap[127] = 66; stripMap[66] = 127;
    stripMap[126] = 67; stripMap[67] = 126; stripMap[125] = 68; stripMap[68] = 125;
    stripMap[124] = 69; stripMap[69] = 124; stripMap[123] = 70; stripMap[70] = 123;
    stripMap[122] = 71; stripMap[71] = 122; stripMap[121] = 72; stripMap[72] = 121;
    stripMap[120] = 73; stripMap[73] = 120; stripMap[119] = 74; stripMap[74] = 119;
    stripMap[118] = 75; stripMap[75] = 118; stripMap[117] = 76; stripMap[76] = 117;
    stripMap[116] = 77; stripMap[77] = 116; stripMap[115] = 78; stripMap[78] = 115;
    stripMap[114] = 79; stripMap[79] = 114; stripMap[113] = 80; stripMap[80] = 113;
    stripMap[112] = 81; stripMap[81] = 112; stripMap[111] = 82; stripMap[82] = 111;
    stripMap[110] = 83; stripMap[83] = 110; stripMap[109] = 84; stripMap[84] = 109;
    stripMap[108] = 85; stripMap[85] = 108; stripMap[107] = 86; stripMap[86] = 107;
    stripMap[106] = 87; stripMap[87] = 106; stripMap[105] = 88; stripMap[88] = 105;
    stripMap[104] = 89; stripMap[89] = 104; stripMap[103] = 90; stripMap[90] = 103;
    stripMap[102] = 91; stripMap[91] = 102; stripMap[101] = 92; stripMap[92] = 101;
    stripMap[100] = 93; stripMap[93] = 100; stripMap[99] = 94; stripMap[94] = 99;
    stripMap[98] = 95; stripMap[95] = 98; stripMap[97] = 96; stripMap[96] = 97;
    
    stripMap[160] = 161; stripMap[161] = 160; stripMap[159] = 162; stripMap[162] = 159;
    stripMap[158] = 163; stripMap[163] = 158; stripMap[157] = 164; stripMap[164] = 157;
    stripMap[156] = 165; stripMap[165] = 156; stripMap[155] = 166; stripMap[166] = 155;
    stripMap[154] = 167; stripMap[167] = 154; stripMap[153] = 168; stripMap[168] = 153;
    stripMap[152] = 169; stripMap[169] = 152; stripMap[151] = 170; stripMap[170] = 151;
    stripMap[150] = 171; stripMap[171] = 150; stripMap[149] = 172; stripMap[172] = 149;
    stripMap[148] = 173; stripMap[173] = 148; stripMap[147] = 174; stripMap[174] = 147;
    stripMap[146] = 175; stripMap[175] = 146; stripMap[145] = 176; stripMap[176] = 145;
    stripMap[144] = 177; stripMap[177] = 144; stripMap[143] = 178; stripMap[178] = 143;
    stripMap[142] = 179; stripMap[179] = 142; stripMap[141] = 180; stripMap[180] = 141;
    stripMap[140] = 181; stripMap[181] = 140; stripMap[139] = 182; stripMap[182] = 139;
    stripMap[138] = 183; stripMap[183] = 138; stripMap[137] = 184; stripMap[184] = 137;
    stripMap[136] = 185; stripMap[185] = 136; stripMap[135] = 186; stripMap[186] = 135;
    stripMap[134] = 187; stripMap[187] = 134; stripMap[133] = 188; stripMap[188] = 133;
    stripMap[132] = 189; stripMap[189] = 132; stripMap[131] = 190; stripMap[190] = 131;
    stripMap[130] = 191; stripMap[191] = 130; stripMap[129] = 192; stripMap[192] = 129;
    
    stripMap[256] = 193; stripMap[193] = 256; stripMap[255] = 194; stripMap[194] = 255;
    stripMap[254] = 195; stripMap[195] = 254; stripMap[253] = 196; stripMap[196] = 253;
    stripMap[252] = 197; stripMap[197] = 252; stripMap[251] = 198; stripMap[198] = 251;
    stripMap[250] = 199; stripMap[199] = 250; stripMap[249] = 200; stripMap[200] = 249;
    stripMap[248] = 201; stripMap[201] = 248; stripMap[247] = 202; stripMap[202] = 247;
    stripMap[246] = 203; stripMap[203] = 246; stripMap[245] = 204; stripMap[204] = 245;
    stripMap[244] = 205; stripMap[205] = 244; stripMap[243] = 206; stripMap[206] = 243;
    stripMap[242] = 207; stripMap[207] = 242; stripMap[241] = 208; stripMap[208] = 241;
    stripMap[240] = 209; stripMap[209] = 240; stripMap[239] = 210; stripMap[210] = 239;
    stripMap[238] = 211; stripMap[211] = 238; stripMap[237] = 212; stripMap[212] = 237;
    stripMap[236] = 213; stripMap[213] = 236; stripMap[235] = 214; stripMap[214] = 235;
    stripMap[234] = 215; stripMap[215] = 234; stripMap[233] = 216; stripMap[216] = 233;
    stripMap[232] = 217; stripMap[217] = 232; stripMap[231] = 218; stripMap[218] = 231;
    stripMap[230] = 219; stripMap[219] = 230; stripMap[229] = 220; stripMap[220] = 229;
    stripMap[228] = 221; stripMap[221] = 228; stripMap[227] = 222; stripMap[222] = 227;
    stripMap[226] = 223; stripMap[223] = 226; stripMap[225] = 224; stripMap[224] = 225;
    
    stripMap[288] = 289; stripMap[289] = 288; stripMap[287] = 290; stripMap[290] = 287;
    stripMap[286] = 291; stripMap[291] = 286; stripMap[285] = 292; stripMap[292] = 285;
    stripMap[284] = 293; stripMap[293] = 284; stripMap[283] = 294; stripMap[294] = 283;
    stripMap[282] = 295; stripMap[295] = 282; stripMap[281] = 296; stripMap[296] = 281;
    stripMap[280] = 297; stripMap[297] = 280; stripMap[279] = 298; stripMap[298] = 279;
    stripMap[278] = 299; stripMap[299] = 278; stripMap[277] = 300; stripMap[300] = 277;
    stripMap[276] = 301; stripMap[301] = 276; stripMap[275] = 302; stripMap[302] = 275;
    stripMap[274] = 303; stripMap[303] = 274; stripMap[273] = 304; stripMap[304] = 273;
    stripMap[272] = 305; stripMap[305] = 272; stripMap[271] = 306; stripMap[306] = 271;
    stripMap[270] = 307; stripMap[307] = 270; stripMap[269] = 308; stripMap[308] = 269;
    stripMap[268] = 309; stripMap[309] = 268; stripMap[267] = 310; stripMap[310] = 267;
    stripMap[266] = 311; stripMap[311] = 266; stripMap[265] = 312; stripMap[312] = 265;
    stripMap[264] = 313; stripMap[313] = 264; stripMap[263] = 314; stripMap[314] = 263;
    stripMap[262] = 315; stripMap[315] = 262; stripMap[261] = 316; stripMap[316] = 261;
    stripMap[260] = 317; stripMap[317] = 260; stripMap[259] = 318; stripMap[318] = 259;
    stripMap[258] = 319; stripMap[319] = 258; stripMap[257] = 320; stripMap[320] = 257;
    
    stripMap[384] = 321; stripMap[321] = 384; stripMap[383] = 322; stripMap[322] = 383;
    stripMap[382] = 323; stripMap[323] = 382; stripMap[381] = 324; stripMap[324] = 381;
    stripMap[380] = 325; stripMap[325] = 380; stripMap[379] = 326; stripMap[326] = 379;
    stripMap[378] = 327; stripMap[327] = 378; stripMap[377] = 328; stripMap[328] = 377;
    stripMap[376] = 329; stripMap[329] = 376; stripMap[375] = 330; stripMap[330] = 375;
    stripMap[374] = 331; stripMap[331] = 374; stripMap[373] = 332; stripMap[332] = 373;
    stripMap[372] = 333; stripMap[333] = 372; stripMap[371] = 334; stripMap[334] = 371;
    stripMap[370] = 335; stripMap[335] = 370; stripMap[369] = 336; stripMap[336] = 369;
    stripMap[368] = 337; stripMap[337] = 368; stripMap[367] = 338; stripMap[338] = 367;
    stripMap[366] = 339; stripMap[339] = 366; stripMap[365] = 340; stripMap[340] = 365;
    stripMap[364] = 341; stripMap[341] = 364; stripMap[363] = 342; stripMap[342] = 363;
    stripMap[362] = 343; stripMap[343] = 362; stripMap[361] = 344; stripMap[344] = 361;
    stripMap[360] = 345; stripMap[345] = 360; stripMap[359] = 346; stripMap[346] = 359;
    stripMap[358] = 347; stripMap[347] = 358; stripMap[357] = 348; stripMap[348] = 357;
    stripMap[356] = 349; stripMap[349] = 356; stripMap[355] = 350; stripMap[350] = 355;
    stripMap[354] = 351; stripMap[351] = 354; stripMap[353] = 352; stripMap[352] = 353;
    
    stripMap[416] = 417; stripMap[417] = 416; stripMap[415] = 418; stripMap[418] = 415;
    stripMap[414] = 419; stripMap[419] = 414; stripMap[413] = 420; stripMap[420] = 413;
    stripMap[412] = 421; stripMap[421] = 412; stripMap[411] = 422; stripMap[422] = 411;
    stripMap[410] = 423; stripMap[423] = 410; stripMap[409] = 424; stripMap[424] = 409;
    stripMap[408] = 425; stripMap[425] = 408; stripMap[407] = 426; stripMap[426] = 407;
    stripMap[406] = 427; stripMap[427] = 406; stripMap[405] = 428; stripMap[428] = 405;
    stripMap[404] = 429; stripMap[429] = 404; stripMap[403] = 430; stripMap[430] = 403;
    stripMap[402] = 431; stripMap[431] = 402; stripMap[401] = 432; stripMap[432] = 401;
    stripMap[400] = 433; stripMap[433] = 400; stripMap[399] = 434; stripMap[434] = 399;
    stripMap[398] = 435; stripMap[435] = 398; stripMap[397] = 436; stripMap[436] = 397;
    stripMap[396] = 437; stripMap[437] = 396; stripMap[395] = 438; stripMap[438] = 395;
    stripMap[394] = 439; stripMap[439] = 394; stripMap[393] = 440; stripMap[440] = 393;
    stripMap[392] = 441; stripMap[441] = 392; stripMap[391] = 442; stripMap[442] = 391;
    stripMap[390] = 443; stripMap[443] = 390; stripMap[389] = 444; stripMap[444] = 389;
    stripMap[388] = 445; stripMap[445] = 388; stripMap[387] = 446; stripMap[446] = 387;
    stripMap[386] = 447; stripMap[447] = 386; stripMap[385] = 448; stripMap[448] = 385;
    
    stripMap[512] = 449; stripMap[449] = 512; stripMap[511] = 450; stripMap[450] = 511;
    stripMap[510] = 451; stripMap[451] = 510; stripMap[509] = 452; stripMap[452] = 509;
    stripMap[508] = 453; stripMap[453] = 508; stripMap[507] = 454; stripMap[454] = 507;
    stripMap[506] = 455; stripMap[455] = 506; stripMap[505] = 456; stripMap[456] = 505;
    stripMap[504] = 457; stripMap[457] = 504; stripMap[503] = 458; stripMap[458] = 503;
    stripMap[502] = 459; stripMap[459] = 502; stripMap[501] = 460; stripMap[460] = 501;
    stripMap[500] = 461; stripMap[461] = 500; stripMap[499] = 462; stripMap[462] = 499;
    stripMap[498] = 463; stripMap[463] = 498; stripMap[497] = 464; stripMap[464] = 497;
    stripMap[496] = 465; stripMap[465] = 496; stripMap[495] = 466; stripMap[466] = 495;
    stripMap[494] = 467; stripMap[467] = 494; stripMap[493] = 468; stripMap[468] = 493;
    stripMap[492] = 469; stripMap[469] = 492; stripMap[491] = 470; stripMap[470] = 491;
    stripMap[490] = 471; stripMap[471] = 490; stripMap[489] = 472; stripMap[472] = 489;
    stripMap[488] = 473; stripMap[473] = 488; stripMap[487] = 474; stripMap[474] = 487;
    stripMap[486] = 475; stripMap[475] = 486; stripMap[485] = 476; stripMap[476] = 485;
    stripMap[484] = 477; stripMap[477] = 484; stripMap[483] = 478; stripMap[478] = 483;
    stripMap[482] = 479; stripMap[479] = 482; stripMap[481] = 480; stripMap[480] = 481;
    
    stripMap[544] = 545; stripMap[545] = 544; stripMap[543] = 546; stripMap[546] = 543;
    stripMap[542] = 547; stripMap[547] = 542; stripMap[541] = 548; stripMap[548] = 541;
    stripMap[540] = 549; stripMap[549] = 540; stripMap[539] = 550; stripMap[550] = 539;
    stripMap[538] = 551; stripMap[551] = 538; stripMap[537] = 552; stripMap[552] = 537;
    stripMap[536] = 553; stripMap[553] = 536; stripMap[535] = 554; stripMap[554] = 535;
    stripMap[534] = 555; stripMap[555] = 534; stripMap[533] = 556; stripMap[556] = 533;
    stripMap[532] = 557; stripMap[557] = 532; stripMap[531] = 558; stripMap[558] = 531;
    stripMap[530] = 559; stripMap[559] = 530; stripMap[529] = 560; stripMap[560] = 529;
    stripMap[528] = 561; stripMap[561] = 528; stripMap[527] = 562; stripMap[562] = 527;
    stripMap[526] = 563; stripMap[563] = 526; stripMap[525] = 564; stripMap[564] = 525;
    stripMap[524] = 565; stripMap[565] = 524; stripMap[523] = 566; stripMap[566] = 523;
    stripMap[522] = 567; stripMap[567] = 522; stripMap[521] = 568; stripMap[568] = 521;
    stripMap[520] = 569; stripMap[569] = 520; stripMap[519] = 570; stripMap[570] = 519;
    stripMap[518] = 571; stripMap[571] = 518; stripMap[517] = 572; stripMap[572] = 517;
    stripMap[516] = 573; stripMap[573] = 516; stripMap[515] = 574; stripMap[574] = 515;
    stripMap[514] = 575; stripMap[575] = 514; stripMap[513] = 576; stripMap[576] = 513;
    
    stripMap[640] = 577; stripMap[577] = 640; stripMap[639] = 578; stripMap[578] = 639;
    stripMap[638] = 579; stripMap[579] = 638; stripMap[637] = 580; stripMap[580] = 637;
    stripMap[636] = 581; stripMap[581] = 636; stripMap[635] = 582; stripMap[582] = 635;
    stripMap[634] = 583; stripMap[583] = 634; stripMap[633] = 584; stripMap[584] = 633;
    stripMap[632] = 585; stripMap[585] = 632; stripMap[631] = 586; stripMap[586] = 631;
    stripMap[630] = 587; stripMap[587] = 630; stripMap[629] = 588; stripMap[588] = 629;
    stripMap[628] = 589; stripMap[589] = 628; stripMap[627] = 590; stripMap[590] = 627;
    stripMap[626] = 591; stripMap[591] = 626; stripMap[625] = 592; stripMap[592] = 625;
    stripMap[624] = 593; stripMap[593] = 624; stripMap[623] = 594; stripMap[594] = 623;
    stripMap[622] = 595; stripMap[595] = 622; stripMap[621] = 596; stripMap[596] = 621;
    stripMap[620] = 597; stripMap[597] = 620; stripMap[619] = 598; stripMap[598] = 619;
    stripMap[618] = 599; stripMap[599] = 618; stripMap[617] = 600; stripMap[600] = 617;
    stripMap[616] = 601; stripMap[601] = 616; stripMap[615] = 602; stripMap[602] = 615;
    stripMap[614] = 603; stripMap[603] = 614; stripMap[613] = 604; stripMap[604] = 613;
    stripMap[612] = 605; stripMap[605] = 612; stripMap[611] = 606; stripMap[606] = 611;
    stripMap[610] = 607; stripMap[607] = 610; stripMap[609] = 608; stripMap[608] = 609;

    mapInitialized = true;
}
*/

// =================================================================
    // --- EXTRACT NOISE FLOOR BEFORE ENDING THE LOOP ---
    // =================================================================

void mmanalys::Loop(bool useEmpirical)

{
    //   In a ROOT session, you can do:
    //      Root > .L mmanalys.cpp
    //      Root > mmanalys t
    //      Root > t.Loop();       // Loop on all entries
    //          Root > t.GetEntry(12); // Fill t data members with entry number 12
    //          Root > t.Show();       // Show values of entry 12
    //          Root > tt.Show(16);     // Read and show values of entry 16
    //
    gErrorIgnoreLevel = kError;

    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(1111);
    TF1 *trackfit = new TF1("trackfit", "pol1", 0, 1000);
    useEmpiricalSigma = useEmpirical;

    // --- EMPIRICAL RESOLUTION CALIBRATION HISTOGRAMS ---
    // X-axis: Cluster Size (1 to 10), Y-axis: Unbiased Residual (mm)
    static TH2F* h_res_vs_sizeX = nullptr;
    static TH2F* h_res_vs_sizeY = nullptr;
    if (!h_res_vs_sizeX) {
    h_res_vs_sizeX = new TH2F("h_res_vs_sizeX", "Unbiased Residual vs Cluster Size (X);Cluster Size;Residual (mm)", 10, 0.5, 10.5, 100, -2.0, 2.0);
    h_res_vs_sizeY = new TH2F("h_res_vs_sizeY", "Unbiased Residual vs Cluster Size (Y);Cluster Size;Residual (mm)", 10, 0.5, 10.5, 100, -2.0, 2.0);
    } else {
    // Second pass: clear old residuals to get fresh ones
    h_res_vs_sizeX->Reset();
    h_res_vs_sizeY->Reset();
    }

    TH2F *residualstrips1X = new TH2F("residualstrips1X", "Residuals vs Strip Number 1X;Strip Number;Residual (mm)", 640, 0, 640, 200, -100, 100);
    TH2F *residualstrips1Y = new TH2F("residualstrips1Y", "Residuals vs Strip Number 1Y;Strip Number;Residual (mm)", 640, 0, 640, 200, -100, 100);
    TH2F *residualstrips2X = new TH2F("residualstrips2X", "Residuals vs Strip Number 2X;Strip Number;Residual (mm)", 640, 0, 640, 200, -100, 100);
    TH2F *residualstrips2Y = new TH2F("residualstrips2Y", "Residuals vs Strip Number 2Y;Strip Number;Residual (mm)", 640, 0, 640, 200, -100, 100);
    TH2F *residualstrips3X = new TH2F("residualstrips3X", "Residuals vs Strip Number 3X;Strip Number;Residual (mm)", 640, 0, 640, 200, -100, 100);
    TH2F *residualstrips3Y = new TH2F("residualstrips3Y", "Residuals vs Strip Number 3Y;Strip Number;Residual (mm)", 640, 0, 640, 200, -100, 100);
    TH2F *residualstrips4X = new TH2F("residualstrips4X", "Residuals vs Strip Number 4X;Strip Number;Residual (mm)", 640, 0, 640, 200, -100, 100);
    TH2F *residualstrips4Y = new TH2F("residualstrips4Y", "Residuals vs Strip Number 4Y;Strip Number;Residual (mm)", 640, 0, 640, 200, -100, 100);

    TH2F *h_stripXY_X1X2 = new TH2F("h_stripXY_X1X2", "Strip Map X1 vs X2 ;Layer 1X Strip;Layer 2X Strip", 640, 0, 640, 640, 0, 640);
    TH2F *h_stripXY_X3X4 = new TH2F("h_stripXY_X3X4", "Strip Map X3 vs X4 ;Layer 3X Strip;Layer 4X Strip", 640, 0, 640, 640, 0, 640);
    // New histograms for angleX using only centroid and only uTPC
    TH1F *h_angleX_centroid = new TH1F("h_angleX_centroid", "Track Angle X (centroid only);Angle X (deg);Counts", 180, 60, 120);
    TH1F *h_angleX_utpc = new TH1F("h_angleX_utpc", "Track Angle X (uTPC only);Angle X (deg);Counts", 180, 60, 120);
    TH1F *h_angleY_centroid = new TH1F("h_angleY_centroid", "Track Angle Y (centroid only);Angle Y (deg);Counts", 180, 60, 120);
    TH1F *h_angleY_utpc = new TH1F("h_angleY_utpc", "Track Angle Y (uTPC only);Angle Y (deg);Counts", 180, 60, 120);

    // ---> TH2F: CHI^2 vs STRIP NUMBER <---
    TH2F *h_chi2_vs_strip_1X = new TH2F("h_chi2_vs_strip_1X", "#chi^{2} vs Strip 1X;Strip #;Track #chi^{2}", 640, 0, 640, 1000, 0, 10000);
    TH2F *h_chi2_vs_strip_2X = new TH2F("h_chi2_vs_strip_2X", "#chi^{2} vs Strip 2X;Strip #;Track #chi^{2}", 640, 0, 640, 1000, 0, 10000);
    TH2F *h_chi2_vs_strip_3X = new TH2F("h_chi2_vs_strip_3X", "#chi^{2} vs Strip 3X;Strip #;Track #chi^{2}", 640, 0, 640, 1000, 0, 10000);
    TH2F *h_chi2_vs_strip_4X = new TH2F("h_chi2_vs_strip_4X", "#chi^{2} vs Strip 4X;Strip #;Track #chi^{2}", 640, 0, 640, 1000, 0, 10000);

    TH2F *h_chi2_vs_strip_1Y = new TH2F("h_chi2_vs_strip_1Y", "#chi^{2} vs Strip 1Y;Strip #;Track #chi^{2}", 640, 0, 640, 1000, 0, 10000);
    TH2F *h_chi2_vs_strip_2Y = new TH2F("h_chi2_vs_strip_2Y", "#chi^{2} vs Strip 2Y;Strip #;Track #chi^{2}", 640, 0, 640, 1000, 0, 10000);
    TH2F *h_chi2_vs_strip_3Y = new TH2F("h_chi2_vs_strip_3Y", "#chi^{2} vs Strip 3Y;Strip #;Track #chi^{2}", 640, 0, 640, 1000, 0, 10000);
    TH2F *h_chi2_vs_strip_4Y = new TH2F("h_chi2_vs_strip_4Y", "#chi^{2} vs Strip 4Y;Strip #;Track #chi^{2}", 640, 0, 640, 1000, 0, 10000);

// ---> CREATE HISTOGRAMS TO TEST SHIFTS FROM -0.5 to +0.5 <---
    TH1F *scan_hist[100];
    for(int i=0; i<=99; i++) {
        std::string name = "scan_hist_" + std::to_string(i);
        scan_hist[i] = new TH1F(name.c_str(), "Scan;Residual;Counts", 200, -50, 50);
    }

    // X: Deviation from 90, Y: Chi2/ndf, Z (Color): Average Residual
    TProfile2D *p_master_graph = new TProfile2D("p_master_graph", 
    "Absolute Residual vs |Angle X - 90#circ| and #chi^{2}/ndf; |Angle X - 90#circ| (deg); #chi^{2}/ndf", 
    20, 0, 20, 500, 0, 10000);

// (Optional) A Profile to see the AVERAGE Chi^2 for each range
TProfile *p_avgChi2_vs_deltaAngleX = new TProfile("p_avgChi2_vs_deltaAngleX", 
    "Average #chi^{2} vs |Angle X - 90#circ|; |Angle X - 90#circ| (deg); Mean #chi^{2}", 
    20, 0, 20);

// --- Ιστογράμματα για την Κατανομή Landau (Συνολικό Φορτίο Cluster) ---
TH1F *h_landau_charge_X1 = new TH1F("h_landau_charge_X1", "Cluster Charge Landau X1;Total Charge;Counts", 500, 0, 5000);
TH1F *h_landau_charge_X2 = new TH1F("h_landau_charge_X2", "Cluster Charge Landau X2;Total Charge;Counts", 500, 0, 5000);
TH1F *h_landau_charge_X3 = new TH1F("h_landau_charge_X3", "Cluster Charge Landau X3;Total Charge;Counts", 500, 0, 5000);
TH1F *h_landau_charge_X4 = new TH1F("h_landau_charge_X4", "Cluster Charge Landau X4;Total Charge;Counts", 500, 0, 5000);

TH1F *h_landau_charge_Y1 = new TH1F("h_landau_charge_Y1", "Cluster Charge Landau Y1;Total Charge;Counts", 500, 0, 5000);
TH1F *h_landau_charge_Y2 = new TH1F("h_landau_charge_Y2", "Cluster Charge Landau Y2;Total Charge;Counts", 500, 0, 5000);
TH1F *h_landau_charge_Y3 = new TH1F("h_landau_charge_Y3", "Cluster Charge Landau Y3;Total Charge;Counts", 500, 0, 5000);
TH1F *h_landau_charge_Y4 = new TH1F("h_landau_charge_Y4", "Cluster Charge Landau Y4;Total Charge;Counts", 500, 0, 5000);


    //  HBooking  from included file
    TH1F *h_singlestrip_1X= new TH1F("h_singlestrip_1X", ";strip#",641,0,641);
    TH1F *h_singlestrip_2X= new TH1F("h_singlestrip_2X", ";strip#",641,0,641);
    TH1F *h_singlestrip_3X= new TH1F("h_singlestrip_3X", ";strip#",641,0,641);
    TH1F *h_singlestrip_4X= new TH1F("h_singlestrip_4X", ";strip#",641,0,641);
    TH1F *h_singlestrip_1Y= new TH1F("h_singlestrip_1Y", ";strip#",641,0,641);
    TH1F *h_singlestrip_2Y= new TH1F("h_singlestrip_2Y", ";strip#",641,0,641);
    TH1F *h_singlestrip_3Y= new TH1F("h_singlestrip_3Y", ";strip#",641,0,641);
    TH1F *h_singlestrip_4Y= new TH1F("h_singlestrip_4Y", ";strip#",641,0,641);

    TH1F *h_singlestripcluster_1X= new TH1F("h_singlestripcluster_1X", ";strip#",640,0,641);
    TH1F *h_singlestripcluster_2X= new TH1F("h_singlestripcluster_2X", ";strip#",640,0,641);
    TH1F *h_singlestripcluster_3X= new TH1F("h_singlestripcluster_3X", ";strip#",640,0,641);
    TH1F *h_singlestripcluster_4X= new TH1F("h_singlestripcluster_4X", ";strip#",640,0,641);

    TH1F *h_stripschi2_1X= new TH1F("h_stripschi2_1X", ";strip#",640,0,641);
    TH1F *h_stripschi2_2X= new TH1F("h_stripschi2_2X", ";strip#",640,0,641);
    TH1F *h_stripschi2_3X= new TH1F("h_stripschi2_3X", ";strip#",640,0,641);
    TH1F *h_stripschi2_4X= new TH1F("h_stripschi2_4X", ";strip#",640,0,641);
    TH1F *h_stripschi2_1Y= new TH1F("h_stripschi2_1Y", ";strip#",640,0,641);
    TH1F *h_stripschi2_2Y= new TH1F("h_stripschi2_2Y", ";strip#",640,0,641);
    TH1F *h_stripschi2_3Y= new TH1F("h_stripschi2_3Y", ";strip#",640,0,641);
    TH1F *h_stripschi2_4Y= new TH1F("h_stripschi2_4Y", ";strip#",640,0,641);

TH2F *h_stripXY_X1Y1cl0 = new TH2F("h_stripXY_X1Y1cl0", "X1 vs Y1 cl0 ;X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_X2Y2cl0 = new TH2F("h_stripXY_X2Y2cl0", "X2 vs Y2 cl0 ;X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_X3Y3cl0 = new TH2F("h_stripXY_X3Y3cl0", "X3 vs Y3 cl0 ;X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_X4Y4cl0 = new TH2F("h_stripXY_X4Y4cl0", "X4 vs Y4 cl0 ;X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);
/*
TH2F *h_stripXY_X1Y1cl1 = new TH2F("h_stripXY_X1Y1cl1", "X1 vs Y1 cl1 ;X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_X2Y2cl1 = new TH2F("h_stripXY_X2Y2cl1", "X2 vs Y2 cl1 ;X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_X3Y3cl1 = new TH2F("h_stripXY_X3Y3cl1", "X3 vs Y3 cl1 ;X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_X4Y4cl1 = new TH2F("h_stripXY_X4Y4cl1", "X4 vs Y4 cl1 ;X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);

TH2F *h_stripXY_X1Y1cl2 = new TH2F("h_stripXY_X1Y1cl2", "X1 vs Y1 cl2 ;X1 strip;Y1 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_X2Y2cl2 = new TH2F("h_stripXY_X2Y2cl2", "X2 vs Y2 cl2 ;X2 strip;Y2 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_X3Y3cl2 = new TH2F("h_stripXY_X3Y3cl2", "X3 vs Y3 cl2 ;X3 strip;Y3 strip", 640, 0, 640, 640, 0, 640);
TH2F *h_stripXY_X4Y4cl2 = new TH2F("h_stripXY_X4Y4cl2", "X4 vs Y4 cl2 ;X4 strip;Y4 strip", 640, 0, 640, 640, 0, 640);
*/

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

/*
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
*/
/*
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
*/
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
/*
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
*/

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

/*
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
*/

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
/*
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
*/
// --- Slope Analysis Histograms ---
TH1F *h_slope_X = new TH1F("h_slope_X", "Track Slope (X direction);Slope (a);Counts", 151, -75, 75);
TH1F *h_slope_Y = new TH1F("h_slope_Y", "Track Slope (Y direction);Slope (a);Counts", 151, -75, 75);
TH1F *h_intercept_X = new TH1F("h_intercept_X", "Track Intercept (X direction);Intercept (b) [mm];Counts", 9500, -4500, 5500);
TH1F *h_intercept_Y = new TH1F("h_intercept_Y", "Track Intercept (Y direction);Intercept (b) [mm];Counts", 9500, -4500, 5500);
TH2F *h_slope_vs_intercept_X = new TH2F("h_slope_vs_intercept_X", "Slope vs Intercept (X);Slope (a);Intercept (b)", 151, -75, 75, 9500, -4500, 5500);
TH2F *h_slope_vs_intercept_Y = new TH2F("h_slope_vs_intercept_Y", "Slope vs Intercept (Y);Slope (a);Intercept (b)", 151, -75, 75, 9500, -4500, 5500);

TH1F *h_slopeX_chi2lt = new TH1F("h_slopeX_chi2lt", "Track Angle X (#chi^{2});Angle X (deg);Counts", 180, 60, 120);
TH1F *h_slopeY_chi2lt = new TH1F("h_slopeY_chi2lt", "Track Angle Y (#chi^{2});Angle Y (deg);Counts", 180, 60, 120);
TH2F *h_slopeXY_chi2lt = new TH2F("h_slopeXY_chi2lt", "Track Angles X vs Y (#chi^{2});Angle X (deg);Angle Y (deg)", 180, 60, 120, 180, 60, 120);
TH1F *h_track_chi2 = new TH1F("h_track_chi2", "Track #chi^{2};#chi^{2};Counts", 200, 0, 5000);
    enum chambers{ ARTEMIS_X, APOLLON_Y, CLEO_X, DIONYSOS_Y, BEATRICE_X, CASTOR_Y, DAPHNE_X, BYRON_Y };
    map<string, int> enumResolver; 


    enumResolver["ARTEMIS_X"] = ARTEMIS_X;
    enumResolver["APOLLON_Y"] = APOLLON_Y;
    enumResolver["CLEO_X"] = CLEO_X;
    enumResolver["DIONYSOS_Y"] = DIONYSOS_Y;
    enumResolver["BEATRICE_X"] = BEATRICE_X;
    enumResolver["CASTOR_Y"] = CASTOR_Y;
    enumResolver["DAPHNE_X"] = DAPHNE_X;
    enumResolver["BYRON_Y"] = BYRON_Y;

    vector<MMHit> HT1X, HT2X, HT3X, HT4X;
    vector<MMHit> HT1Y, HT2Y, HT3Y, HT4Y;
    MMHit  hit1;
    

    MMCluster C1X, C2X, C3X, C4X, C1Y, C2Y, C3Y, C4Y;
 
    MMCluster CΜΧ1X, CΜΧ2X, CΜΧ3X, CΜΧ4X, CΜΧ1Y, CΜΧ2Y, CΜΧ3Y, CΜΧ4Y;   // leading cluster with max charge

    vector<MMCluster*> VClustersX1, VClustersX2, VClustersX3, VClustersX4;
    vector<MMCluster*> VClustersY1, VClustersY2, VClustersY3, VClustersY4;

    MMPacmanAlgo  clusterAlgo;
    
 

    
    int   event_typeX=0, event_typeY=0  ;
    int   clm1X=0, clm2X=0,clm3X=0, clm4X=0;
    int   clm1Y=0, clm2Y=0,clm3Y=0, clm4Y=0;

    int   clpev1x=0,  clpev2x=0, clpev3x=0, clpev4x=0;
    int   clpev1y=0,  clpev2y=0, clpev3y=0, clpev4y=0;
    int   leadcl=0;
//  corr_pos2X= -0.3718,  corr_pos3X= -0.1
//    float       corr_pos1X= 0.164+0.08519,   corr_pos2X= 0.11-0.01605 ,  corr_pos3X= -0.075-0.2096  ,       corr_pos4X = 0.02+0.2615;
//    float       corr_pos1X= 0.425,   corr_pos2X= -0.668,  corr_pos3X= 0.217 ,       corr_pos4X = -0.002; //.350;
//float       corr_pos1X= 0.164,   corr_pos2X= 0.0,  corr_pos3X= -0.344014 ,       corr_pos4X = 0.02;
float       corr_pos1Y= -0.410 ,  corr_pos2Y= 0.100 , corr_pos3Y= 0.30 ,     corr_pos4Y = -0.260;
//    float       corr_pos1X= -5.0752,   corr_pos2X= 4.118,  corr_pos3X= 5.694 ,       corr_pos4X = -4.07339; //.350;
//    float       corr_pos1Y= 0.020 ,  corr_pos2Y= 0.046 , corr_pos3Y= -0.154 ,     corr_pos4Y = 0.086;
//    float       pos1Z= 210.4,        pos2Z=  139.677,     pos3Z=  70.326,    pos4Z = 0.5;
//float       corr_pos1Y= 0.0,   corr_pos2Y= 0.00,  corr_pos3Y= 0.00 ,       corr_pos4Y = 0.00;
//float       corr_pos1X= 0.0,   corr_pos2X= 0.00,  corr_pos3X= 0.00 ,       corr_pos4X = 0.00;
float corr_pos1X = 0.359, corr_pos2X = -0.487, corr_pos3X = 0.237, corr_pos4X = 0.034;
//float       corr_pos1Y= 0.00 ,  corr_pos2Y= 0.000 , corr_pos3Y= 0.00 ,     corr_pos4Y = 0.00;
    float       pos1Z= 1000.0,        pos2Z=  750.0,     pos3Z=  250.0,    pos4Z = 0.0;
    float       posZ=0;

    float pos1X=0, pos2X=0, pos3X=0, pos4X=0 ;
    float pos1Y=0, pos2Y=0, pos3Y=0, pos4Y=0 ;

    int leadclX1=0, leadclX2=0, leadclX3=0, leadclX4=0;
    
    //float t01=200.0, t02=200.0, t03=200.0, t04=200.0;
    float t01=180.1, t02=187.8, t03=185.0, t04=189.6;
    //float t01=25.0, t02=25.0, t03=25.0, t04=25.0;
    float tmx1=226.971,tmx2=263.993,tmx3=270.719,tmx4= 284.496;
    float v1=0.047, v2=0.047, v3 = 0.047, v4=0.047;        //drift velocity in mm/nanosec
    float v12=0.047; // average drift velocity in mm/nanosec
 
    float mapX1[400000], mapY1[400000];
    int mi1=0;
   
    
    int countstrips=0;
    
    float algnX1,algnX2, algnY1, algnY2;
    float xef=0, zef=0;
  
    int  eff2_all=0, eff2_hit=0, eff3_all=0, eff3_hit=0 ;
    float eff3, eff4,eff2, eff1;
    int eff3_hit1234=0, eff3_hit124=0, eff3_hit134=0, eff3_hit14=0, eff3_hit234=0, eff3_hit24=0;
    
    float coefX=0, intrcX=0, coefY=0, intrcY=0,a=0, b=0, angleX=0, angleY=0, angle24=0,angle14=0, angle12=0, angle23=0;
    float slope14=0,  slope24=0, slope12=0, slope23=0;
    float ang1=0, ang2=0, ang3=0, ang4=0;
    
    float xpos[20], ypos[20], zpos[20];
    float residual1, residual2, residual3, residual4;
    float Rp2=0,  pp2=0;
    
    float a1gr=0, a0gr=0, sa1=0, sa0=0, a1=0, a0=0, as0=0, as1=0, aa1=0, aa0=0, aa1er=0, aa0er=0, chi2=0, ref_x2=0, ref_ndf=0, ref_Px2;
    float pitch1=0.6;                //pitch width in mm
    float chisqr;
    float chi2X=0, chi2Y=0;
    int entries=0;
    char  chr;
    float prj_X=0;
    
    // CUTS
    int cl_minstrp=0 ,cl_maxstrp=10;  // min and max cluster size (in strips)

    
    TH1F *h_event_cluster= new TH1F("h_event_cluster", ";strips",150,0,150);

    // ==========================================
    // --- ML DATA BRIDGE: OPEN CSV EXPORT ---
    // ==========================================
    std::ofstream ml_data;
    if (useEmpirical) {
        ml_data.open("ml_real_data.csv");
        // Write the header row for Python to read
        ml_data << "EventID,X1,Y1,Z1,Q1,Size1,"
                << "X2,Y2,Z2,Q2,Size2,"
                << "X3,Y3,Z3,Q3,Size3,"
                << "X4,Y4,Z4,Q4,Size4,"
                << "AngleX,AngleY\n";
    }

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    std::cout << nentries << "  =  N entries"  << std::endl;
    Long64_t nbytes = 0, nb = 0;
// for (Long64_t jentry=0; jentry< 10000 ;jentry++) {
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
     
        /*
        int nstrips  = mmStrip->size();
        h_event_cluster->Fill(nstrips);
    
      for(int k=0; k<nstrips; k++){
          hit1.SetStrip(mmStrip->at(k));
          hit1.sety(((FD_time->at(k))-t01)*v12);
          hit1.SetCharge(max_q->at(k));
          hit1.setyer(max_q->at(k));
          switch(enumResolver[mmChamber->at(k)])
          {
              case 0 :
                  C1X.AddHit(hit1);                break;
              case 1 :
                  C1Y.AddHit(hit1);                break;
              case 2 :
                  C2X.AddHit(hit1);                break;
              case 3 :
                  C2Y.AddHit(hit1);                break;
              case 4 :
                  C3X.AddHit(hit1);                break;
              case 5 :
                  C3Y.AddHit(hit1);                break;
              case 6 :
                  C4X.AddHit(hit1);                break;
              case 7 :
                  C4Y.AddHit(hit1);                break;
              default : ;
          }
       } //  loop over the strips in each event
*/
    
        int nstrips  = mmStrip->size();
        h_event_cluster->Fill(nstrips);
        
   for(int k=0; k<nstrips; k++){
          float current_charge = max_q->at(k);
          
          //if (current_charge > max_allowed_charge) {
          //  current_charge = max_allowed_charge;
          //}

          // 1. Call your initialization function (it runs once)
          if (!mapInitialized) initializeMap();

          // 2. Lookup the corrected strip number
          int raw_strip = mmStrip->at(k);
          int mapped_strip = raw_strip;

          // Pull the right value from your global array
          if (raw_strip >= 0 && raw_strip < 800) {
              mapped_strip = stripMap[raw_strip];
          }

          
          int chamber_idx = enumResolver[mmChamber->at(k)];
          int actual_threshold = 80 ;
          if (current_charge <= actual_threshold) {
            continue; 
          }
          // 3. Apply Time, Charge, and Y-coordinate scaling
          hit1.SetStrip(mapped_strip);        
          hit1.sety(((FD_time->at(k)) - t01) * v12);
          hit1.SetCharge(current_charge);
          hit1.SetTime(FD_time->at(k));
          hit1.setyer(current_charge);

          
          switch(chamber_idx)
          {
              case 0 :  C1X.AddHit(hit1); break;
              case 1 :  C1Y.AddHit(hit1); break;
              case 2 :  C2X.AddHit(hit1); break;
              case 3 :  C2Y.AddHit(hit1); break;
              case 4 :  C3X.AddHit(hit1); break;
              case 5 :  C3Y.AddHit(hit1); break;
              case 6 :  C4X.AddHit(hit1); break;
              case 7 :  C4Y.AddHit(hit1); break;
              default : ;
          }
       };
        /*
// Define your charge cap to suppress delta electrons
      float max_allowed_charge = 10000.0; 

      for(int k=0; k<nstrips; k++){
          float current_charge = max_q->at(k);
          
          // Truncate (cap) the charge if it exceeds the limit
          if (current_charge > max_allowed_charge) {
              current_charge = max_allowed_charge;
          }

          hit1.SetStrip(mmStrip->at(k));
          hit1.sety(((FD_time->at(k))-t01)*v12);
          hit1.SetCharge(current_charge);
          hit1.SetTime(FD_time->at(k));
          hit1.setyer(current_charge);
          
          int chamber_idx = enumResolver[mmChamber->at(k)];

          switch(chamber_idx)
          {
              case 0 :  C1X.AddHit(hit1); break;
              case 1 :  C1Y.AddHit(hit1); break;
              case 2 :  C2X.AddHit(hit1); break;
              case 3 :  C2Y.AddHit(hit1); break;
              case 4 :  C3X.AddHit(hit1); break;
              case 5 :  C3Y.AddHit(hit1); break;
              case 6 :  C4X.AddHit(hit1); break;
              case 7 :  C4Y.AddHit(hit1); break;
              default : ;
          }
       }
        */  
         /*
         for(int k=0; k<nstrips; k++){
          
          // --- APV HARDWARE-TO-SOFTWARE MAPPING FIX ---
          int raw_strip = mmStrip->at(k);
          int mapped_strip = raw_strip;

          // If the strip is in the LAST APV (512-639), shift it to the FIRST APV (0-127)
          if (raw_strip >= 512 && raw_strip <= 639) {
              mapped_strip = raw_strip - 512; 
          } 
          // If the strip is in the FIRST APV (0-127), shift it to the LAST APV (512-639)
          else if (raw_strip >= 0 && raw_strip <= 127) {
              mapped_strip = raw_strip + 512;
          }

          hit1.SetStrip(mapped_strip); // Feed the fixed strip into the reconstruction!
          // --------------------------------------------

          hit1.sety(((FD_time->at(k))-t01)*v12);
          hit1.SetCharge(max_q->at(k));
          hit1.SetTime(FD_time->at(k));
        
    */
// Define your charge cap to suppress delta electrons
      //float max_allowed_charge = 1600.0; 
/*
      for(int k=0; k<nstrips; k++){
          float current_charge = max_q->at(k);
          
          // Truncate (cap) the charge if it exceeds the limit
          if (current_charge > max_allowed_charge) {
              current_charge = max_allowed_charge;
          }

          hit1.SetStrip(mmStrip->at(k));
          hit1.sety(((FD_time->at(k))-t01)*v12);
          hit1.SetCharge(current_charge);
          hit1.SetTime(FD_time->at(k));
          hit1.setyer(current_charge);
          */
/*
        // ==========================================================
          // --- FULL APV BLOCK SWAP (EXPLICIT STRIP RANGES) ---
          // ==========================================================
          int raw_strip = mmStrip->at(k);
          int mapped_strip = raw_strip;
        
          // 1. Swap APV 0 (0 to 127) with APV 4 (512 to 639)
          if (raw_strip >= 0 && raw_strip <= 127) {
              mapped_strip = raw_strip + 512;
          } 
          else if (raw_strip >= 512 && raw_strip <= 639) {
              mapped_strip = raw_strip - 512;
          }
          
          // 2. Swap APV 1 (128 to 255) with APV 3 (384 to 511)
          else if (raw_strip >= 128 && raw_strip <= 255) {
              mapped_strip = raw_strip + 256;
          } 
          else if (raw_strip >= 384 && raw_strip <= 511) {
              mapped_strip = raw_strip - 256;
          }
          
          // Note: APV 2 (256 to 383) is the middle APV. 
          // It does not trigger any of these 'if' statements, 
          // so it stays perfectly in place!

          hit1.SetStrip(mapped_strip); 
          // ==========================================================
*/
          /*
int chamber_idx = enumResolver[mmChamber->at(k)];
          int raw_strip = mmStrip->at(k);
          int mapped_strip = raw_strip;
          
          
          // ==========================================================
          // --- SURGICAL APV MAPPING FIX (ALL LAYERS) ---
          // ==========================================================

          // --- LAYER 1X (chamber_idx == 0) ---
          if (chamber_idx == 0) {
              if (raw_strip >= 0 && raw_strip <= 63) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 64 && raw_strip <= 127) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 128 && raw_strip <= 191) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 192 && raw_strip <= 255) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 256 && raw_strip <= 319) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 320 && raw_strip <= 383) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 384 && raw_strip <= 447) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 448 && raw_strip <= 511) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 512 && raw_strip <= 575) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 576 && raw_strip <= 639) mapped_strip = raw_strip - 64;
          }
          // --- LAYER 2X (chamber_idx == 2) ---
            else if (chamber_idx == 2) {
              if (raw_strip >= 0 && raw_strip <= 63) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 64 && raw_strip <= 127) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 128 && raw_strip <= 191) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 192 && raw_strip <= 255) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 256 && raw_strip <= 319) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 320 && raw_strip <= 383) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 384 && raw_strip <= 447) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 448 && raw_strip <= 511) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 512 && raw_strip <= 575) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 576 && raw_strip <= 639) mapped_strip = raw_strip - 64;
          }
          // --- LAYER 3X (chamber_idx == 4) ---
          else if (chamber_idx == 4) {
              if (raw_strip >= 0 && raw_strip <= 63) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 64 && raw_strip <= 127) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 128 && raw_strip <= 191) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 192 && raw_strip <= 255) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 256 && raw_strip <= 319) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 320 && raw_strip <= 383) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 384 && raw_strip <= 447) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 448 && raw_strip <= 511) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 512 && raw_strip <= 575) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 576 && raw_strip <= 639) mapped_strip = raw_strip - 64;   
          }
          // --- LAYER 4X (chamber_idx == 6) ---
          else if (chamber_idx == 6) {
              if (raw_strip >= 0 && raw_strip <= 63) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 64 && raw_strip <= 127) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 128 && raw_strip <= 191) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 192 && raw_strip <= 255) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 256 && raw_strip <= 319) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 320 && raw_strip <= 383) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 384 && raw_strip <= 447) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 448 && raw_strip <= 511) mapped_strip = raw_strip - 64;
              else if (raw_strip >= 512 && raw_strip <= 575) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 576 && raw_strip <= 639) mapped_strip = raw_strip - 64;
          }
          
          hit1.SetStrip(mapped_strip); */
          // ==========================================================
        //
          // Map fix attempt succeded6
        
          /*
          int chamber_idx = enumResolver[mmChamber->at(k)];
          
          // ==========================================================
          // --- SURGICAL APV MAPPING FIX (ALL LAYERS) ---
          // ==========================================================
          int raw_strip = mmStrip->at(k);
          int mapped_strip = raw_strip;
          
          // --- LAYER 1X (chamber_idx == 0) ---
          if (chamber_idx == 0) {
              if (raw_strip >= 512 && raw_strip <= 575) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 576 && raw_strip <= 639) mapped_strip = raw_strip - 64;
          }
          // --- LAYER 2X (chamber_idx == 2) ---
          else if (chamber_idx == 2) {
              // 0-64 <-> 64-128
              if (raw_strip >= 0 && raw_strip <= 63) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 64 && raw_strip <= 127) mapped_strip = raw_strip - 64;
              // 384-448 <-> 448-512
              else if (raw_strip >= 384 && raw_strip <= 447) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 448 && raw_strip <= 511) mapped_strip = raw_strip - 64;
              // 512-576 <-> 576-640
              else if (raw_strip >= 512 && raw_strip <= 575) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 576 && raw_strip <= 639) mapped_strip = raw_strip - 64;
          }
          // --- LAYER 3X (chamber_idx == 4) ---
          else if (chamber_idx == 4) {
              // 0-64 <-> 64-128
              if (raw_strip >= 0 && raw_strip <= 63) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 64 && raw_strip <= 127) mapped_strip = raw_strip - 64;
              // 384-448 <-> 448-512
              else if (raw_strip >= 384 && raw_strip <= 447) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 448 && raw_strip <= 511) mapped_strip = raw_strip - 64;
              // 512-576 <-> 576-640
              else if (raw_strip >= 512 && raw_strip <= 575) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 576 && raw_strip <= 639) mapped_strip = raw_strip - 64;
          }
          // --- LAYER 4X (chamber_idx == 6) ---
          else if (chamber_idx == 6) {
              // 512-576 <-> 576-640
              if (raw_strip >= 512 && raw_strip <= 575) mapped_strip = raw_strip + 64;
              else if (raw_strip >= 576 && raw_strip <= 639) mapped_strip = raw_strip - 64;
          }
          
          hit1.SetStrip(mapped_strip); 
          // ==========================================================
*/
          /*
        int chamber_idx = enumResolver[mmChamber->at(k)];
        float v_drift = 0.047; // CORRECTED DRIFT VELOCITY
        // ==========================================
        // --- SURGICAL APV MAPPING FIX ---
        // ==========================================
        int raw_strip = mmStrip->at(k);
        int mapped_strip = raw_strip;
          
        // Layer 1X (0) and Layer 4X (6) have APV 1 and APV 4 swapped
        if (chamber_idx == 0 || chamber_idx == 6) { 
            if (raw_strip >= 128 && raw_strip <= 255) {
                mapped_strip = raw_strip + 384; 
            } 
            else if (raw_strip >= 512 && raw_strip <= 639) {
                mapped_strip = raw_strip - 384; 
            }
        }
        
        hit1.SetStrip(mapped_strip); 
        // ==========================================

        // Find the correct t0 for each specific chamber
        float t0_correct = t01;
        if (chamber_idx == 2 || chamber_idx == 3) t0_correct = t02;      
        else if (chamber_idx == 4 || chamber_idx == 5) t0_correct = t03; 
        else if (chamber_idx == 6 || chamber_idx == 7) t0_correct = t04; 

        // Apply proper Y coordinate scaling using correct t0 and velocity
        hit1.sety(((FD_time->at(k)) - t0_correct) * v_drift);
        hit1.SetCharge(current_charge);
        hit1.SetTime(FD_time->at(k));
        hit1.setyer(current_charge);
        */
          // I thought I had fixed the chamber mapping issue, but it seems to persist. this code was just a test to see if the mapping was correct, but the strip distributions still look wrong. I will need to investigate further.
          //int chamber_idx = enumResolver[mmChamber->at(k)];
          
          // ==========================================
          // --- MAPPING FIX: REVERSE LAYER 2 & 4 ---
          // ==========================================
          //int raw_strip = mmStrip->at(k);
          //int mapped_strip = raw_strip;
          /*
          // chamber_idx 0 is ARTEMIS_X (1X),chamber_idx 2 is BEATRICE_X (2X), chamber_idx 4 is CLEO_X(3X), chamber_idx 6 is DAPHNE_X (4X)
            if (chamber_idx == 2 || chamber_idx == 4) { // If it's layer 2 or 4, reverse the strip numbering
                // Assuming strips are numbered 1 to 640. If they are 0 to 639, use (639 - raw_strip)
                mapped_strip = 641 - raw_strip; 
          }


          hit1.SetStrip(mapped_strip); 
          // ==========================================
          */
         // V2 of the mapping fix: reverse the order of the APV cables for layers 1X which are the ones that look wrong. This is a more complex mapping than just reversing the strip numbers, because each APV has 128 strips and they are swapped in groups of 128.
        /* int chamber_idx = enumResolver[mmChamber->at(k)];
          
          // ==========================================
          // --- MAPPING FIX: REVERSE APV CABLE ORDER ---
          // ==========================================
          int raw_strip = mmStrip->at(k);
          int mapped_strip = raw_strip;
          
          // Assuming Layer 1X (0) and Layer 4X (6) have the swapped cables
          if (chamber_idx == 0 || chamber_idx == 6) { 
              // 1. Find out which of the 5 APVs this strip belongs to (0, 1, 2, 3, or 4)
              int apv_num = raw_strip / 128; 
              
              // 2. Find the strip's local position inside that specific APV (0 to 127)
              int strip_in_apv = raw_strip % 128; 
              
              // 3. Mirror the APV order (0 becomes 4, 1 becomes 3, 2 stays 2, etc.)
              int new_apv_num = 4 - apv_num; 
              
              // 4. Reconstruct the global strip number
              mapped_strip = (new_apv_num * 128) + strip_in_apv; 
          }
          
          hit1.SetStrip(mapped_strip); 
          // ==========================================*/
         /*
         // V3 of the mapping fix: swap only APV 0 and APV 4 for layers 1X and 4X, leaving APVs 1, 2, and 3 unchanged. This is based on the observation that only the first and last APVs seem to be swapped in the data, while the middle ones are correct.
          int chamber_idx = enumResolver[mmChamber->at(k)];
          
          // ==========================================
          // --- MAPPING FIX: SWAP ONLY APV 0 AND APV 4 ---
          // ==========================================
          int raw_strip = mmStrip->at(k);
          int mapped_strip = raw_strip;
          
          // Apply only to Layer 1X (0) and Layer 4X (6)
          if (chamber_idx == 0 || chamber_idx == 6) { 
              if (raw_strip >= 512 && raw_strip <= 639) {
                  // Move the last APV to the first slot
                  mapped_strip = raw_strip - 512; 
              } 
              else if (raw_strip >= 0 && raw_strip <= 127) {
                  // Move the first APV to the last slot
                  mapped_strip = raw_strip + 512;
              }
              // Notice we do NOTHING to strips between 128 and 511. 
              // They stay exactly where they are!
          }
          
          hit1.SetStrip(mapped_strip); 
          // ==========================================
         */
        /*
        //V4 of the mapping fix: swap APV 0 with APV 4 for layer 1X, and swap APV 1 with APV 3 for layer 4X, based on the specific patterns observed in the data. This is a more surgical fix that targets only the APVs that seem to be swapped in each layer, rather than applying a blanket reversal to all APVs in those layers.
        int chamber_idx = enumResolver[mmChamber->at(k)];
          
          // ==========================================
          // --- SURGICAL APV MAPPING FIX ---
          // ==========================================
          int raw_strip = mmStrip->at(k);
          int mapped_strip = raw_strip;
          
          if (chamber_idx == 0) { // LAYER 1X
              // Swapped: APV 1 (128-255) and APV 4 (512-639)
              if (raw_strip >= 128 && raw_strip <= 255) {
                  mapped_strip = raw_strip + 384; 
              } else if (raw_strip >= 512 && raw_strip <= 639) {
                  mapped_strip = raw_strip - 384; 
              }
          } 
          else if (chamber_idx == 6) { // LAYER 4X
              // Swapped: APV 1 (128-255) and APV 3 (384-511)
              if (raw_strip >= 128 && raw_strip <= 255) {
                  mapped_strip = raw_strip + 256; 
              } else if (raw_strip >= 384 && raw_strip <= 511) {
                  mapped_strip = raw_strip - 256; 
              }
          }
          
          hit1.SetStrip(mapped_strip); 
          // ==========================================
          */
         
          /*
         for(int k=0; k<nstrips; k++){
          
          // --- APV HARDWARE-TO-SOFTWARE MAPPING FIX ---
          int raw_strip = mmStrip->at(k);
          int mapped_strip = raw_strip;

          // If the strip is in the LAST APV (512-639), shift it to the FIRST APV (0-127)
          if (raw_strip >= 512 && raw_strip <= 639) {
              mapped_strip = raw_strip - 512; 
          } 
          // If the strip is in the FIRST APV (0-127), shift it to the LAST APV (512-639)
          else if (raw_strip >= 0 && raw_strip <= 127) {
              mapped_strip = raw_strip + 512;
          }

          hit1.SetStrip(mapped_strip); // Feed the fixed strip into the reconstruction!
          // --------------------------------------------

          hit1.sety(((FD_time->at(k))-t01)*v12);
          hit1.SetCharge(max_q->at(k));
          hit1.SetTime(FD_time->at(k));
         */
        /*
         for(int k=0; k<nstrips; k++){
          hit1.SetStrip(mmStrip->at(k));
          hit1.sety(((FD_time->at(k))-t01)*v12);
          hit1.SetCharge(max_q->at(k));
          hit1.SetTime(FD_time->at(k));
          hit1.setyer(max_q->at(k));
        */
       /*
          for(int k=0; k<nstrips; k++){
          float current_charge = max_q->at(k);
          
          if (current_charge > max_allowed_charge) {
              current_charge = max_allowed_charge;
          }

          // 1. Initialize your custom map (runs only once)
          if (!mapInitialized) initializeMap();

          int raw_strip = mmStrip->at(k);
          int mapped_strip = raw_strip;


          // 3. APPLY THE FULL APV CABLE SWAP (1<->5, 2<->4)
          int final_strip = mapped_strip;


          // 4. Set the fully corrected strip
          hit1.SetStrip(final_strip); 

          // 5. Apply physical variables
          hit1.sety(((FD_time->at(k)) - t01) * v12);
          hit1.SetCharge(current_charge);
          hit1.SetTime(FD_time->at(k));
          hit1.setyer(current_charge);
    
          int chamber_idx = enumResolver[mmChamber->at(k)];

          switch(enumResolver[mmChamber->at(k)])
          {
              case 0 :  C1X.AddHit(hit1); break;
              case 1 :  C1Y.AddHit(hit1); break;
              case 2 :  C2X.AddHit(hit1); break;
              case 3 :  C2Y.AddHit(hit1); break;
              case 4 :  C3X.AddHit(hit1); break;
              case 5 :  C3Y.AddHit(hit1); break;
              case 6 :  C4X.AddHit(hit1); break;
              case 7 :  C4Y.AddHit(hit1); break;
              default : ;
          }
       } // end of loop over the strips
    */
   
    if(C1X.size()==1)  h_singlestrip_1X->Fill(C1X.Get(0).Strip());
    if(C2X.size()==1)  h_singlestrip_2X->Fill(C2X.Get(0).Strip());
    if(C3X.size()==1)  h_singlestrip_3X->Fill(C3X.Get(0).Strip());
    if(C4X.size()==1)  h_singlestrip_4X->Fill(C4X.Get(0).Strip());
    if(C1Y.size()==1)  h_singlestrip_1Y->Fill(C1Y.Get(0).Strip());
    if(C2Y.size()==1)  h_singlestrip_2Y->Fill(C2Y.Get(0).Strip());
    if(C3Y.size()==1)  h_singlestrip_3Y->Fill(C3Y.Get(0).Strip());
    if(C4Y.size()==1)  h_singlestrip_4Y->Fill(C4Y.Get(0).Strip());




        
  //      cout <<"\n\t\t\t\tNEW event  # of strips " << C1X.size() << endl;
        if(C1X.size()>1)
            for(int k=0; k<C1X.size(); k++)  h_strps_Tmm1X->Fill( C1X.Get(k).Strip());

        if(C2X.size()>1)
            for(int k=0; k<C2X.size(); k++)  h_strps_Tmm2X->Fill( C2X.Get(k).Strip());

        if(C3X.size()>1)
            for(int k=0; k<C3X.size(); k++)  h_strps_Tmm3X->Fill( C3X.Get(k).Strip());

        if(C4X.size()>1)
            for(int k=0; k<C4X.size(); k++)  h_strps_Tmm4X->Fill( C4X.Get(k).Strip());

        if(C1Y.size()>1)
            for(int k=0; k<C1Y.size(); k++)  h_strps_Tmm1Y->Fill( C1Y.Get(k).Strip());

        if(C2Y.size()>1)
            for(int k=0; k<C2Y.size(); k++)  h_strps_Tmm2Y->Fill( C2Y.Get(k).Strip());

        if(C3Y.size()>1)
            for(int k=0; k<C3Y.size(); k++)  h_strps_Tmm3Y->Fill( C3Y.Get(k).Strip());

        if(C4Y.size()>1)
            for(int k=0; k<C4Y.size(); k++)  h_strps_Tmm4Y->Fill( C4Y.Get(k).Strip());

       
    ///   CLUSTERING
        //  STEP #1 :  Sort the hits in every Layer
       C1X.SortHits();   C2X.SortHits();   C3X.SortHits();   C4X.SortHits();
       C1Y.SortHits();   C2Y.SortHits();   C3Y.SortHits();   C4Y.SortHits();
        
       //  STEP #2 : Put the hits in vector<MMHit> for the clustering method
       for (int i=0; i< C1X.size(); i++)   HT1X.push_back(C1X.Get(i));
       for (int i=0; i< C2X.size(); i++)   HT2X.push_back(C2X.Get(i));
       for (int i=0; i< C3X.size(); i++)   HT3X.push_back(C3X.Get(i));
       for (int i=0; i< C4X.size(); i++)   HT4X.push_back(C4X.Get(i));
       for (int i=0; i< C1Y.size(); i++)   HT1Y.push_back(C1Y.Get(i));
       for (int i=0; i< C2Y.size(); i++)   HT2Y.push_back(C2Y.Get(i));
       for (int i=0; i< C3Y.size(); i++)   HT3Y.push_back(C3Y.Get(i));
       for (int i=0; i< C4Y.size(); i++)   HT4Y.push_back(C4Y.Get(i));

       //  STEP #3 : The real clustering, create Cluster (in a vector<MMCluster>)
       VClustersX1 = clusterAlgo.Clustering( HT1X); VClustersY1 = clusterAlgo.Clustering( HT1Y);
       VClustersX2 = clusterAlgo.Clustering( HT2X); VClustersY2 = clusterAlgo.Clustering( HT2Y);
       VClustersX3 = clusterAlgo.Clustering( HT3X); VClustersY3 = clusterAlgo.Clustering( HT3Y);
       VClustersX4 = clusterAlgo.Clustering( HT4X); VClustersY4 = clusterAlgo.Clustering( HT4Y);

       // ==========================================
       // --- CLUSTER CHARGE SUPPRESSION ---
       // ==========================================
       float max_cluster_charge = 10000000.0; // Set your desired cluster maximum here
/*
       auto suppressClusterCharge = [&](vector<MMCluster*> &clusters) {
           for (auto* cluster : clusters) {
               float total_charge = cluster->Charge();
               
               if (total_charge > max_cluster_charge) {
                   // Calculate how much we need to scale down the cluster
                   float scale_factor = max_cluster_charge / total_charge;
                   
                   // Scale down every strip inside this cluster proportionally
                   for (int j = 0; j < (int)cluster->size(); j++) {
                       // Note: We use const_cast in case Get() returns a const reference
                       MMHit& hit = const_cast<MMHit&>(cluster->Get(j)); 
                       hit.SetCharge(hit.Charge() * scale_factor);
                   }
               }
           }
       };

       suppressClusterCharge(VClustersX1); suppressClusterCharge(VClustersY1);
       suppressClusterCharge(VClustersX2); suppressClusterCharge(VClustersY2);
       suppressClusterCharge(VClustersX3); suppressClusterCharge(VClustersY3);
       suppressClusterCharge(VClustersX4); suppressClusterCharge(VClustersY4);
       // ==========================================
*/
       auto applyClusterSizeCuts = [&](vector<MMCluster*> &v) {
           for (auto it = v.begin(); it != v.end(); ) {
               const int sz = (int)(*it)->size();
               if (sz < cl_minstrp || sz > cl_maxstrp) {
                   delete *it;
                   it = v.erase(it);
               } else {
                   ++it;
               }
           }
       };

       applyClusterSizeCuts(VClustersX1); applyClusterSizeCuts(VClustersY1);
       applyClusterSizeCuts(VClustersX2); applyClusterSizeCuts(VClustersY2);
       applyClusterSizeCuts(VClustersX3); applyClusterSizeCuts(VClustersY3);
       applyClusterSizeCuts(VClustersX4); applyClusterSizeCuts(VClustersY4);
       
// ==========================================
// --- 2. FILLING (ΜΕΣΑ ΣΤΟ LOOP) ---
// ==========================================

// Filling για τα Strips (X/Y) using the true leading cluster
#define FILL_LEADING_CLUSTER_STRIPS(VCLUS, HIST)                                          \
do {                                                                                       \
    if(!(VCLUS).empty()) {                                                                 \
        int lead = GetLeadingCluster((VCLUS));                                             \
        if(lead >= 0 && lead < (int)(VCLUS).size()) {                                      \
            MMCluster *c = (VCLUS)[lead];                                                  \
            for(int j=0; j < (int)c->size(); j++) {                                        \
                (HIST)->Fill(c->Get(j).Strip(), c->Get(j).Charge());                       \
            }                                                                              \
        }                                                                                  \
    }                                                                                      \
} while(0)

FILL_LEADING_CLUSTER_STRIPS(VClustersX1, h_leadClu_Charge_vs_Strips1X);
FILL_LEADING_CLUSTER_STRIPS(VClustersX2, h_leadClu_Charge_vs_Strips2X);
FILL_LEADING_CLUSTER_STRIPS(VClustersX3, h_leadClu_Charge_vs_Strips3X);
FILL_LEADING_CLUSTER_STRIPS(VClustersX4, h_leadClu_Charge_vs_Strips4X);

FILL_LEADING_CLUSTER_STRIPS(VClustersY1, h_leadClu_Charge_vs_Strips1Y);
FILL_LEADING_CLUSTER_STRIPS(VClustersY2, h_leadClu_Charge_vs_Strips2Y);
FILL_LEADING_CLUSTER_STRIPS(VClustersY3, h_leadClu_Charge_vs_Strips3Y);
FILL_LEADING_CLUSTER_STRIPS(VClustersY4, h_leadClu_Charge_vs_Strips4Y);

#undef FILL_LEADING_CLUSTER_STRIPS
       
       // --- Fill 1D Charge-vs-Strip using only the peak strip of the leading cluster ---
#define FILL_PEAK_STRIP_CHARGE(VCLUS, HIST, MINSIZE)                                      \
do {                                                                                       \
    if(!(VCLUS).empty()) {                                                                 \
        int lead = GetLeadingCluster((VCLUS));                                             \
        if(lead >= 0 && lead < (int)(VCLUS).size()) {                                      \
            MMCluster *c = (VCLUS)[lead];                                                  \
            if((int)c->size() >= (MINSIZE)) {                                              \
                int bestStrip = -1;                                                        \
                double bestCharge = -1.0;                                                  \
                for(int j=0; j < (int)c->size(); j++) {                                    \
                    const double q = c->Get(j).Charge();                                   \
                    if(q > bestCharge) {                                                   \
                        bestCharge = q;                                                    \
                        bestStrip = c->Get(j).Strip();                                     \
                    }                                                                      \
                }                                                                          \
                if(bestStrip >= 0) (HIST)->Fill(bestStrip, bestCharge);                   \
            }                                                                              \
        }                                                                                  \
    }                                                                                      \
} while(0)


// --- Γέμισμα για Κατανομή Landau ---
            // Χρησιμοποιούμε το Leading Cluster (που βρίσκεται στη θέση 0 μετά το swap)
            if (VClustersX1.size() > 0) h_landau_charge_X1->Fill(VClustersX1.at(0)->Charge());
            if (VClustersX2.size() > 0) h_landau_charge_X2->Fill(VClustersX2.at(0)->Charge());
            if (VClustersX3.size() > 0) h_landau_charge_X3->Fill(VClustersX3.at(0)->Charge());
            if (VClustersX4.size() > 0) h_landau_charge_X4->Fill(VClustersX4.at(0)->Charge());

            if (VClustersY1.size() > 0) h_landau_charge_Y1->Fill(VClustersY1.at(0)->Charge());
            if (VClustersY2.size() > 0) h_landau_charge_Y2->Fill(VClustersY2.at(0)->Charge());
            if (VClustersY3.size() > 0) h_landau_charge_Y3->Fill(VClustersY3.at(0)->Charge());
            if (VClustersY4.size() > 0) h_landau_charge_Y4->Fill(VClustersY4.at(0)->Charge());


// minClusterSize = 1
FILL_PEAK_STRIP_CHARGE(VClustersX1, h_charge_strip_1X1, 1);
FILL_PEAK_STRIP_CHARGE(VClustersX2, h_charge_strip_1X2, 1);
FILL_PEAK_STRIP_CHARGE(VClustersX3, h_charge_strip_1X3, 1);
FILL_PEAK_STRIP_CHARGE(VClustersX4, h_charge_strip_1X4, 1);
FILL_PEAK_STRIP_CHARGE(VClustersY1, h_charge_strip_1Y1, 1);
FILL_PEAK_STRIP_CHARGE(VClustersY2, h_charge_strip_1Y2, 1);
FILL_PEAK_STRIP_CHARGE(VClustersY3, h_charge_strip_1Y3, 1);
FILL_PEAK_STRIP_CHARGE(VClustersY4, h_charge_strip_1Y4, 1);
/*
// minClusterSize = 2
FILL_PEAK_STRIP_CHARGE(VClustersX1, h_charge_strip_2X1, 2);
FILL_PEAK_STRIP_CHARGE(VClustersX2, h_charge_strip_2X2, 2);
FILL_PEAK_STRIP_CHARGE(VClustersX3, h_charge_strip_2X3, 2);
FILL_PEAK_STRIP_CHARGE(VClustersX4, h_charge_strip_2X4, 2);
FILL_PEAK_STRIP_CHARGE(VClustersY1, h_charge_strip_2Y1, 2);
FILL_PEAK_STRIP_CHARGE(VClustersY2, h_charge_strip_2Y2, 2);
FILL_PEAK_STRIP_CHARGE(VClustersY3, h_charge_strip_2Y3, 2);
FILL_PEAK_STRIP_CHARGE(VClustersY4, h_charge_strip_2Y4, 2);

// minClusterSize = 3
FILL_PEAK_STRIP_CHARGE(VClustersX1, h_charge_strip_3X1, 3);
FILL_PEAK_STRIP_CHARGE(VClustersX2, h_charge_strip_3X2, 3);
FILL_PEAK_STRIP_CHARGE(VClustersX3, h_charge_strip_3X3, 3);
FILL_PEAK_STRIP_CHARGE(VClustersX4, h_charge_strip_3X4, 3);
FILL_PEAK_STRIP_CHARGE(VClustersY1, h_charge_strip_3Y1, 3);
FILL_PEAK_STRIP_CHARGE(VClustersY2, h_charge_strip_3Y2, 3);
FILL_PEAK_STRIP_CHARGE(VClustersY3, h_charge_strip_3Y3, 3);
FILL_PEAK_STRIP_CHARGE(VClustersY4, h_charge_strip_3Y4, 3);

#undef FILL_PEAK_STRIP_CHARGE
*/

  // --- ΓΕΜΙΣΜΑ STRIP ENTRIES ΑΝΑ ΓΩΝΙΑ (Angle-based) ---
if(angleY >= 80.0 && angleY < 88.0) {
    for(int i=0; i < (int)VClustersY1.size(); i++) {
        for(int j=0; j < (int)VClustersY1.at(i)->size(); j++) {
            h_strips_ang_80_88->Fill(VClustersY1.at(i)->Get(j).Strip());
        }
    }
}
if(angleY >= 88.0 && angleY <= 92.0) {
    for(int i=0; i < (int)VClustersY1.size(); i++) {
        for(int j=0; j < (int)VClustersY1.at(i)->size(); j++) {
            h_strips_ang_88_92->Fill(VClustersY1.at(i)->Get(j).Strip());
        }
    }
}
if(angleY > 92.0 && angleY <= 100.0) {
    for(int i=0; i < (int)VClustersY1.size(); i++) {
        for(int j=0; j < (int)VClustersY1.at(i)->size(); j++) {
            h_strips_ang_92_100->Fill(VClustersY1.at(i)->Get(j).Strip());
        }
    }
}

// Fill lead charge/strips using the true leading cluster, not index 0.
#define FILL_LEADING_CHARGE_AND_STRIPS(VCLUS, MINCOUNT, HCHARGE, HSTRIPS)                \
do {                                                                                       \
    if((int)(VCLUS).size() > (MINCOUNT)) {                                                 \
        int lead = GetLeadingCluster((VCLUS));                                             \
        if(lead >= 0 && lead < (int)(VCLUS).size()) {                                      \
            MMCluster *c = (VCLUS).at(lead);                                               \
            (HCHARGE)->Fill(c->Charge());                                                  \
            for(int j=0; j < (int)c->size(); j++) {                                        \
                (HSTRIPS)->Fill(c->Get(j).Strip());                                        \
            }                                                                              \
        }                                                                                  \
    }                                                                                      \
} while(0)
/*
// ---2 FILLING (requires >2 clusters) ---
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX1, 2, h_lead2_charge_X1, h_lead2_strips_X1);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX2, 2, h_lead2_charge_X2, h_lead2_strips_X2);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX3, 2, h_lead2_charge_X3, h_lead2_strips_X3);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX4, 2, h_lead2_charge_X4, h_lead2_strips_X4);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY1, 2, h_lead2_charge_Y1, h_lead2_strips_Y1);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY2, 2, h_lead2_charge_Y2, h_lead2_strips_Y2);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY3, 2, h_lead2_charge_Y3, h_lead2_strips_Y3);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY4, 2, h_lead2_charge_Y4, h_lead2_strips_Y4);

// ---1 FILLING (requires >1 clusters) ---
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX1, 1, h_lead1_charge_X1, h_lead1_strips_X1);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX2, 1, h_lead1_charge_X2, h_lead1_strips_X2);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX3, 1, h_lead1_charge_X3, h_lead1_strips_X3);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX4, 1, h_lead1_charge_X4, h_lead1_strips_X4);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY1, 1, h_lead1_charge_Y1, h_lead1_strips_Y1);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY2, 1, h_lead1_charge_Y2, h_lead1_strips_Y2);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY3, 1, h_lead1_charge_Y3, h_lead1_strips_Y3);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY4, 1, h_lead1_charge_Y4, h_lead1_strips_Y4);

// ---0 FILLING (requires >0 clusters) ---
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX1, 0, h_lead0_charge_X1, h_lead0_strips_X1);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX2, 0, h_lead0_charge_X2, h_lead0_strips_X2);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX3, 0, h_lead0_charge_X3, h_lead0_strips_X3);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersX4, 0, h_lead0_charge_X4, h_lead0_strips_X4);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY1, 0, h_lead0_charge_Y1, h_lead0_strips_Y1);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY2, 0, h_lead0_charge_Y2, h_lead0_strips_Y2);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY3, 0, h_lead0_charge_Y3, h_lead0_strips_Y3);
FILL_LEADING_CHARGE_AND_STRIPS(VClustersY4, 0, h_lead0_charge_Y4, h_lead0_strips_Y4);
*/
#undef FILL_LEADING_CHARGE_AND_STRIPS

// MIXED (remove if you want separation)
//h_clus_angle_X->Fill(angleX, VClustersX1.size() + VClustersX2.size() + VClustersX3.size() + VClustersX4.size());
//h_clus_angle_Y->Fill(angleY, VClustersY1.size() + VClustersY2.size() + VClustersY3.size() + VClustersY4.size());
/*
for (auto *c : VClustersX1) h_charge_angle_X->Fill(angleX, c->Charge());
for (auto *c : VClustersX2) h_charge_angle_X->Fill(angleX, c->Charge());
for (auto *c : VClustersX3) h_charge_angle_X->Fill(angleX, c->Charge());
for (auto *c : VClustersX4) h_charge_angle_X->Fill(angleX, c->Charge());

for (auto *c : VClustersY1) h_charge_angle_Y->Fill(angleY, c->Charge());
for (auto *c : VClustersY2) h_charge_angle_Y->Fill(angleY, c->Charge());
for (auto *c : VClustersY3) h_charge_angle_Y->Fill(angleY, c->Charge());
for (auto *c : VClustersY4) h_charge_angle_Y->Fill(angleY, c->Charge());
*/
// --- Filling Angle Histograms (80-120 degrees) ---


       // Cluster / event
        if(VClustersX1.size()>0) {clu_per_eventTmm1X->Fill(VClustersX1.size()); clm1X=1; clpev1x=VClustersX1.size();}else clm1X=0;
        if(VClustersX2.size()>0) {clu_per_eventTmm2X->Fill(VClustersX2.size()); clm2X=2; clpev2x=VClustersX2.size();}else clm2X=0;
        if(VClustersX3.size()>0) {clu_per_eventTmm3X->Fill(VClustersX3.size()); clm3X=4; clpev3x=VClustersX3.size();}else clm3X=0;
        if(VClustersX4.size()>0) {clu_per_eventTmm4X->Fill(VClustersX4.size()); clm4X=8; clpev4x=VClustersX4.size();}else clm4X=0;
        if(VClustersY1.size()>0) {clu_per_eventTmm1Y->Fill(VClustersY1.size()); clm1Y=1; clpev1y=VClustersY1.size();}else clm1Y=0;
        if(VClustersY2.size()>0) {clu_per_eventTmm2Y->Fill(VClustersY2.size()); clm2Y=2; clpev2y=VClustersY2.size();}else clm2Y=0;
        if(VClustersY3.size()>0) {clu_per_eventTmm3Y->Fill(VClustersY3.size()); clm3Y=4; clpev3y=VClustersY3.size();}else clm3Y=0;
        if(VClustersY4.size()>0) {clu_per_eventTmm4Y->Fill(VClustersY4.size()); clm4Y=8; clpev4y=VClustersY4.size();}else clm4Y=0;
        
      
       if(VClustersX1.size()>0)  for (int i=0; i< VClustersX1.size(); i++) {
          if(VClustersX1.at(i)->size() ==1 ) h_singlestripcluster_1X->Fill(VClustersX1.at(i)->Get(0).Strip());   }
       if(VClustersX2.size()>0)  for (int i=0; i< VClustersX2.size(); i++) {
          if(VClustersX2.at(i)->size() ==1 ) h_singlestripcluster_2X->Fill(VClustersX2.at(i)->Get(0).Strip());   }
           if(VClustersX3.size()>0)  for (int i=0; i< VClustersX3.size(); i++) {
          if(VClustersX3.at(i)->size() ==1 ) h_singlestripcluster_3X->Fill(VClustersX3.at(i)->Get(0).Strip());   }
           if(VClustersX4.size()>0)  for (int i=0; i< VClustersX4.size(); i++) {
          if(VClustersX4.at(i)->size() ==1 ) h_singlestripcluster_4X->Fill(VClustersX4.at(i)->Get(0).Strip());   }
      
      
 
        event_typeX=clm1X+clm2X+clm3X+clm4X;                 //EVENT TYPE
        if(event_typeX>0) h_Event_typeX->Fill(event_typeX);
        event_typeY=clm1Y+clm2Y+clm3Y+clm4Y;                 //EVENT TYPE
        if(event_typeY>0) h_Event_typeY->Fill(event_typeY);
      
      if (event_typeX>0 && event_typeY>0) h_Event_typeXY->Fill(event_typeX, event_typeY);

       // Cluster Charge & Cluster Length
        for (int i=0; i< VClustersX1.size(); i++) {cl_Charge_Tmm1X->Fill(VClustersX1.at(i)->Charge()); cluster_length_Tmm1X->Fill(VClustersX1.at(i)->size());}
        for (int i=0; i< VClustersX2.size(); i++) {cl_Charge_Tmm2X->Fill(VClustersX2.at(i)->Charge()); cluster_length_Tmm2X->Fill(VClustersX2.at(i)->size());}
        for (int i=0; i< VClustersX3.size(); i++) {cl_Charge_Tmm3X->Fill(VClustersX3.at(i)->Charge()); cluster_length_Tmm3X->Fill(VClustersX3.at(i)->size());}
        for (int i=0; i< VClustersX4.size(); i++) {cl_Charge_Tmm4X->Fill(VClustersX4.at(i)->Charge()); cluster_length_Tmm4X->Fill(VClustersX4.at(i)->size());}
        for (int i=0; i< VClustersY1.size(); i++) {cl_Charge_Tmm1Y->Fill(VClustersY1.at(i)->Charge()); cluster_length_Tmm1Y->Fill(VClustersY1.at(i)->size());}
        for (int i=0; i< VClustersY2.size(); i++) {cl_Charge_Tmm2Y->Fill(VClustersY2.at(i)->Charge()); cluster_length_Tmm2Y->Fill(VClustersY2.at(i)->size());}
        for (int i=0; i< VClustersY3.size(); i++) {cl_Charge_Tmm3Y->Fill(VClustersY3.at(i)->Charge()); cluster_length_Tmm3Y->Fill(VClustersY3.at(i)->size());}
        for (int i=0; i< VClustersY4.size(); i++) {cl_Charge_Tmm4Y->Fill(VClustersY4.at(i)->Charge()); cluster_length_Tmm4Y->Fill(VClustersY4.at(i)->size());}

 // find the cluster with largest charge and put it at the first position of the vector VClusters??
     
     if(VClustersX1.size()>1){ 
        leadcl=GetLeadingCluster(VClustersX1); 
        if(leadcl>0) {
            swap(VClustersX1[0], VClustersX1[leadcl]); 
        }
     }
    //     cout<< jentry << "\t"<< leadcl <<"\t"<< VClustersX1.at(leadcl)->Charge() <<"\t"<< VClustersX1.size()  <<"\t" << VClustersX1.at(0)->Charge() <<endl;
     if(VClustersX2.size()>1){ 
        leadcl=GetLeadingCluster(VClustersX2); 
        if(leadcl>0){
            swap(VClustersX2[0], VClustersX2[leadcl]); 
        }
     }
     if(VClustersX3.size()>1){ 
        leadcl=GetLeadingCluster(VClustersX3); 
        if(leadcl>0) {
            swap(VClustersX3[0], VClustersX3[leadcl]); 
        }
    }
     if(VClustersX4.size()>1){ 
        leadcl=GetLeadingCluster(VClustersX4); 
        if(leadcl>0){
            swap(VClustersX4[0], VClustersX4[leadcl]); 
        }
    }
     if(VClustersY1.size()>1){ 
        leadcl=GetLeadingCluster(VClustersY1); 
        if(leadcl>0){
            swap(VClustersY1[0], VClustersY1[leadcl]);
        }
    }
    //     cout<< jentry << "\t"<< leadcl <<"\t"<< VClustersX1.at(leadcl)->Charge() <<"\t"<< VClustersX1.size()  <<"\t" << VClustersX1.at(0)->Charge() <<endl;
     if(VClustersY2.size()>1){ 
        leadcl=GetLeadingCluster(VClustersY2); 
        if(leadcl>0) {
            swap(VClustersY2[0], VClustersY2[leadcl]); 
        }
    }
     if(VClustersY3.size()>1){ 
        leadcl=GetLeadingCluster(VClustersY3); 
        if(leadcl>0) {
            swap(VClustersY3[0], VClustersY3[leadcl]); 
        }
     }
     if(VClustersY4.size()>1){ 
        leadcl=GetLeadingCluster(VClustersY4); 
        if(leadcl>0) {
            swap(VClustersY4[0], VClustersY4[leadcl]); 
        }
    }
/*
     //  if(clpev1x*clpev1y==1)    // 1 cluster at every X Layer
     if (clpev1x>0 && clpev4x>0 && clpev2x>0 && clpev3x>0 &&
    VClustersX1.at(0)->size()>1 && VClustersX1.at(0)->size()<8 &&
    VClustersX2.at(0)->size()>1 && VClustersX2.at(0)->size()<8 &&
    VClustersX3.at(0)->size()>1 && VClustersX3.at(0)->size()<8 &&
    VClustersX4.at(0)->size()>1 && VClustersX4.at(0)->size()<8)

    if (clpev1y>0 && clpev4y>0 && clpev2y>0 && clpev3y>0 &&
    VClustersY1.at(0)->size()>1 && VClustersY1.at(0)->size()<8 &&
    VClustersY2.at(0)->size()>1 && VClustersY2.at(0)->size()<8 &&
    VClustersY3.at(0)->size()>1 && VClustersY3.at(0)->size()<8 &&
    VClustersY4.at(0)->size()>1 && VClustersY4.at(0)->size()<8)
        {    // 1 cluster at every X Layer
 
            pos1X=VClustersX1.at(0)->Pos(0.600);
            pos2X=VClustersX2.at(0)->Pos(0.600);
            pos3X=VClustersX3.at(0)->Pos(0.600);
            pos4X=VClustersX4.at(0)->Pos(0.600);
  
            // 1 cluster at every Y Layer
                pos1Y=VClustersY1.at(0)->Pos(0.600);
                pos2Y=VClustersY2.at(0)->Pos(0.600);
                pos3Y=VClustersY3.at(0)->Pos(0.600);
                pos4Y=VClustersY4.at(0)->Pos(0.600);
*/

    if (clpev1x>0 && clpev4x>0 && clpev2x>0 && clpev3x>0 &&
    VClustersX1.at(0)->size()>1 && VClustersX1.at(0)->size()<8 &&
    VClustersX2.at(0)->size()>1 && VClustersX2.at(0)->size()<8 &&
    VClustersX3.at(0)->size()>1 && VClustersX3.at(0)->size()<8 &&
    VClustersX4.at(0)->size()>1 && VClustersX4.at(0)->size()<8)

    if (clpev1y>0 && clpev4y>0 && clpev2y>0 && clpev3y>0 &&
    VClustersY1.at(0)->size()>1 && VClustersY1.at(0)->size()<8 &&
    VClustersY2.at(0)->size()>1 && VClustersY2.at(0)->size()<8 &&
    VClustersY3.at(0)->size()>1 && VClustersY3.at(0)->size()<8 &&
    VClustersY4.at(0)->size()>1 && VClustersY4.at(0)->size()<8)
        {    // 1 cluster at every X and Y Layer
            // --- OLD CENTROID CALCULATOR ---
            // --- CUSTOM CENTROID CALCULATOR (Delta-Electron Capping) ---
/*            auto getCappedPos = [](MMCluster* cluster, float pitch, float capLimit) {
                if (!cluster || cluster->size() == 0) return 0.0;
                double sumCharge = 0.0;
                double sumPosCharge = 0.0;
                for (int j = 0; j < (int)cluster->size(); j++) {
                    double q = cluster->Get(j).Charge();
                    // Lower the delta electron pull ONLY for position math
                    if (q > capLimit) q = capLimit; 
                    
                    double stripPos = cluster->Get(j).Strip() * pitch;
                    sumPosCharge += (stripPos * q);
                    sumCharge += q;
                }
                if (sumCharge == 0) return 0.0;
                return sumPosCharge / sumCharge;
            };
*/

// ---------------------------------------------------------------------
// NEW: Simple centroid for small clusters (size < 4)
//      For larger clusters, use the more robust peak-relative method.
// ---------------------------------------------------------------------
auto getCappedPos = [&](MMCluster* cluster, float pitch, float chargeCap) {
    if (!cluster || cluster->size() == 0) return 0.0;
    
    int size = cluster->size();

    // ---------------------------
    // CASE 1: Small cluster (< 4 strips)
    // ---------------------------
    if (size < 4) {
        double sum_q  = 0.0;
        double sum_qx = 0.0;
        for (int j = 0; j < size; ++j) {
            double q = cluster->Get(j).Charge();
            // Optional mild cap to avoid delta-electron bias
            if (q > chargeCap) q = chargeCap; 
            double x = (cluster->Get(j).Strip() + 0.5) * pitch;
            sum_q  += q;
            sum_qx += q * x;
        }
        if (sum_q == 0.0) return (cluster->Get(0).Strip() + 0.5) * pitch;
        return sum_qx / sum_q;
    }

    // ---------------------------
    // CASE 2: Large cluster (≥ 4 strips) – original “peak-relative” method
    // ---------------------------
    // 1. Find the peak strip
    double max_q = -1.0;
    int peak_idx = -1;
    for (int j = 0; j < size; ++j) {
        double q = cluster->Get(j).Charge();
        if (q > max_q) {
            max_q = q;
            peak_idx = j;
        }
    }

    // 2. Define a tail threshold (ignore strips with < 5% of peak charge)
    double tail_threshold = max_q * 0.05;
    double sum_q  = 0.0;
    double sum_qx = 0.0;

    for (int j = 0; j < size; ++j) {
        double q = cluster->Get(j).Charge();
        if (q < tail_threshold) continue;          // reject weak tails
        if (q > chargeCap) q = chargeCap;          // cap delta spikes
        double x = (cluster->Get(j).Strip() + 0.5) * pitch;
        sum_q  += q;
        sum_qx += q * x;
    }

    if (sum_q == 0.0) return (cluster->Get(peak_idx).Strip() + 0.5) * pitch;
    return sum_qx / sum_q;
};

/*
            // =========================================================================
            // --- UPGRADED "PEAK-RELATIVE" CENTROID CALCULATOR ---
            // =========================================================================
            auto getCappedPos = [](MMCluster* cluster, float pitch, float capLimit) {
                if (!cluster || cluster->size() == 0) return 0.0;
                
                int size = cluster->size();
                
                // 1. If it's a 1-strip cluster, return its exact physical center
                if (size == 1) {
                    return (cluster->Get(0).Strip() + 0.5) * pitch;
                }

                // 2. Find the peak strip (The absolute core of the avalanche)
                double max_q = -1.0;
                int peak_idx = -1;
                for (int j = 0; j < size; j++) {
                    double q = cluster->Get(j).Charge();
                    if (q > max_q) {
                        max_q = q;
                        peak_idx = j;
                    }
                }

                // 3. Define a dynamic threshold to ignore delta-electron tails
                // e.g., Ignore strips that have less than 5% of the core's charge
                double tail_threshold = max_q * 0.05; 

                double sum_qx = 0.0;
                double sum_q = 0.0;

                // 4. Calculate the standard Charge Centroid on the clean strips
                for (int j = 0; j < size; j++) {
                    double q = cluster->Get(j).Charge();
                    
                    // CUT 1: Remove weak asymmetric diffusion tails
                    if (q < tail_threshold) continue;
                    
                    // CUT 2: Cap massive delta-electron spikes to prevent heavy anchoring
                    if (q > capLimit) q = capLimit;

                    double x = (cluster->Get(j).Strip() + 0.5) * pitch;
                    
                    sum_qx += q * x;
                    sum_q  += q;
                }

                // Safety fallback if math fails
                if (sum_q == 0.0) return (cluster->Get(peak_idx).Strip() + 0.5) * pitch;

                return sum_qx / sum_q;
            };
*/
            float charge_cap = 2000.0; // Threshold to suppress delta ray pull
            float pitch_val = 0.600;   // Original pitch used for Pos()
            //old way of calculating positions using the standard centroid
            /*
            // Calculate positions using the delta-resistant centroid
            pos1X = getCappedPos(VClustersX1.at(0), pitch_val, charge_cap);
            pos2X = getCappedPos(VClustersX2.at(0), pitch_val, charge_cap);
            pos3X = getCappedPos(VClustersX3.at(0), pitch_val, charge_cap);
            pos4X = getCappedPos(VClustersX4.at(0), pitch_val, charge_cap);
  
            pos1Y = getCappedPos(VClustersY1.at(0), pitch_val, charge_cap);
            pos2Y = getCappedPos(VClustersY2.at(0), pitch_val, charge_cap);
            pos3Y = getCappedPos(VClustersY3.at(0), pitch_val, charge_cap);
            pos4Y = getCappedPos(VClustersY4.at(0), pitch_val, charge_cap);
            */
        


            /*
            // ==========================================
            // --- 1. TRACK SEEDING (Anchor on Layers 1 & 4) ---
            // ==========================================
            pos1X = getCappedPos(VClustersX1.at(0), pitch_val, charge_cap);
            pos4X = getCappedPos(VClustersX4.at(0), pitch_val, charge_cap);
            pos1Y = getCappedPos(VClustersY1.at(0), pitch_val, charge_cap);
            pos4Y = getCappedPos(VClustersY4.at(0), pitch_val, charge_cap);

            // --- 2. Define the theoretical "road" (Straight Line) ---
            float seed_slopeX = (pos4X - pos1X) / (pos4Z - pos1Z);
            float seed_intrcX = pos1X - (seed_slopeX * pos1Z);
            
            float seed_slopeY = (pos4Y - pos1Y) / (pos4Z - pos1Z);
            float seed_intrcY = pos1Y - (seed_slopeY * pos1Z);

            float exp_pos2X = (seed_slopeX * pos2Z) + seed_intrcX;
            float exp_pos3X = (seed_slopeX * pos3Z) + seed_intrcX;
            float exp_pos2Y = (seed_slopeY * pos2Z) + seed_intrcY;
            float exp_pos3Y = (seed_slopeY * pos3Z) + seed_intrcY;

            // --- 3. Search for the best clusters in the inner layers ---
            auto findBestCluster = [&](vector<MMCluster*>& clusters, float expected_pos) {
                float best_pos = 0;
                float min_dist = 9999.0;
                for(int j = 0; j < (int)clusters.size(); j++) {
                    float test_pos = getCappedPos(clusters.at(j), pitch_val, charge_cap);
                    float dist = fabs(test_pos - expected_pos);
                    if(dist < min_dist) {
                        min_dist = dist;
                        best_pos = test_pos;
                    }
                }
                return best_pos; // Returns the position of the cluster CLOSEST to the road!
            };

            pos2X = findBestCluster(VClustersX2, exp_pos2X);
            pos3X = findBestCluster(VClustersX3, exp_pos3X);
            pos2Y = findBestCluster(VClustersY2, exp_pos2Y);
            pos3Y = findBestCluster(VClustersY3, exp_pos3Y);
            */
        /*    
            // --- 4. PRE-FIT GEOMETRY CUT (The Garbage Collector) ---
            // If the best matching cluster is STILL more than 15mm away from the line, 
            // this is a fake track made of noise. Skip it entirely to save CPU time!
           bool valid_track = true;
            float max_tolerance = 3.0; 
            if (fabs(pos2X - exp_pos2X) > max_tolerance || fabs(pos3X - exp_pos3X) > max_tolerance ||
                fabs(pos2Y - exp_pos2Y) > max_tolerance || fabs(pos3Y - exp_pos3Y) > max_tolerance) {
                valid_track = false; 
            }
            */

            // Wrap EVERYTHING below this (down to line 1390) inside this IF statement!
             
            zpos[0]=pos1Z; zpos[1]=pos2Z; zpos[2]=pos3Z; zpos[3]=pos4Z;

            xpos[0]=pos1X+corr_pos1X;
            xpos[1]=pos2X+corr_pos2X;
            xpos[2]=pos3X+corr_pos3X;
            xpos[3]=pos4X+corr_pos4X;

            ypos[0]=pos1Y+corr_pos1Y;
            ypos[1]=pos2Y+corr_pos2Y;
            ypos[2]=pos3Y+corr_pos3Y;
            ypos[3]=pos4Y+corr_pos4Y;
/* graphx->Delete
            zpos[0]=pos1Z;zpos[1]=pos2Z;zpos[2]=pos3Z;zpos[3]=pos4Z;

           xpos[0]=pos1X+corr_pos1X;
            xpos[1]=pos2X+corr_pos2X;
            xpos[2]=pos3X+corr_pos3X;
            xpos[3]=pos4X+corr_pos4X;

            ypos[0]=pos1Y+corr_pos1Y;
            ypos[1]=pos2Y+corr_pos2Y;
            ypos[2]=pos3Y+corr_pos3Y;
            ypos[3]=pos4Y+corr_pos4Y;
*/
/*
//1st attempt at uTPC position correction (linear fit of strip vs drift)
            auto computeUtpcPosition = [&](MMCluster* cluster, double centroidPosition) {
                if (cluster == 0 || cluster->size() < 2) {
                    return centroidPosition;
                }

                const int hitCount = (int)cluster->size();
                double sumX = 0.0;
                double sumY = 0.0;
                double sumXX = 0.0;
                double sumXY = 0.0;

                for (int hitIndex = 0; hitIndex < hitCount; ++hitIndex) {
                    const double stripPositionMm = cluster->Get(hitIndex).Strip() * pitch1;
                    const double driftPosition = cluster->Get(hitIndex).gety();
                    sumX += stripPositionMm;
                    sumY += driftPosition;
                    sumXX += stripPositionMm * stripPositionMm;
                    sumXY += stripPositionMm * driftPosition;
                }

                const double denom = hitCount * sumXX - sumX * sumX;
                if (fabs(denom) < 1e-9) {
                    return centroidPosition;
                }
                // Fit X as a function of Y: X = m*Y + c
                const double m_utpc = (hitCount * sumXY - sumX * sumY) / denom;
                const double c_utpc = (sumX - m_utpc * sumY) / hitCount;

                // Evaluate X at the middle of the 5.0 mm drift gap (Y = 2.5 mm)
                return m_utpc * 2.5 + c_utpc;
                };
*/
//2nd attempt at uTPC position correction (linear fit of strip vs drift, but with the correct math this time)
/* 
//den pairnw mesaio strip alla to akriano poy einai lathos
auto computeUtpcPosition = [&](MMCluster* cluster, double centroidPosition) {
                if (cluster == 0 || cluster->size() < 2) {
                    return centroidPosition;
                }

                const int hitCount = (int)cluster->size();
                double sumX = 0.0, sumY = 0.0, sumYY = 0.0, sumXY = 0.0; // FIXED: Using sumYY

                for (int hitIndex = 0; hitIndex < hitCount; ++hitIndex) {
                    const double X = cluster->Get(hitIndex).Strip() * pitch1;
                    const double Y = cluster->Get(hitIndex).gety();
                    sumX += X;
                    sumY += Y;
                    sumYY += Y * Y;   // FIXED: Variance of Y
                    sumXY += X * Y;
                }

                const double denom = hitCount * sumYY - sumY * sumY; // FIXED: Denominator is Y
                if (fabs(denom) < 1e-9) {
                    return centroidPosition;
                }
                
                // Fit X as a function of Y: X = m*Y + c
                const double m_utpc = (hitCount * sumXY - sumX * sumY) / denom;
                const double c_utpc = (sumX - m_utpc * sumY) / hitCount;

                // Evaluate X at the middle of the 5.0 mm drift gap (Y = 2.5 mm)
                double utpcPos = m_utpc * 2.5 + c_utpc;

                // --- THE SAFETY NET ---
                // If the math ever tries to throw the position more than 1.5mm away, reject it!
                if (fabs(utpcPos - centroidPosition) > 1.5) {
                    return centroidPosition;
                }

                return utpcPos;
            };
*/
/*
//3rd attempt at uTPC position correction (core-based fit, more robust against noise and delta electrons)
//but is making angle closer to 75-105 instead of 70-110, so maybe the safety net is too tight? 

auto computeUtpcPosition = [&](MMCluster* cluster) {
    if (cluster == 0) return 0.0;
    if (cluster->size() < 2) {
        return (cluster->Get(0).Strip() + 0.5) * pitch1;
    }

    const int hitCount = (int)cluster->size();
    double maxCharge = -1.0;
    double coreX = 0.0;
    double coreY = 0.0;
    
    // 1. Find Core Strip
    for (int i = 0; i < hitCount; ++i) {
        if (cluster->Get(i).Charge() > maxCharge) {
            maxCharge = cluster->Get(i).Charge();
            coreX = (cluster->Get(i).Strip() + 0.5) * pitch1; 
            coreY = cluster->Get(i).gety();
        }
    }

    // 2. THE SQUISHER BUG (Constrained Fit)
    double sumYY = 0.0;
    double sumXY = 0.0;
    for (int i = 0; i < hitCount; ++i) {
        // By shifting everything relative to coreX, it forces the line through the center
        const double shiftedX = ((cluster->Get(i).Strip() + 0.5) * pitch1) - coreX;
        const double shiftedY = cluster->Get(i).gety() - coreY;
        const double weight = cluster->Get(i).Charge(); 
        
        sumYY += weight * shiftedY * shiftedY;
        sumXY += weight * shiftedX * shiftedY;
    }

    if (fabs(sumYY) < 1e-9) return coreX;
    const double m_utpc = sumXY / sumYY;
    
    // Safety check that also censors angles!
    double utpcAngle = 180.0 * atan2(1.0, m_utpc) / TMath::Pi();
    if (utpcAngle < 0.0) utpcAngle += 180.0;
    if (utpcAngle < 70.0 || utpcAngle > 110.0) return coreX;

    return coreX + m_utpc * (2.5 - coreY);
};
     */

        // =========================================================================
            // --- UPGRADED CHARGE-WEIGHTED uTPC (DYNAMIC MIDPOINT EVALUATION) ---
            // =========================================================================
            auto computeUtpcPosition = [&](MMCluster* cluster) {
                if (cluster == 0) return 0.0;
                
                const int hitCount = (int)cluster->size();
                
                // uTPC requires at least 2 points to draw a line. Fallback to centroid!
                if (hitCount < 2) {
                    return getCappedPos(cluster, pitch_val, charge_cap);
                }

                double Sw = 0.0;   // Sum of weights
                double Sz = 0.0;   // Sum of Z (drift distance, stored in gety())
                double Sx = 0.0;   // Sum of Strip Position (X or Y depending on the layer)
                double Szz = 0.0;  // Sum of Z^2
                double Szx = 0.0;  // Sum of Z*Strip

                // Variables to track the Z boundaries of the electron cloud
                double minZ = 9999.0;
                double maxZ = -9999.0;

                for (int i = 0; i < hitCount; ++i) {
                    double strip_pos = (cluster->Get(i).Strip() + 0.5) * pitch_val;
                    double z_drift = cluster->Get(i).gety(); 
                    double w_charge = cluster->Get(i).Charge();

                    // Apply the delta-electron cap to the weight!
                    if (w_charge > charge_cap) w_charge = charge_cap;

                    // Track the physical top and bottom of the avalanche
                    if (z_drift < minZ) minZ = z_drift;
                    if (z_drift > maxZ) maxZ = z_drift;

                    Sw  += w_charge;
                    Sz  += w_charge * z_drift;
                    Sx  += w_charge * strip_pos;
                    Szz += w_charge * z_drift * z_drift;
                    Szx += w_charge * z_drift * strip_pos;
                }

                double denom = (Sw * Szz) - (Sz * Sz);
                
                // Safety: If denom is ~0, the track is perfectly vertical. Fallback to centroid!
                if (fabs(denom) < 1e-6) {
                    return getCappedPos(cluster, pitch_val, charge_cap); 
                }

                // Fit Strip(Z) as a function of Z (Drift): Strip(Z) = m*Z + c
                double m_utpc = ((Sw * Szx) - (Sz * Sx)) / denom;
                double c_utpc = (Sx - m_utpc * Sz) / Sw;

                // --- THE NEW DYNAMIC MIDPOINT LOGIC ---
                // 1. Calculate Zutpc (The true physical center of the drift cloud)
                double z_mid = (maxZ + minZ) / 2.0;

                // 2. Calculate Xutpc or Yutpc (The strip coordinate exactly at Zutpc)
                double utpcPos = (m_utpc * z_mid) + c_utpc;

                // --- THE EXTRAPOLATION SAFETY NET ---
                double centroidPos = getCappedPos(cluster, pitch_val, charge_cap);
                double clusterWidth = hitCount * pitch_val;
                
                // If the timing data is highly distorted by noise, the line might extrapolate 
                // the hit outside the physical bounds of the cluster. Trust the centroid instead!
                if (fabs(utpcPos - centroidPos) > clusterWidth) {
                    return centroidPos; 
                }

                return utpcPos;
            };

     /*
// ==========================================
// --- CORRECTED CHARGE-WEIGHTED uTPC ---
// ==========================================
auto computeUtpcPosition = [&](MMCluster* cluster) {
    if (cluster == 0) return 0.0;
    //if (angleX < 80 || angleX > 100 || angleY < 80 || angleY > 100 )
    {
    const int hitCount = (int)cluster->size();
    
    if (hitCount < 2) {
        return (cluster->Get(0).Strip() + 0.5) * pitch1;
    }

    double sumX = 0.0, sumY = 0.0, sumYY = 0.0, sumXY = 0.0, sumW = 0.0;
    double maxCharge = -1.0;
    double fallbackX = 0.0;
    
    // Track the physical boundaries of the electron cloud
    double minY = 9999.0;
    double maxY = -9999.0;

    for (int i = 0; i < hitCount; ++i) {
        const double X = (cluster->Get(i).Strip() + 0.5) * pitch1;
        const double Y = cluster->Get(i).gety();
        const double W = cluster->Get(i).Charge(); 

        if (W > maxCharge) {
            maxCharge = W;
            fallbackX = X;
        }

        if (Y < minY) minY = Y;
        if (Y > maxY) maxY = Y;

        sumX  += W * X;
        sumY  += W * Y;
        sumYY += W * Y * Y;
        sumXY += W * X * Y;
        sumW  += W;
    }

    const double denom = sumW * sumYY - sumY * sumY;
    
    if (fabs(denom) < 1e-9) {
        return fallbackX; 
    }

    const double m_utpc = (sumW * sumXY - sumX * sumY) / denom;
    const double c_utpc = (sumX - m_utpc * sumY) / sumW;

    double utpcAngle = 180.0 * atan2(1.0, m_utpc) / TMath::Pi();
    if (utpcAngle < 0.0) utpcAngle += 180.0;

    if (utpcAngle < 65.0 || utpcAngle > 115.0) {
        return fallbackX; 
    }

    // --- THE FIX: DYNAMIC MIDPOINT ---
    // Instead of forcing 2.5, find the physical middle of the hits.
    // This prevents "Extrapolation Blowout" on tracks that clip the gap!
    double true_mid_Y = (maxY + minY) / 2.0;
    
    return (m_utpc * true_mid_Y) + c_utpc; 
};};*/
/*
// ==========================================
// --- OLS CORRECTED uTPC (Y vs X) ---
// ==========================================
auto computeUtpcPosition = [&](MMCluster* cluster) {
    if (cluster == 0) return 0.0;
    
    const int hitCount = (int)cluster->size();
    double fallbackX = (cluster->Get(0).Strip() + 0.5) * pitch1; 

    if (hitCount < 2) {
        return fallbackX;
    }

    double sumX = 0.0, sumY = 0.0, sumXX = 0.0, sumXY = 0.0;
    
    // Track the physical boundaries of the electron cloud
    double minY = 9999.0;
    double maxY = -9999.0;

    for (int i = 0; i < hitCount; ++i) {
        const double X = (cluster->Get(i).Strip() + 0.5) * pitch1;
        const double Y = cluster->Get(i).gety();

        if (Y < minY) minY = Y;
        if (Y > maxY) maxY = Y;

        // UNWEIGHTED SUMS - Restores the full geometric span
        sumX  += X;
        sumY  += Y;
        sumXX += X * X;
        sumXY += X * Y;
    }

    // FITTING Y vs X (Resolves Attenuation Bias)
    const double denom = hitCount * sumXX - sumX * sumX;
    
    // If denom is 0, the track hit exactly one X coordinate (90 degrees vertical)
    if (fabs(denom) < 1e-9) {
        return sumX / hitCount; 
    }

    // Fit Y = M * X + C
    const double M = (hitCount * sumXY - sumX * sumY) / denom;
    const double C = (sumY - M * sumX) / hitCount;

    // Prevent division by zero if track is perfectly horizontal (impossible for cosmics)
    if (fabs(M) < 1e-9) {
        return fallbackX;
    }

    // M is now dY/dX = tan(theta)
    double utpcAngle = 180.0 * atan2(M, 1.0) / TMath::Pi();
    if (utpcAngle < 0.0) utpcAngle += 180.0;

    // The track is now free to expand to its natural limits!
    if (utpcAngle < 65.0 || utpcAngle > 115.0) {
        return fallbackX; 
    }

    // DYNAMIC MIDPOINT EVALUATION
    // We solve X = (Y - C) / M to find the true center of the avalanche
    double true_mid_Y = (maxY + minY) / 2.0;
    double reconstructed_X = (true_mid_Y - C) / M;
    
    return reconstructed_X; 
};*/
/*        // Preliminary centroid-only fit to decide the reconstruction mode.
            TGraph* prelimGraphX = new TGraph(4, zpos, xpos);
            trackfit->SetParameters(0.0, 0.0);
            prelimGraphX->Fit("trackfit","q");
            const double prelimIntrcX = trackfit->GetParameter(0); // <--- ADD THIS
            const double prelimCoefX = trackfit->GetParameter(1);
            double prelimAngleX = 180.0 * atan2(1.0, prelimCoefX) / TMath::Pi();
            if (prelimAngleX < 0) prelimAngleX = 180 + prelimAngleX;
            delete prelimGraphX;

            TGraph* prelimGraphY = new TGraph(4, zpos, ypos);
            trackfit->SetParameters(0.0, 0.0);
            prelimGraphY->Fit("trackfit","q");
            const double prelimIntrcY = trackfit->GetParameter(0); // <--- ADD THIS
            const double prelimCoefY = trackfit->GetParameter(1);
            double prelimAngleY = 180.0 * atan2(1.0, prelimCoefY) / TMath::Pi();
            if (prelimAngleY < 0) prelimAngleY = 180 + prelimAngleY;
            
            delete prelimGraphY;

            const bool useCentroidX = (prelimAngleX >= 85.0 && prelimAngleX <= 95.0); 
            const bool useCentroidY = (prelimAngleY >= 85.0 && prelimAngleY <= 95.0);
*/
            /*
            const double recoPos1X = useCentroidX ? pos1X : computeUtpcPosition(VClustersX1.at(0), pos1X);
            const double recoPos2X = useCentroidX ? pos2X : computeUtpcPosition(VClustersX2.at(0), pos2X);
            const double recoPos3X = useCentroidX ? pos3X : computeUtpcPosition(VClustersX3.at(0), pos3X);
            const double recoPos4X = useCentroidX ? pos4X : computeUtpcPosition(VClustersX4.at(0), pos4X);

            const double recoPos1Y = useCentroidY ? pos1Y : computeUtpcPosition(VClustersY1.at(0), pos1Y);
            const double recoPos2Y = useCentroidY ? pos2Y : computeUtpcPosition(VClustersY2.at(0), pos2Y);
            const double recoPos3Y = useCentroidY ? pos3Y : computeUtpcPosition(VClustersY3.at(0), pos3Y);
            const double recoPos4Y = useCentroidY ? pos4Y : computeUtpcPosition(VClustersY4.at(0), pos4Y);
*/
// ==========================================
            // --- STRICT uTPC RECONSTRUCTION METHOD ---
            // ==========================================
            
            // pairnw to akrives shmeio oxi ton meso oro

/* 
// ==========================================
// --- DYNAMIC HYBRID RECONSTRUCTION ---
// ==========================================
// Use Centroid for narrow clusters (1-2 strips), uTPC for wide ones (3+ strips)

auto getHybridPos = [&](MMCluster* cluster, double centroidPos) {
    if (cluster == 0) return centroidPos; // Fallback
    
    if (cluster->size() <= 2) {
        return centroidPos; // Centroid is highly accurate for small clouds
    } else {
        return computeUtpcPosition(cluster); // uTPC resolves the angle
    }
};

// Calculate hybrid positions for X
const double recoPos1X = getHybridPos(VClustersX1.at(0), pos1X);
const double recoPos2X = getHybridPos(VClustersX2.at(0), pos2X);
const double recoPos3X = getHybridPos(VClustersX3.at(0), pos3X);
const double recoPos4X = getHybridPos(VClustersX4.at(0), pos4X);

// Calculate hybrid positions for Y
const double recoPos1Y = getHybridPos(VClustersY1.at(0), pos1Y);
const double recoPos2Y = getHybridPos(VClustersY2.at(0), pos2Y);
const double recoPos3Y = getHybridPos(VClustersY3.at(0), pos3Y);
const double recoPos4Y = getHybridPos(VClustersY4.at(0), pos4Y);
*/
/*
            // ==========================================
            // --- CENTROID-MAXIMIZED HYBRID RECONSTRUCTION ---
            // ==========================================

            // X-Axis Logic
            auto getHybridPosX = [&](MMCluster* cluster, double centroidPos) {
                if (cluster == 0) return centroidPos; 
                
                int size = cluster->size();
                
                // 1. MASSIVE POWER TO CENTROID: Always use for sizes 1, 2, and 3!
                if (size <= 3) return centroidPos;
                
                // 2. Even for Size 4+, if the track is mostly vertical, the cloud 
                // is symmetrical. The Centroid is highly accurate here. Centroid wins!
                if (prelimAngleX >= 80.0 && prelimAngleX <= 100.0) {
                    return centroidPos;
                }
                
                // 3. ONLY if the cluster is huge (4+) AND the angle is steep, 
                // the Centroid is physically biased. Allow uTPC to fix it.
                return computeUtpcPosition(cluster);
            };

            // Y-Axis Logic
            auto getHybridPosY = [&](MMCluster* cluster, double centroidPos) {
                if (cluster == 0) return centroidPos; 
                
                int size = cluster->size();
                
                // 1. MASSIVE POWER TO CENTROID
                if (size <= 3) return centroidPos;
                
                // 2. Symmetric cloud check
                if (prelimAngleY >= 80.0 && prelimAngleY <= 100.0) {
                    return centroidPos;
                }
                
                // 3. uTPC only for steep, wide clusters
                return computeUtpcPosition(cluster);
            };

            // Calculate hybrid positions for X using the newly optimized clusters
const double recoPos1X = getHybridPosX(clusX1, pos1X);
const double recoPos2X = getHybridPosX(clusX2, pos2X);
const double recoPos3X = getHybridPosX(clusX3, pos3X);
const double recoPos4X = getHybridPosX(clusX4, pos4X);

// Calculate hybrid positions for Y using the newly optimized clusters
const double recoPos1Y = getHybridPosY(clusY1, pos1Y);
const double recoPos2Y = getHybridPosY(clusY2, pos2Y);
const double recoPos3Y = getHybridPosY(clusY3, pos3Y);
const double recoPos4Y = getHybridPosY(clusY4, pos4Y);
*/
            // =========================================================================
            // --- 1. COMBINATORIAL TRACK FINDING (ITERATIVE HIT SWAPPING) ---
            // =========================================================================
            
            auto optimizeHitSelection = [&](vector<MMCluster*>& v1, vector<MMCluster*>& v2, vector<MMCluster*>& v3, vector<MMCluster*>& v4) {
                std::vector<int> sel = {0, 0, 0, 0}; 
                
                auto evalFit = [&](int i1, int i2, int i3, int i4, double& max_res, int& worst_layer) {
                    double x[4] = {
                        getCappedPos(v1.at(i1), pitch_val, charge_cap),
                        getCappedPos(v2.at(i2), pitch_val, charge_cap),
                        getCappedPos(v3.at(i3), pitch_val, charge_cap),
                        getCappedPos(v4.at(i4), pitch_val, charge_cap)
                    };
                    double z[4] = {pos1Z, pos2Z, pos3Z, pos4Z};
                    
                    double sumX=0, sumZ=0, sumZZ=0, sumXZ=0;
                    for(int i=0; i<4; ++i) { sumX+=x[i]; sumZ+=z[i]; sumZZ+=z[i]*z[i]; sumXZ+=x[i]*z[i]; }
                    double denom = 4.0 * sumZZ - sumZ * sumZ;
                    
                    if(fabs(denom) < 1e-9) { max_res = 999.0; worst_layer = -1; return; }
                    
                    double m = (4.0 * sumXZ - sumX * sumZ) / denom;
                    double c = (sumX - m * sumZ) / 4.0;
                    
                    max_res = -1.0;
                    worst_layer = -1;
                    for(int i=0; i<4; ++i) {
                        double res = fabs(x[i] - (m * z[i] + c));
                        if (res > max_res) { max_res = res; worst_layer = i; }
                    }
                };

                double current_max_res;
                int current_worst_layer;
                evalFit(sel[0], sel[1], sel[2], sel[3], current_max_res, current_worst_layer);

                int iterations = 0;
                while (current_max_res > 10.0 && current_worst_layer != -1 && iterations < 4) {
                    bool improved = false;
                    int best_j = sel[current_worst_layer];
                    double best_res = current_max_res;
                    
                    int num_clusters = 0;
                    if (current_worst_layer == 0) num_clusters = v1.size();
                    else if (current_worst_layer == 1) num_clusters = v2.size();
                    else if (current_worst_layer == 2) num_clusters = v3.size();
                    else if (current_worst_layer == 3) num_clusters = v4.size();

                    for (int j = 0; j < num_clusters; ++j) {
                        if (j == sel[current_worst_layer]) continue; 
                        
                        int test_sel[4] = {sel[0], sel[1], sel[2], sel[3]};
                        test_sel[current_worst_layer] = j;
                        
                        double test_res;
                        int test_worst;
                        evalFit(test_sel[0], test_sel[1], test_sel[2], test_sel[3], test_res, test_worst);
                        
                        if (test_res < best_res) {
                            best_res = test_res;
                            best_j = j;
                            improved = true;
                        }
                    }

                    if (improved) {
                        sel[current_worst_layer] = best_j;
                        evalFit(sel[0], sel[1], sel[2], sel[3], current_max_res, current_worst_layer);
                    } else { break; }
                    iterations++;
                }
                return sel;
            };

            std::vector<int> bestX = optimizeHitSelection(VClustersX1, VClustersX2, VClustersX3, VClustersX4);
            std::vector<int> bestY = optimizeHitSelection(VClustersY1, VClustersY2, VClustersY3, VClustersY4);

            MMCluster* clusX1 = VClustersX1.at(bestX[0]);
            MMCluster* clusX2 = VClustersX2.at(bestX[1]);
            MMCluster* clusX3 = VClustersX3.at(bestX[2]);
            MMCluster* clusX4 = VClustersX4.at(bestX[3]);

            MMCluster* clusY1 = VClustersY1.at(bestY[0]);
            MMCluster* clusY2 = VClustersY2.at(bestY[1]);
            MMCluster* clusY3 = VClustersY3.at(bestY[2]);
            MMCluster* clusY4 = VClustersY4.at(bestY[3]);

            pos1X = getCappedPos(clusX1, pitch_val, charge_cap);
            pos2X = getCappedPos(clusX2, pitch_val, charge_cap);
            pos3X = getCappedPos(clusX3, pitch_val, charge_cap);
            pos4X = getCappedPos(clusX4, pitch_val, charge_cap);

            pos1Y = getCappedPos(clusY1, pitch_val, charge_cap);
            pos2Y = getCappedPos(clusY2, pitch_val, charge_cap);
            pos3Y = getCappedPos(clusY3, pitch_val, charge_cap);
            pos4Y = getCappedPos(clusY4, pitch_val, charge_cap);

            // =========================================================================
            // --- 2. PRELIMINARY FIT (TO FIND THE ANGLE) ---
            // =========================================================================
            
            // Note: We use pos1X here, not xpos, because xpos doesn't have the corr_pos shifts yet!
            double temp_z[4] = {pos1Z, pos2Z, pos3Z, pos4Z};
            double temp_x[4] = {pos1X, pos2X, pos3X, pos4X};
            double temp_y[4] = {pos1Y, pos2Y, pos3Y, pos4Y};

            TGraph* prelimGraphX = new TGraph(4, temp_z, temp_x);
            trackfit->SetParameters(0.0, 0.0);
            prelimGraphX->Fit("trackfit","q");
            double prelimAngleX = 180.0 * atan2(1.0, trackfit->GetParameter(1)) / TMath::Pi();
            if (prelimAngleX < 0) prelimAngleX = 180 + prelimAngleX;
            delete prelimGraphX;

            TGraph* prelimGraphY = new TGraph(4, temp_z, temp_y);
            trackfit->SetParameters(0.0, 0.0);
            prelimGraphY->Fit("trackfit","q");
            double prelimAngleY = 180.0 * atan2(1.0, trackfit->GetParameter(1)) / TMath::Pi();
            if (prelimAngleY < 0) prelimAngleY = 180 + prelimAngleY;
            delete prelimGraphY;

            // =========================================================================
            // --- 3. CENTROID-MAXIMIZED HYBRID RECONSTRUCTION ---
            // =========================================================================

            auto getHybridPosX = [&](MMCluster* cluster, double centroidPos) {
                if (cluster == 0) return centroidPos; 
                int size = cluster->size();
                if (size <= 3) return centroidPos;
                if (prelimAngleX >= 80.0 && prelimAngleX <= 100.0) return centroidPos;
                return computeUtpcPosition(cluster);
            };

            auto getHybridPosY = [&](MMCluster* cluster, double centroidPos) {
                if (cluster == 0) return centroidPos; 
                int size = cluster->size();
                if (size <= 3) return centroidPos;
                if (prelimAngleY >= 80.0 && prelimAngleY <= 100.0) return centroidPos;
                return computeUtpcPosition(cluster);
            };

            const double recoPos1X = getHybridPosX(clusX1, pos1X);
            const double recoPos2X = getHybridPosX(clusX2, pos2X);
            const double recoPos3X = getHybridPosX(clusX3, pos3X);
            const double recoPos4X = getHybridPosX(clusX4, pos4X);

            const double recoPos1Y = getHybridPosY(clusY1, pos1Y);
            const double recoPos2Y = getHybridPosY(clusY2, pos2Y);
            const double recoPos3Y = getHybridPosY(clusY3, pos3Y);
            const double recoPos4Y = getHybridPosY(clusY4, pos4Y);
/*
            const double recoPos1X = computeUtpcPosition(VClustersX1.at(0));
            const double recoPos2X = computeUtpcPosition(VClustersX2.at(0));
            const double recoPos3X = computeUtpcPosition(VClustersX3.at(0));
            const double recoPos4X = computeUtpcPosition(VClustersX4.at(0));

            const double recoPos1Y = computeUtpcPosition(VClustersY1.at(0));
            const double recoPos2Y = computeUtpcPosition(VClustersY2.at(0));
            const double recoPos3Y = computeUtpcPosition(VClustersY3.at(0));
            const double recoPos4Y = computeUtpcPosition(VClustersY4.at(0));
            // ==========================================
*/
            /*
            // ==========================================
            // --- DYNAMIC RECONSTRUCTION METHOD ---
            // ==========================================
            // This function decides the best math based on the physical size of the electron cloud!

            auto getRecoPos = [&](MMCluster* cluster, double centroidPos) {
                if (cluster == 0) return centroidPos;
                
                if (cluster->size() <= 2) {
                    // Cloud is narrow. Charge sharing (Centroid) is highly accurate!
                    return centroidPos;
                } else {
                    // Cloud is wide (3+ strips). Track is angled. Use time-drift (uTPC)!
                    return computeUtpcPosition(cluster, centroidPos);
                }
            };
            // Calculate the final track coordinates based on cluster size
            const double recoPos1X = getRecoPos(VClustersX1.at(0), pos1X);
            const double recoPos2X = getRecoPos(VClustersX2.at(0), pos2X);
            const double recoPos3X = getRecoPos(VClustersX3.at(0), pos3X);
            const double recoPos4X = getRecoPos(VClustersX4.at(0), pos4X);

            const double recoPos1Y = getRecoPos(VClustersY1.at(0), pos1Y);
            const double recoPos2Y = getRecoPos(VClustersY2.at(0), pos2Y);
            const double recoPos3Y = getRecoPos(VClustersY3.at(0), pos3Y);
            const double recoPos4Y = getRecoPos(VClustersY4.at(0), pos4Y);
            // ==========================================
*/

            // ==========================================
            // --- DYNAMIC RECONSTRUCTION METHOD ---
            // ==========================================
            /*
            //pairnoun ton meso oro enw egw prepei na to kanw gia shmeia :)
            // X-Axis Logic
            auto getRecoPosX = [&](MMCluster* cluster, double centroidPos) {
                if (cluster == 0) return centroidPos;
                
                // 1. Force Centroid for extremely angled tracks (Time-Clipping Protection)
                if (prelimAngleX < 70.0 || prelimAngleX > 100.0) return centroidPos;
                
                // 2. Standard Logic: Size 1-2 use Centroid, Size 3+ use uTPC
                if (cluster->size() <= 2) return centroidPos;
                
                return computeUtpcPosition(cluster, centroidPos);
            };

            // Y-Axis Logic
            auto getRecoPosY = [&](MMCluster* cluster, double centroidPos) {
                if (cluster == 0) return centroidPos;
                
                // 1. Force Centroid for extremely angled tracks (Time-Clipping Protection)
                if (prelimAngleY < 70.0 || prelimAngleY > 100.0) return centroidPos;
                
                // 2. Standard Logic: Size 1-2 use Centroid, Size 3+ use uTPC
                if (cluster->size() <= 2) return centroidPos;
                
                return computeUtpcPosition(cluster, centroidPos);
            };

            // Calculate the final track coordinates
            const double recoPos1X = getRecoPosX(VClustersX1.at(0), pos1X);
            const double recoPos2X = getRecoPosX(VClustersX2.at(0), pos2X);
            const double recoPos3X = getRecoPosX(VClustersX3.at(0), pos3X);
            const double recoPos4X = getRecoPosX(VClustersX4.at(0), pos4X);

            const double recoPos1Y = getRecoPosY(VClustersY1.at(0), pos1Y);
            const double recoPos2Y = getRecoPosY(VClustersY2.at(0), pos2Y);
            const double recoPos3Y = getRecoPosY(VClustersY3.at(0), pos3Y);
            const double recoPos4Y = getRecoPosY(VClustersY4.at(0), pos4Y);
            // ==========================================
*/
        double zpos_d[4];
        double centroidX[4] = {pos1X, pos2X, pos3X, pos4X}; // posX variables hold the pure centroid positions
        double utpcX[4] = {
            computeUtpcPosition(VClustersX1.at(0)),
            computeUtpcPosition(VClustersX2.at(0)),
            computeUtpcPosition(VClustersX3.at(0)),
            computeUtpcPosition(VClustersX4.at(0))
        };
        double centroidY[4] = {pos1Y, pos2Y, pos3Y, pos4Y}; // posY variables hold the pure centroid positions
        double utpcY[4] = {
            computeUtpcPosition(VClustersY1.at(0)),
            computeUtpcPosition(VClustersY2.at(0)),
            computeUtpcPosition(VClustersY3.at(0)),
            computeUtpcPosition(VClustersY4.at(0))
        };


            for (int i = 0; i < 4; ++i) zpos_d[i] = static_cast<double>(zpos[i]);
            TGraph *gCentroidX = new TGraph(4, zpos_d, centroidX);
            trackfit->SetParameters(0.0, 0.0);
            gCentroidX->Fit("trackfit","q");
            double aCentroid = trackfit->GetParameter(1);
            double angleCentroid = 180.0 * atan2(1.0, aCentroid) / TMath::Pi();
            if (angleCentroid < 0) angleCentroid = 180 + angleCentroid;
            h_angleX_centroid->Fill(angleCentroid);
            delete gCentroidX;

            TGraph *gUtpcX = new TGraph(4, zpos_d, utpcX);
            trackfit->SetParameters(0.0, 0.0);
            gUtpcX->Fit("trackfit","q");
            double aUtpc = trackfit->GetParameter(1);
            double angleUtpc = 180.0 * atan2(1.0, aUtpc) / TMath::Pi();
            if (angleUtpc < 0) angleUtpc = 180 + angleUtpc;
            h_angleX_utpc->Fill(angleUtpc);
            delete gUtpcX;
    
            TGraph *gCentroidY = new TGraph(4, zpos_d, centroidY);
            trackfit->SetParameters(0.0, 0.0);
            gCentroidY->Fit("trackfit","q");
            double aCentroidY = trackfit->GetParameter(1);
            double angleCentroidY = 180.0 * atan2(1.0, aCentroidY) / TMath::Pi();
            if (angleCentroidY < 0) angleCentroidY = 180 + angleCentroidY;
            h_angleY_centroid->Fill(angleCentroidY);
            delete gCentroidY;

            TGraph *gUtpcY = new TGraph(4, zpos_d, utpcY);
            trackfit->SetParameters(0.0, 0.0);
            gUtpcY->Fit("trackfit","q");
            double aUtpcY = trackfit->GetParameter(1);
            double angleUtpcY = 180.0 * atan2(1.0, aUtpcY) / TMath::Pi();
            if (angleUtpcY < 0) angleUtpcY = 180 + angleUtpcY;
            h_angleY_utpc->Fill(angleUtpcY);
            delete gUtpcY;

            xpos[0]=recoPos1X+corr_pos1X;
            xpos[1]=recoPos2X+corr_pos2X;
            xpos[2]=recoPos3X+corr_pos3X;
            xpos[3]=recoPos4X+corr_pos4X;

            ypos[0]=recoPos1Y+corr_pos1Y;
            ypos[1]=recoPos2Y+corr_pos2Y;
            ypos[2]=recoPos3Y+corr_pos3Y;
            ypos[3]=recoPos4Y+corr_pos4Y;

            h_map_Tmm1_XY ->Fill(xpos[0], ypos[0]);                       // Hit MAP  @ Tmm1
            h_map_Tmm2_XY ->Fill(xpos[1], ypos[1]);                         // Hit MAP  @ Tmm2
            h_map_Tmm3_XY ->Fill(xpos[2], ypos[2]);                         // Hit MAP  @ Tmm3
            h_map_Tmm4_XY ->Fill(xpos[3], ypos[3]);                         // Hit MAP  @ Tmm4
            h_map_XY->Fill(xpos[0], ypos[0]);
            /*
            // --- A way to assign spatial resolution ---
            // 1-strip clusters get the mathematical 72 microns.
            // 2-strip clusters are usually very precise (charge sharing).
            // 3+ strip clusters get worse resolution because the centroid smears.
            auto getDynamicRes = [&](MMCluster* cluster) {
                    if (cluster == 0) return 0.017;
                    int size = cluster->size();
                    if (size == 1) return 0.017; // Pitch / sqrt(12)
                    if (size == 2) return 0.034; // Charge sharing yields excellent resolution!
                return 0.100 + (size * 0.02); // Wider clusters get progressively worse errors
            };

            // Now, calculate the specific error for each of the 4 layers!
            float errX1 = getDynamicRes(VClustersX1.at(0));
            float errX2 = getDynamicRes(VClustersX2.at(0));
            float errX3 = getDynamicRes(VClustersX3.at(0));
            float errX4 = getDynamicRes(VClustersX4.at(0));

            //Calculate errors for the 4 Y layers
            float errY1 = getDynamicRes(VClustersY1.at(0));
            float errY2 = getDynamicRes(VClustersY2.at(0));
            float errY3 = getDynamicRes(VClustersY3.at(0));
            float errY4 = getDynamicRes(VClustersY4.at(0));

            float z_err[4] = {0.0, 0.0, 0.0, 0.0}; 
            float xy_err[4] = {errX1, errX2, errX3, errX4}; // Feed dynamic errors into the fit!
            // Define your spatial resolution as a FLOAT
            float spatial_res = 0.072; 


            // Arrays for the 3-layer residual fits (change 'double' to 'float')
            float z_err3[3] = {0.0, 0.0, 0.0};
            float xy_err3[3] = {spatial_res, spatial_res, spatial_res};
            */
/*
if (cluster > 1) {
                int Xrec = 0, qtot = 0, Sqstrip = 0;
                Xrec = getCappedPos
                
                }
*/
/**/
            // --- A smarter way to assign spatial resolution ---
            auto getDynamicRes = [&](MMCluster* cluster) {
                //if (cluster == 0) return 0.060;
                int size = cluster->size();
                if (size == 0) return 10.173; // If no cluster, assign a large error to effectively remove it from the fit
                
                //sq = 16 adc
                // NEW DIMENSION: Size 1 means electron cloud stayed inside the safe zone!
                if (size == 1) return 0.173; // If only 1 strip, it's likely a grazing track with a tiny cloud. Assign a large error to reflect the uncertainty.
                //if (size == 1) return 3; 

                if (size == 2) return 0.35; // Charge sharing yields excellent resolution!
                if (size == 3) return 0.52; // 3-strip clusters are wider, so assign a worse resolution.
                if (size == 4) return 0.7;
                if (size == 5) return 0.87;
                if (size == 6) return 1.04;
                if (size == 7) return 1.21;
                if (size == 8) return 1.39;
                if (size == 9) return 1.57;
                //if (size == 2) return 3;
                //if (size >= 3) return 3; // For wide clusters, the centroid is unreliable. Assign a large error to effectively remove it from the fit.
                if (size >= 10) return 2.0; // For very wide clusters, the centroid is extremely unreliable. Assign a very large error to effectively remove it from the fit.
                //return 0.100 + (size * 0.02);
            };

            // ==========================================
// --- EMPIRICAL / ANALYTIC ERROR ASSIGNMENT ---
// ==========================================
// Helper lambda for analytic fallback (copy of your original formula)
auto analyticSigma = [&](MMCluster* cluster) {
    if (!cluster) return 10.0;
    int size = cluster->size();
    float pitch = 0.600;
    if (size == 1) return pitch / sqrt(12.0);
    double Q_tot = 0.0;
    double x_rec = getCappedPos(cluster, pitch, charge_cap);
    double sigma_q = 15.0;   // electronic noise
    for (int j = 0; j < size; ++j) Q_tot += cluster->Get(j).Charge();
    double sum_sq_diff = 0.0;
    for (int j = 0; j < size; ++j) {
        double x_i = (cluster->Get(j).Strip() + 0.5) * pitch;
        sum_sq_diff += (x_i - x_rec) * (x_i - x_rec);
    }
    double sigma_x = (sigma_q / Q_tot) * sqrt(sum_sq_diff);
    if (sigma_x < 0.050) sigma_x = 0.050;
    if (sigma_x > 2.0) sigma_x = 2.0;
    return sigma_x;
};

// --- X LAYERS ---
float errX1, errX2, errX3, errX4;
auto getErrorX = [&](vector<MMCluster*>& v, int layer) {
    if (v.size() == 0) return 10.0f;
    int sz = v.at(0)->size();
    if (sz < 1) sz = 1; else if (sz > 10) sz = 10;
    if (useEmpiricalSigma && empSigmaX[sz] > 0) {
        return (float)empSigmaX[sz];
    } else {
        return (float)analyticSigma(v.at(0));
    }
};
errX1 = getErrorX(VClustersX1, 1);
errX2 = getErrorX(VClustersX2, 2);
errX3 = getErrorX(VClustersX3, 3);
errX4 = getErrorX(VClustersX4, 4);

// --- Y LAYERS ---
float errY1, errY2, errY3, errY4;
auto getErrorY = [&](vector<MMCluster*>& v, int layer) {
    if (v.size() == 0) return 10.0f;
    int sz = v.at(0)->size();
    if (sz < 1) sz = 1; else if (sz > 10) sz = 10;
    if (useEmpiricalSigma && empSigmaY[sz] > 0) {
        return (float)empSigmaY[sz];
    } else {
        return (float)analyticSigma(v.at(0));
    }
};
errY1 = getErrorY(VClustersY1, 1);
errY2 = getErrorY(VClustersY2, 2);
errY3 = getErrorY(VClustersY3, 3);
errY4 = getErrorY(VClustersY4, 4);

// --- These lines already exist in your code; keep them as they are ---
float z_err[4] = {0.0, 0.0, 0.0, 0.0};
float xy_err[4] = {errX1, errX2, errX3, errX4};   // for X fit
// For Y fit, you will later use {errY1, errY2, errY3, errY4}
float y_err[4] = {errY1, errY2, errY3, errY4};
            
            // For the 3-layer unbiased fit later on (update the baseline!)
            float spatial_res = 0.060; // Base resolution updated
            float z_err3[3] = {0.0, 0.0, 0.0};
            float xy_err3[3] = {spatial_res, spatial_res, spatial_res}; // Uniform error for 3-layer fits

/*
// =========================================================================
            // --- ANALYTICAL SPATIAL ERROR ASSIGNMENT (CONSTANT NOISE) ---
            // =========================================================================
            auto getDynamicRes = [&](MMCluster* cluster) {
                // 1. Missing Hit: Assign massive error to ignore it
                if (cluster == 0) return 10.0; 
                
                int size = cluster->size();
                float pitch = 0.600; // mm
                
                // 2. Special Case: Single-strip clusters
                if (size == 1) {
                    return pitch / sqrt(12.0); // Exactly 0.173 mm
                }
                
                // 3. Prepare variables
                double Q_tot = 0.0;
                double x_rec = getCappedPos(cluster, pitch, charge_cap); 
                
                // Use a flat estimate for baseline electronic noise (e.g., 15 ADC)
                double sigma_q = 15.0; 
                
                // Calculate Total Charge
                for (int j = 0; j < size; j++) {
                    Q_tot += cluster->Get(j).Charge();
                }
                
                // Calculate the sum of squared differences
                double sum_sq_diff = 0.0;
                for (int j = 0; j < size; j++) {
                    double x_i = (cluster->Get(j).Strip() + 0.5) * pitch;
                    sum_sq_diff += (x_i - x_rec) * (x_i - x_rec);
                }
                
                // 4. THE EQUATION:
                double sigma_x = (sigma_q / Q_tot) * sqrt(sum_sq_diff);
                
                // 5. Physics Limits (Safety Clamps)
                if (sigma_x < 0.050) sigma_x = 0.050; // Minimum limit ~50 microns
                if (sigma_x > 2.0) sigma_x = 2.0;     // Maximum cap
                
                return sigma_x;
            };
            // Calculate the specific analytical error for each layer
            float errX1 = getDynamicRes(VClustersX1.at(0));
            float errX2 = getDynamicRes(VClustersX2.at(0));
            float errX3 = getDynamicRes(VClustersX3.at(0));
            float errX4 = getDynamicRes(VClustersX4.at(0));

            float errY1 = getDynamicRes(VClustersY1.at(0));
            float errY2 = getDynamicRes(VClustersY2.at(0));
            float errY3 = getDynamicRes(VClustersY3.at(0));
            float errY4 = getDynamicRes(VClustersY4.at(0));

            float z_err[4] = {0.0, 0.0, 0.0, 0.0}; 
            float xy_err[4] = {errX1, errX2, errX3, errX4}; 
            float z_err3[3] = {0.0, 0.0, 0.0};
*/
            // for X
            int limit = 100000; // chi2 limit for slope histograms
            // Use TGraphErrors and add z_err and xy_err at the end
            TGraphErrors* graph0 = new TGraphErrors(4, zpos, xpos, z_err, xy_err);
            trackfit->SetParameters(0.0, 0.0);
            graph0->Fit("trackfit","q");
            b=trackfit->GetParameter(0);      a=trackfit->GetParameter(1);
                        chi2X=trackfit->GetChisquare();
                        chisqr=chi2X;  h_track_chi2->Fill(chisqr);
            coefX=a;       intrcX=b;
            angleX = 180.0 * atan2(1.0, a) / TMath::Pi();
            if (angleX<0) angleX=180+angleX;
            slopeX->Fill(angleX);
                        if (chi2X < limit) h_slopeX_chi2lt->Fill(angleX);
            // 1. Calculate the distance from 90 degrees
            double deltaAngleX = fabs(angleX - 90.0);

            // 2. Fill the new graphs
            p_avgChi2_vs_deltaAngleX->Fill(deltaAngleX, chi2X);
delete graph0;

            // for Y
            TGraphErrors* graph1 = new TGraphErrors(4, zpos, ypos, z_err, y_err);
            trackfit->SetParameters(0.0, 0.0);
            graph1->Fit("trackfit","q");
            b=trackfit->GetParameter(0);      a=trackfit->GetParameter(1);
                        chi2Y=trackfit->GetChisquare();
                        chisqr=chi2Y;  h_track_chi2->Fill(chisqr);
            coefY=a;       intrcY=b;
            angleY = 180.0 * atan2(1.0, a) / TMath::Pi();
            if (angleY<0) angleY=180+angleY;
            slopeY->Fill(angleY);
                        if (chi2Y < limit) 
                        h_slopeY_chi2lt->Fill(angleY);
                        
            slopeXY->Fill(angleX, angleY);

            auto computeAngleFromPoints = [&](const double* xvals, const double* zvals, int count) {
                double sumX = 0.0;
                double sumZ = 0.0;
                double sumZZ = 0.0;  // <-- CHANGED FROM sumXX
                double sumXZ = 0.0;
                for (int index = 0; index < count; ++index) {
                    sumX += xvals[index];
                    sumZ += zvals[index];
                    sumZZ += zvals[index] * zvals[index]; // <-- CHANGED TO USE Z
                    sumXZ += xvals[index] * zvals[index];
                }

                // Treat Z as the independent variable
                const double denom = count * sumZZ - sumZ * sumZ; 
                if (fabs(denom) < 1e-9) {
                    return 90.0;
                }

                const double localSlope = (count * sumXZ - sumX * sumZ) / denom;
                // Use same formula as angleX: angle measured from vertical for consistency
                double localAngle = 180.0 * atan2(1.0, localSlope) / TMath::Pi();
                if (localAngle < 0.0) {
                    localAngle += 180.0;
                }
                return localAngle;
            };

            auto projectTrackX = [&](double zValue) {
                if (fabs(coefX) < 1e-9) {
                    return 0.0;
                }
                return (zValue - intrcX) / coefX;
            };

            const double x432[3] = {xpos[1], xpos[2], xpos[3]};
            const double z432[3] = {zpos[1], zpos[2], zpos[3]};
            const double x321[3] = {xpos[0], xpos[1], xpos[2]};
            const double z321[3] = {zpos[0], zpos[1], zpos[2]};
            const double x43[2] = {xpos[2], xpos[3]};
            const double z43[2] = {zpos[2], zpos[3]};
            const double x32[2] = {xpos[1], xpos[2]};
            const double z32[2] = {zpos[1], zpos[2]};
            const double x21[2] = {xpos[0], xpos[1]};
            const double z21[2] = {zpos[0], zpos[1]};

            const double angle4321Value = angleX;
            const double angle432Value = computeAngleFromPoints(x432, z432, 3);
            const double angle43Value = computeAngleFromPoints(x43, z43, 2);
            const double angle321Value = computeAngleFromPoints(x321, z321, 3);
            const double angle32Value = computeAngleFromPoints(x32, z32, 2);
            const double angle21Value = computeAngleFromPoints(x21, z21, 2);

            hslope4321->Fill(angle4321Value);
            hslope432->Fill(angle432Value);
            hslope43->Fill(angle43Value);
            hslope321->Fill(angle321Value);
            hslope32->Fill(angle32Value);
            hslope21->Fill(angle21Value);

            const double projection1Value = projectTrackX(pos1Z);
            const double projection2Value = projectTrackX(pos2Z);
            const double projection3Value = projectTrackX(pos3Z);
            const double projection4Value = projectTrackX(pos4Z);

            projection1->Fill(projection1Value);
            projection2->Fill(projection2Value);
            projection3->Fill(projection3Value);
            projection4->Fill(projection4Value);

            const double deltaPosValue1 = projection1Value - xpos[0];
            const double deltaPosValue2 = projection2Value - xpos[1];
            const double deltaPosValue3 = projection3Value - xpos[2];
            const double deltaPosValue4 = projection4Value - xpos[3];

            Delta_pos1->Fill(deltaPosValue1);
            Delta_pos2->Fill(deltaPosValue2);
            Delta_pos3->Fill(deltaPosValue3);
            Delta_pos4->Fill(deltaPosValue4);

            Delta_pos_angle1->Fill(deltaPosValue1, angle4321Value);
            Delta_pos_angle2->Fill(deltaPosValue2, angle4321Value);
            Delta_pos_angle3->Fill(deltaPosValue3, angle4321Value);
            Delta_pos_angle4->Fill(deltaPosValue4, angle4321Value);

                        // Add these right after prelimAngleX and prelimAngleY are calculated (around line 828)


// Add these right after deltaAngleX and deltaAngleY are filled (around line 935)
//delete graph0;
delete graph1;

            // ==========================================
            // DELTA ANGLE CALCULATION (CLEANED)
            // ==========================================
            // Only calculate delta angles if the 4-layer track is physically valid!
            if (angleX > 60.0 && angleX < 120.0) {
                
                const double x432[3] = {xpos[1], xpos[2], xpos[3]};
                const double z432[3] = {zpos[1], zpos[2], zpos[3]};
                const double x321[3] = {xpos[0], xpos[1], xpos[2]};
                const double z321[3] = {zpos[0], zpos[1], zpos[2]};
                const double x43[2]  = {xpos[2], xpos[3]};
                const double z43[2]  = {zpos[2], zpos[3]};
                const double x32[2]  = {xpos[1], xpos[2]};
                const double z32[2]  = {zpos[1], zpos[2]};
                const double x21[2]  = {xpos[0], xpos[1]};
                const double z21[2]  = {zpos[0], zpos[1]};

                const double angle4321Value = angleX;
                const double angle432Value = computeAngleFromPoints(x432, z432, 3);
                const double angle43Value  = computeAngleFromPoints(x43, z43, 2);
                const double angle321Value = computeAngleFromPoints(x321, z321, 3);
                const double angle32Value  = computeAngleFromPoints(x32, z32, 2);
                const double angle21Value  = computeAngleFromPoints(x21, z21, 2);

                hslope4321->Fill(angle4321Value);
                hslope432->Fill(angle432Value);
                hslope43->Fill(angle43Value);
                hslope321->Fill(angle321Value);
                hslope32->Fill(angle32Value);
                hslope21->Fill(angle21Value);

                delta_angle1->Fill(angle4321Value - angle432Value);
                delta_angle2->Fill(angle4321Value - angle43Value);
                delta_angle3->Fill(angle4321Value - angle321Value);
                delta_angle4->Fill(angle4321Value - angle21Value);
                
                // You can also move your Delta_pos and projection fills in here 
                // to completely clean up your spatial residuals!
            }

            /*
            delta_angle1->Fill(angle4321Value - angle432Value);
            delta_angle2->Fill(angle4321Value - angle43Value);
            delta_angle3->Fill(angle4321Value - angle321Value);
            delta_angle4->Fill(angle4321Value - angle21Value);

            auto fillLayerDiagnostics = [&](MMCluster* cluster, TH2F* stripChargeHist, TH1F* timeHist, TH2F* slopeAngleHist, TH2F* centroidHist, TH1F* utpcSlopeHist) {
                if (cluster == 0) {
                    return;
                }

                const int hitCount = (int)cluster->size();
                double sumX = 0.0;
                double sumY = 0.0;
                double sumXX = 0.0;
                double sumXY = 0.0;

                for (int hitIndex = 0; hitIndex < hitCount; ++hitIndex) {
                    const double stripPosition = cluster->Get(hitIndex).Strip();
                    const double stripPositionMm = stripPosition * pitch1;
                    const double driftPosition = cluster->Get(hitIndex).gety();

                    stripChargeHist->Fill(stripPosition, cluster->Get(hitIndex).Charge());
                    timeHist->Fill(cluster->Get(hitIndex).Time());

                    sumX += stripPositionMm;
                    sumY += driftPosition;
                    sumXX += stripPositionMm * stripPositionMm;
                    sumXY += stripPositionMm * driftPosition;
                }

                if (hitCount < 2) {
                    return;
                }

                const double denom = hitCount * sumXX - sumX * sumX;
                if (fabs(denom) < 1e-9) {
                    return;
                }

                const double utpcSlope = (hitCount * sumXY - sumX * sumY) / denom;
                const double utpcIntercept = (sumY - utpcSlope * sumX) / hitCount;
                double utpcAngle = 180.0 * atan(utpcSlope) / TMath::Pi();
                if (utpcAngle < 0.0) {
                    utpcAngle += 180.0;
                }

                const double centroidPosition = cluster->Pos(pitch1);
                const double utpcPosition = fabs(utpcSlope) < 1e-9 ? centroidPosition : (-utpcIntercept / utpcSlope);

                slopeAngleHist->Fill(angle4321Value, utpcAngle);
                centroidHist->Fill(centroidPosition, utpcPosition);
                utpcSlopeHist->Fill(utpcAngle);
            };
*/

auto fillLayerDiagnostics = [&](MMCluster* cluster, TH2F* stripChargeHist, TH1F* timeHist, TH2F* slopeAngleHist, TH2F* centroidHist, TH1F* utpcSlopeHist) {
                if (cluster == 0) return;

                const int hitCount = (int)cluster->size();
                double sumX = 0.0, sumY = 0.0, sumYY = 0.0, sumXY = 0.0;

                for (int hitIndex = 0; hitIndex < hitCount; ++hitIndex) {
                    const double X = cluster->Get(hitIndex).Strip() * pitch1;
                    const double Y = cluster->Get(hitIndex).gety();

                    stripChargeHist->Fill(cluster->Get(hitIndex).Strip(), cluster->Get(hitIndex).Charge());
                    timeHist->Fill(cluster->Get(hitIndex).Time());

                    sumX += X;
                    sumY += Y;
                    sumYY += Y * Y;
                    sumXY += X * Y;
                }

                if (hitCount < 2) return;

                const double denom = hitCount * sumYY - sumY * sumY;
                if (fabs(denom) < 1e-9) return;

                // Fit X as a function of Y: X = m*Y + c
                const double m_utpc = (hitCount * sumXY - sumX * sumY) / denom;
                const double c_utpc = (sumX - m_utpc * sumY) / hitCount;

                // Calculate angle relative to horizontal (Matches global track fit exactly)
                double utpcAngle = 180.0 * atan2(1.0, m_utpc) / TMath::Pi();
                if (utpcAngle < 0.0) utpcAngle += 180.0;

                const double centroidPosition = cluster->Pos(pitch1);
                const double utpcPosition = m_utpc * 2.5 + c_utpc;

                slopeAngleHist->Fill(angle4321Value, utpcAngle);
                centroidHist->Fill(centroidPosition, utpcPosition);
                utpcSlopeHist->Fill(utpcAngle);
            };

            fillLayerDiagnostics(VClustersX1.at(0), charge_per_strip_tmm1, timeT1, slope_tpc_ang1, centroid_utpc1, utpc_slop1);
            fillLayerDiagnostics(VClustersX2.at(0), charge_per_strip_tmm2, timeT2, slope_tpc_ang2, centroid_utpc2, utpc_slop2);
            fillLayerDiagnostics(VClustersX3.at(0), charge_per_strip_tmm3, timeT3, slope_tpc_ang3, centroid_utpc3, utpc_slop3);
            fillLayerDiagnostics(VClustersX4.at(0), charge_per_strip_tmm4, timeT4, slope_tpc_ang4, centroid_utpc4, utpc_slop4);

            const double angle234Value = angle432Value - 90.0;
            h1_angle_vs_q_234->Fill(angle234Value, VClustersX1.at(0)->Charge());
            h2_angle_vs_q_234->Fill(angle234Value, VClustersX2.at(0)->Charge());
            h3_angle_vs_q_234->Fill(angle234Value, VClustersX3.at(0)->Charge());
            h4_angle_vs_q_234->Fill(angle234Value, VClustersX4.at(0)->Charge());

            
            // ---> FILL TH2F FOR X LAYERS <---
            if (chi2X < 10000) {
                if(VClustersX1.size() > 0) 
                    for(int j=0; j < (int)VClustersX1.at(0)->size(); j++) h_chi2_vs_strip_1X->Fill(VClustersX1.at(0)->Get(j).Strip(), chi2X);
                if(VClustersX2.size() > 0) 
                    for(int j=0; j < (int)VClustersX2.at(0)->size(); j++) h_chi2_vs_strip_2X->Fill(VClustersX2.at(0)->Get(j).Strip(), chi2X);
                if(VClustersX3.size() > 0) 
                    for(int j=0; j < (int)VClustersX3.at(0)->size(); j++) h_chi2_vs_strip_3X->Fill(VClustersX3.at(0)->Get(j).Strip(), chi2X);
                if(VClustersX4.size() > 0) 
                    for(int j=0; j < (int)VClustersX4.at(0)->size(); j++) h_chi2_vs_strip_4X->Fill(VClustersX4.at(0)->Get(j).Strip(), chi2X);
            }

            // ---> FILL TH2F FOR Y LAYERS <---
            if (chi2Y < 10000) {
                if(VClustersY1.size() > 0) 
                    for(int j=0; j < (int)VClustersY1.at(0)->size(); j++) h_chi2_vs_strip_1Y->Fill(VClustersY1.at(0)->Get(j).Strip(), chi2Y);
                if(VClustersY2.size() > 0) 
                    for(int j=0; j < (int)VClustersY2.at(0)->size(); j++) h_chi2_vs_strip_2Y->Fill(VClustersY2.at(0)->Get(j).Strip(), chi2Y);
                if(VClustersY3.size() > 0) 
                    for(int j=0; j < (int)VClustersY3.at(0)->size(); j++) h_chi2_vs_strip_3Y->Fill(VClustersY3.at(0)->Get(j).Strip(), chi2Y);
                if(VClustersY4.size() > 0) 
                    for(int j=0; j < (int)VClustersY4.at(0)->size(); j++) h_chi2_vs_strip_4Y->Fill(VClustersY4.at(0)->Get(j).Strip(), chi2Y);
            }

if (angleX >= 60.0 && angleX <= 120.0) {
    double chargeX1 = 0.0;
    double chargeX2 = 0.0;
    double chargeX3 = 0.0;
    double chargeX4 = 0.0;
    for (auto *c : VClustersX1) chargeX1 += c->Charge();
    for (auto *c : VClustersX2) chargeX2 += c->Charge();
    for (auto *c : VClustersX3) chargeX3 += c->Charge();
    for (auto *c : VClustersX4) chargeX4 += c->Charge();
    
    // nclus>0
    if (VClustersX1.size() > 0 && VClustersX2.size() > 0 && VClustersX3.size() > 0 && VClustersX4.size() > 0) {
        h_clus_angle_X_cl0->Fill(angleX, VClustersX1.size());
        h_clus_angle_X_cl0->Fill(angleX, VClustersX2.size());
        h_clus_angle_X_cl0->Fill(angleX, VClustersX3.size());
        h_clus_angle_X_cl0->Fill(angleX, VClustersX4.size());
        h_charge_angle_X_cl0->Fill(angleX, chargeX1);
        h_charge_angle_X_cl0->Fill(angleX, chargeX2);
        h_charge_angle_X_cl0->Fill(angleX, chargeX3);
        h_charge_angle_X_cl0->Fill(angleX, chargeX4);
    }
    // nclus>1
    if (VClustersX1.size() > 1 && VClustersX2.size() > 1 && VClustersX3.size() > 1 && VClustersX4.size() > 1) {
        h_clus_angle_X_cl1->Fill(angleX, VClustersX1.size());
        h_clus_angle_X_cl1->Fill(angleX, VClustersX2.size());
        h_clus_angle_X_cl1->Fill(angleX, VClustersX3.size());
        h_clus_angle_X_cl1->Fill(angleX, VClustersX4.size());
        h_charge_angle_X_cl1->Fill(angleX, chargeX1);
        h_charge_angle_X_cl1->Fill(angleX, chargeX2);
        h_charge_angle_X_cl1->Fill(angleX, chargeX3);
        h_charge_angle_X_cl1->Fill(angleX, chargeX4);
    }
    // nclus>2
    if (VClustersX1.size() > 2 && VClustersX2.size() > 2 && VClustersX3.size() > 2 && VClustersX4.size() > 2) {
        h_clus_angle_X_cl2->Fill(angleX, VClustersX1.size());
        h_clus_angle_X_cl2->Fill(angleX, VClustersX2.size());
        h_clus_angle_X_cl2->Fill(angleX, VClustersX3.size());
        h_clus_angle_X_cl2->Fill(angleX, VClustersX4.size());
        h_charge_angle_X_cl2->Fill(angleX, chargeX1);
        h_charge_angle_X_cl2->Fill(angleX, chargeX2);
        h_charge_angle_X_cl2->Fill(angleX, chargeX3);
        h_charge_angle_X_cl2->Fill(angleX, chargeX4);
    }
}

if (angleY >= 60.0 && angleY <= 120.0) {
    double chargeY1 = 0.0;
    double chargeY2 = 0.0;
    double chargeY3 = 0.0;
    double chargeY4 = 0.0;
    for (auto *c : VClustersY1) chargeY1 += c->Charge();
    for (auto *c : VClustersY2) chargeY2 += c->Charge();
    for (auto *c : VClustersY3) chargeY3 += c->Charge();
    for (auto *c : VClustersY4) chargeY4 += c->Charge();
    
    // nclus>0
    if (VClustersY1.size() > 0 && VClustersY2.size() > 0 && VClustersY3.size() > 0 && VClustersY4.size() > 0) {
        h_clus_angle_Y_cl0->Fill(angleY, VClustersY1.size());
        h_clus_angle_Y_cl0->Fill(angleY, VClustersY2.size());
        h_clus_angle_Y_cl0->Fill(angleY, VClustersY3.size());
        h_clus_angle_Y_cl0->Fill(angleY, VClustersY4.size());
        h_charge_angle_Y_cl0->Fill(angleY, chargeY1);
        h_charge_angle_Y_cl0->Fill(angleY, chargeY2);
        h_charge_angle_Y_cl0->Fill(angleY, chargeY3);
        h_charge_angle_Y_cl0->Fill(angleY, chargeY4);
    }
    // nclus>1
    if (VClustersY1.size() > 1 && VClustersY2.size() > 1 && VClustersY3.size() > 1 && VClustersY4.size() > 1) {
        h_clus_angle_Y_cl1->Fill(angleY, VClustersY1.size());
        h_clus_angle_Y_cl1->Fill(angleY, VClustersY2.size());
        h_clus_angle_Y_cl1->Fill(angleY, VClustersY3.size());
        h_clus_angle_Y_cl1->Fill(angleY, VClustersY4.size());
        h_charge_angle_Y_cl1->Fill(angleY, chargeY1);
        h_charge_angle_Y_cl1->Fill(angleY, chargeY2);
        h_charge_angle_Y_cl1->Fill(angleY, chargeY3);
        h_charge_angle_Y_cl1->Fill(angleY, chargeY4);
    }
    // nclus>2
    if (VClustersY1.size() > 2 && VClustersY2.size() > 2 && VClustersY3.size() > 2 && VClustersY4.size() > 2) {
        h_clus_angle_Y_cl2->Fill(angleY, VClustersY1.size());
        h_clus_angle_Y_cl2->Fill(angleY, VClustersY2.size());
        h_clus_angle_Y_cl2->Fill(angleY, VClustersY3.size());
        h_clus_angle_Y_cl2->Fill(angleY, VClustersY4.size());
        h_charge_angle_Y_cl2->Fill(angleY, chargeY1);
        h_charge_angle_Y_cl2->Fill(angleY, chargeY2);
        h_charge_angle_Y_cl2->Fill(angleY, chargeY3);
        h_charge_angle_Y_cl2->Fill(angleY, chargeY4);
    }
}

            int deg = 0;
            int dif = 5;
            int lim = 5000; // chi2 limit for residual histograms

/*
            // 
            if(angleX > 60.0 && angleX < 120.0 ){//&& chi2X<lim){
        
            residual1 = (coefX * zpos[0] + intrcX) - xpos[0];
            residual2 = (coefX * zpos[1] + intrcX) - xpos[1];
            residual3 = (coefX * zpos[2] + intrcX) - xpos[2];
            residual4 = (coefX * zpos[3] + intrcX) - xpos[3];
            Residual1->Fill(residual1);
            Residual2->Fill(residual2);
            Residual3->Fill(residual3);
            Residual4->Fill(residual4);
                
            if (VClustersX1.size() > 0) residualstrips1X->Fill(VClustersX1.at(0)->Get(0).Strip(), residual1);
            if (VClustersX2.size() > 0) residualstrips2X->Fill(VClustersX2.at(0)->Get(0).Strip(), residual2);
            if (VClustersX3.size() > 0) residualstrips3X->Fill(VClustersX3.at(0)->Get(0).Strip(), residual3);
            if (VClustersX4.size() > 0) residualstrips4X->Fill(VClustersX4.at(0)->Get(0).Strip(), residual4);

            }
            if(angleY >60.0 && angleY<120.0 ){//&& chi2Y<lim){
            
            //  hslope4321->Fill(angleX);
            residual1 = (coefY * zpos[0] + intrcY) - ypos[0];
            residual2 = (coefY * zpos[1] + intrcY) - ypos[1];
            residual3 = (coefY * zpos[2] + intrcY) - ypos[2];
            residual4 = (coefY * zpos[3] + intrcY) - ypos[3];                   
            Residual1Y->Fill(residual1);
            Residual2Y->Fill(residual2);
            Residual3Y->Fill(residual3);
            Residual4Y->Fill(residual4);
            }
        */
        
        // --- UNBIASED X-DIRECTION RESIDUALS ---
            if(angleX > 60.0 && angleX < 120.0){
                hslope4321->Fill(angleX);

                // Residual 1X
                float z_excl1[3] = {zpos[1], zpos[2], zpos[3]};
                float x_excl1[3] = {xpos[1], xpos[2], xpos[3]};
                float err_excl1[3] = {errX2, errX3, errX4}; 
                TGraphErrors* fit1 = new TGraphErrors(3, z_excl1, x_excl1, z_err3, err_excl1);
                if(fit1->Fit("pol1", "Q") == 0) { 
                    float res1X = fit1->GetFunction("pol1")->Eval(zpos[0]) - xpos[0];
                    Residual1->Fill(res1X);
                    // NEW: Record residual against cluster size
                    if (VClustersX1.size() > 0) h_res_vs_sizeX->Fill(VClustersX1.at(0)->size(), res1X);
                    
                    double deltaAngleX = fabs(angleX - 90.0);
                    double ndfX = trackfit->GetNDF();
                    double chi2_ndf_X = (ndfX > 0) ? chi2X / ndfX : 0;
                    p_master_graph->Fill(deltaAngleX, chi2_ndf_X, fabs(res1X));
                }
                delete fit1;

                // Residual 2X
                float z_excl2[3] = {zpos[0], zpos[2], zpos[3]};
                float x_excl2[3] = {xpos[0], xpos[2], xpos[3]};
                float err_excl2[3] = {errX1, errX3, errX4}; 
                TGraphErrors* fit2 = new TGraphErrors(3, z_excl2, x_excl2, z_err3, err_excl2);
                if(fit2->Fit("pol1", "Q") == 0) { 
                    float res2X = fit2->GetFunction("pol1")->Eval(zpos[1]) - xpos[1];
                    Residual2->Fill(res2X);
                    // NEW: Record residual against cluster size
                    if (VClustersX2.size() > 0) h_res_vs_sizeX->Fill(VClustersX2.at(0)->size(), res2X);
                }
                delete fit2;

                // Residual 3X 
                float z_excl3[3] = {zpos[0], zpos[1], zpos[3]};
                float x_excl3[3] = {xpos[0], xpos[1], xpos[3]};
                float err_excl3[3] = {errX1, errX2, errX4}; 
                TGraphErrors* fit3 = new TGraphErrors(3, z_excl3, x_excl3, z_err3, err_excl3);
                if(fit3->Fit("pol1", "Q") == 0) {
                    float res3X = fit3->GetFunction("pol1")->Eval(zpos[2]) - xpos[2];
                    Residual3->Fill(res3X);
                    // NEW: Record residual against cluster size
                    if (VClustersX3.size() > 0) h_res_vs_sizeX->Fill(VClustersX3.at(0)->size(), res3X);
                }
                delete fit3;

                // Residual 4X
                float z_excl4[3] = {zpos[0], zpos[1], zpos[2]};
                float x_excl4[3] = {xpos[0], xpos[1], xpos[2]};
                float err_excl4[3] = {errX1, errX2, errX3};
                TGraphErrors* fit4 = new TGraphErrors(3, z_excl4, x_excl4, z_err3, err_excl4);
                if(fit4->Fit("pol1", "Q") == 0) {
                    float res4X = fit4->GetFunction("pol1")->Eval(zpos[3]) - xpos[3];
                    Residual4->Fill(res4X);
                    // NEW: Record residual against cluster size
                    if (VClustersX4.size() > 0) h_res_vs_sizeX->Fill(VClustersX4.at(0)->size(), res4X);
                }
                delete fit4;
            }

            // --- UNBIASED X-DIRECTION RESIDUALS ---
            if(angleX > 60.0 && angleX < 120.0){
                hslope4321->Fill(angleX);

                // Residual 1Y
                float z_excl1[3] = {zpos[1], zpos[2], zpos[3]};
                float y_excl1[3] = {ypos[1], ypos[2], ypos[3]};
                float err_excl1[3] = {errY2, errY3, errY4}; 
                TGraphErrors* fit1 = new TGraphErrors(3, z_excl1, y_excl1, z_err3, err_excl1);
                if(fit1->Fit("pol1", "Q") == 0) { 
                    float res1Y = fit1->GetFunction("pol1")->Eval(zpos[0]) - ypos[0];
                    Residual1Y->Fill(res1Y);
                    // NEW: Record residual against cluster size
                    if (VClustersY1.size() > 0) h_res_vs_sizeY->Fill(VClustersY1.at(0)->size(), res1Y);
                    
                    double deltaAngleY = fabs(angleY - 90.0);
                    double ndfY = trackfit->GetNDF();
                    double chi2_ndf_Y = (ndfY > 0) ? chi2Y / ndfY : 0;
                    p_master_graph->Fill(deltaAngleY, chi2_ndf_Y, fabs(res1Y));
                }
                delete fit1;

                // Residual 2Y
                float z_excl2[3] = {zpos[0], zpos[2], zpos[3]};
                float y_excl2[3] = {ypos[0], ypos[2], ypos[3]};
                float err_excl2[3] = {errY1, errY3, errY4}; 
                TGraphErrors* fit2 = new TGraphErrors(3, z_excl2, y_excl2, z_err3, err_excl2);
                if(fit2->Fit("pol1", "Q") == 0) { 
                    float res2Y = fit2->GetFunction("pol1")->Eval(zpos[1]) - ypos[1];
                    Residual2Y->Fill(res2Y);
                    // NEW: Record residual against cluster size
                    if (VClustersY2.size() > 0) h_res_vs_sizeY->Fill(VClustersY2.at(0)->size(), res2Y);
                }
                delete fit2;

                // Residual 3Y 
                float z_excl3[3] = {zpos[0], zpos[1], zpos[3]};
                float y_excl3[3] = {ypos[0], ypos[1], ypos[3]};
                float err_excl3[3] = {errX1, errX2, errX4}; 
                TGraphErrors* fit3 = new TGraphErrors(3, z_excl3, y_excl3, z_err3, err_excl3);
                if(fit3->Fit("pol1", "Q") == 0) {
                    float res3Y = fit3->GetFunction("pol1")->Eval(zpos[2]) - ypos[2];
                    Residual3Y->Fill(res3Y);
                    // NEW: Record residual against cluster size
                    if (VClustersY3.size() > 0) h_res_vs_sizeY->Fill(VClustersY3.at(0)->size(), res3Y);
                }
                delete fit3;

                // Residual 4Y
                float z_excl4[3] = {zpos[0], zpos[1], zpos[2]};
                float y_excl4[3] = {ypos[0], ypos[1], ypos[2]};
                float err_excl4[3] = {errX1, errX2, errX3};
                TGraphErrors* fit4 = new TGraphErrors(3, z_excl4, y_excl4, z_err3, err_excl4);
                if(fit4->Fit("pol1", "Q") == 0) {
                    float res4Y = fit4->GetFunction("pol1")->Eval(zpos[3]) - ypos[3];
                    Residual4Y->Fill(res4Y);
                    // NEW: Record residual against cluster size
                    if (VClustersY4.size() > 0) h_res_vs_sizeY->Fill(VClustersY4.at(0)->size(), res4Y);
                }
                delete fit4;
            }
        
        /*
            // --- UNBIASED X-DIRECTION RESIDUALS ---
                       
            if(angleX > 60.0 && angleX < 120.0){// && chi2X < lim) {
                hslope4321->Fill(angleX);

                // Residual 1X
                float z_excl1[3] = {zpos[1], zpos[2], zpos[3]};
                float x_excl1[3] = {xpos[1], xpos[2], xpos[3]};
                float err_excl1[3] = {errX2, errX3, errX4}; 
                
                TGraphErrors* fit1 = new TGraphErrors(3, z_excl1, x_excl1, z_err3, err_excl1);
                if(fit1->Fit("pol1", "Q") == 0) { 
                    float res1X = fit1->GetFunction("pol1")->Eval(zpos[0]) - xpos[0];
                    Residual1->Fill(res1X);
                    
                    double deltaAngleX = fabs(angleX - 90.0);
                    double ndfX = trackfit->GetNDF();
                    double chi2_ndf_X = 0;
                    if (ndfX > 0) chi2_ndf_X = chi2X / ndfX;
                    p_master_graph->Fill(deltaAngleX, chi2_ndf_X, fabs(res1X));
                }
                delete fit1;

                // Residual 2X
                float z_excl2[3] = {zpos[0], zpos[2], zpos[3]};
                float x_excl2[3] = {xpos[0], xpos[2], xpos[3]};
                float err_excl2[3] = {errX1, errX3, errX4}; 
                
                TGraphErrors* fit2 = new TGraphErrors(3, z_excl2, x_excl2, z_err3, err_excl2);
                if(fit2->Fit("pol1", "Q") == 0) { 
                    Residual2->Fill(fit2->GetFunction("pol1")->Eval(zpos[1]) - xpos[1]);

                    for(int i=0; i<=99; i++) {
                        float test_shift = -0.5 + (i * 0.05); 
                        float test_xpos = pos2X + test_shift;    
                        float test_res = fit2->GetFunction("pol1")->Eval(zpos[1]) - test_xpos;
                        scan_hist[i]->Fill(test_res);
                    }
                }
                delete fit2;

                // Residual 3X 
                float z_excl3[3] = {zpos[0], zpos[1], zpos[3]};
                float x_excl3[3] = {xpos[0], xpos[1], xpos[3]};
                float err_excl3[3] = {errX1, errX2, errX4}; 
                
                TGraphErrors* fit3 = new TGraphErrors(3, z_excl3, x_excl3, z_err3, err_excl3);
                if(fit3->Fit("pol1", "Q") == 0) {
                    Residual3->Fill(fit3->GetFunction("pol1")->Eval(zpos[2]) - xpos[2]);
                }
                delete fit3;

                // Residual 4X
                float z_excl4[3] = {zpos[0], zpos[1], zpos[2]};
                float x_excl4[3] = {xpos[0], xpos[1], xpos[2]};
                float err_excl4[3] = {errX1, errX2, errX3};
                
                TGraphErrors* fit4 = new TGraphErrors(3, z_excl4, x_excl4, z_err3, err_excl4);
                if(fit4->Fit("pol1", "Q") == 0) {
                    Residual4->Fill(fit4->GetFunction("pol1")->Eval(zpos[3]) - xpos[3]);
                }
                delete fit4;
            }
        
        
            // --- UNBIASED Y-DIRECTION RESIDUALS ---
            if(angleY > 60.0 && angleY < 120.0 ){//&& chi2Y < lim) { //&& chi2Y < lim) {
                hslope4321->Fill(angleY);

                // Residual 1X (Excludes Layer 1, uses errors from 2, 3, 4)
                float z_excl1[3] = {zpos[1], zpos[2], zpos[3]};
                float y_excl1[3] = {ypos[1], ypos[2], ypos[3]};
                float err_excl1[3] = {errY2, errY3, errY4}; // <--- Dynamic!
                
                TGraphErrors* fit1 = new TGraphErrors(3, z_excl1, y_excl1, z_err3, err_excl1);
                if(fit1->Fit("pol1", "Q") == 0) { 
                    float res1Y = fit1->GetFunction("pol1")->Eval(zpos[0]) - ypos[0];
                    Residual1Y->Fill(res1Y);
                    
                    // --- FILL THE ULTIMATE 3D GRAPH ---
                    double deltaAngleY = fabs(angleY - 90.0);
                    double ndfY = trackfit->GetNDF();
                    double chi2_ndf_Y = 0;
                    if (ndfY > 0) chi2_ndf_Y = chi2Y / ndfY;
                    p_master_graph->Fill(deltaAngleY, chi2_ndf_Y, fabs(res1Y));
                    // ----------------------------------
                }
                delete fit1;

                // Residual 2X (Excludes Layer 2, uses errors from 1, 3, 4)
                float z_excl2[3] = {zpos[0], zpos[2], zpos[3]};
                float y_excl2[3] = {ypos[0], ypos[2], ypos[3]};
                float err_excl2[3] = {errY1, errY3, errY4}; // <--- Dynamic!
                
                TGraphErrors* fit2 = new TGraphErrors(3, z_excl2, y_excl2, z_err3, err_excl2);
                if(fit2->Fit("pol1", "Q") == 0) {
                    Residual2Y->Fill(fit2->GetFunction("pol1")->Eval(zpos[1]) - ypos[1]);
                }
                delete fit2;

                // Residual 3X (Excludes Layer 3, uses errors from 1, 2, 4)
                float z_excl3[3] = {zpos[0], zpos[1], zpos[3]};
                float y_excl3[3] = {ypos[0], ypos[1], ypos[3]};
                float err_excl3[3] = {errY1, errY2, errY4}; // <--- Dynamic!
                
                TGraphErrors* fit3 = new TGraphErrors(3, z_excl3, y_excl3, z_err3, err_excl3);
                if(fit3->Fit("pol1", "Q") == 0) {
                    Residual3Y->Fill(fit3->GetFunction("pol1")->Eval(zpos[2]) - ypos[2]);
                }
                delete fit3;

                // Residual 4X (Excludes Layer 4, uses errors from 1, 2, 3)
                float z_excl4[3] = {zpos[0], zpos[1], zpos[2]};
                float y_excl4[3] = {ypos[0], ypos[1], ypos[2]};
                float err_excl4[3] = {errY1, errY2, errY3}; // <--- Dynamic!
                
                TGraphErrors* fit4 = new TGraphErrors(3, z_excl4, y_excl4, z_err3, err_excl4);
                if(fit4->Fit("pol1", "Q") == 0) {
                    Residual4Y->Fill(fit4->GetFunction("pol1")->Eval(zpos[3]) - ypos[3]);
                }
                delete fit4;
            }*/
               
            
           /*
            // for slopeX
            TGraphErrors* graphx = new TGraphErrors(4, zpos, xpos, z_err, xy_err);
            trackfit->SetParameters(0.0, 0.0);
            graphx->Fit("trackfit","q");
            b=trackfit->GetParameter(0);      a=trackfit->GetParameter(1);
            coefX=a;       intrcX=b;
//            angleX = 180.0 * atan2(1.0, a) / TMath::Pi();
//            if (angleX<0) angleX=180+angleX;
            slopeX->Fill(angleX);
            // Fill slope histograms
            h_slope_X->Fill(a);
            h_intercept_X->Fill(b);
            h_slope_vs_intercept_X->Fill(a, b);

            // for slopeY
            TGraphErrors* graphy = new TGraphErrors(4, zpos, ypos, z_err, xy_err);
            trackfit->SetParameters(0.0, 0.0);
            graphy->Fit("trackfit","q");
            b=trackfit->GetParameter(0);      a=trackfit->GetParameter(1);
            coefY=a;       intrcY=b;
//            angleY = 180.0 * atan2(1.0, a) / TMath::Pi();
//            if (angleY<0) angleY=180+angleY;
            slopeY->Fill(angleY);
            // Fill slope histograms
            h_slope_Y->Fill(a);
            h_intercept_Y->Fill(b);
            h_slope_vs_intercept_Y->Fill(a, b);
*/


// Apply your geometry corrections
xpos[0] = recoPos1X + corr_pos1X;
xpos[1] = recoPos2X + corr_pos2X;
xpos[2] = recoPos3X + corr_pos3X;
xpos[3] = recoPos4X + corr_pos4X;

ypos[0] = recoPos1Y + corr_pos1Y;
ypos[1] = recoPos2Y + corr_pos2Y;
ypos[2] = recoPos3Y + corr_pos3Y;
ypos[3] = recoPos4Y + corr_pos4Y;

// ==========================================
// --- NOW PERFORM YOUR HYBRID FIT ---
// ==========================================

TGraphErrors* graphx = new TGraphErrors(4, zpos, xpos, z_err, xy_err);
trackfit->SetParameters(0.0, 0.0);
graphx->Fit("trackfit","q");
b = trackfit->GetParameter(0);      
a = trackfit->GetParameter(1);
coefX = a;       
intrcX = b;

// angleX is calculated globally earlier, but if you want to recalculate it based on this new fit:
angleX = 180.0 * atan2(1.0, a) / TMath::Pi();
if (angleX < 0) angleX = 180.0 + angleX;

slopeX->Fill(angleX);
h_slope_X->Fill(a);
h_intercept_X->Fill(b);
h_slope_vs_intercept_X->Fill(a, b);

delete graphx;

TGraphErrors* graphy = new TGraphErrors(4, zpos, ypos, z_err, xy_err);
trackfit->SetParameters(0.0, 0.0);
graphy->Fit("trackfit","q");
b = trackfit->GetParameter(0);      
a = trackfit->GetParameter(1);
coefY = a;       
intrcY = b;

// angleY is calculated globally earlier, but if you want to recalculate it based on this new fit:
angleY = 180.0 * atan2(1.0, a) / TMath::Pi();
if (angleY < 0) angleY = 180.0 + angleY;

slopeY->Fill(angleY);
h_slope_Y->Fill(a);
h_intercept_Y->Fill(b);
h_slope_vs_intercept_Y->Fill(a, b);

delete graphy;
            // Fill strip correlation maps using all clusters from each XY pair.
            auto fillStripPairMapAll = [&](const vector<MMCluster*>& vx, const vector<MMCluster*>& vy, TH2F* hmap, int minClusters) {
                if ((int)vx.size() <= minClusters || (int)vy.size() <= minClusters) return;
                for (auto* cx : vx) {
                    for (auto* cy : vy) {
                        for (int ix = 0; ix < (int)cx->size(); ++ix) {
                            for (int iy = 0; iy < (int)cy->size(); ++iy) {
                                hmap->Fill(cx->Get(ix).Strip(), cy->Get(iy).Strip());
                            }
                        }
                    }
                }
            };
/*
            // Require both fitted angles in the requested window for cleaner selection.
            if (angleX >= 60.0 && angleX < 80.0 && angleY >= 60.0 && angleY < 80.0) {
                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_60_80_X1Y1, 0);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_60_80_X2Y2, 0);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_60_80_X3Y3, 0);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_60_80_X4Y4, 0);

                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_60_80_X1Y1_cl1, 1);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_60_80_X2Y2_cl1, 1);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_60_80_X3Y3_cl1, 1);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_60_80_X4Y4_cl1, 1);

                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_60_80_X1Y1_cl2, 2);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_60_80_X2Y2_cl2, 2);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_60_80_X3Y3_cl2, 2);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_60_80_X4Y4_cl2, 2);
            }
*/
            if (angleX >= 85.0 && angleX <= 95.0 && angleY >= 85.0 && angleY <= 95.0) {
                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_85_95_X1Y1, 0);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_85_95_X2Y2, 0);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_85_95_X3Y3, 0);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_85_95_X4Y4, 0);

                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_85_95_X1Y1_cl1, 1);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_85_95_X2Y2_cl1, 1);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_85_95_X3Y3_cl1, 1);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_85_95_X4Y4_cl1, 1);

                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_85_95_X1Y1_cl2, 2);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_85_95_X2Y2_cl2, 2);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_85_95_X3Y3_cl2, 2);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_85_95_X4Y4_cl2, 2);
            }

            if (angleX >= 100.0 && angleX <= 120.0 && angleY >= 100.0 && angleY <= 120.0) {
                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_100_120_X1Y1, 0);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_100_120_X2Y2, 0);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_100_120_X3Y3, 0);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_100_120_X4Y4, 0);

                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_100_120_X1Y1_cl1, 1);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_100_120_X2Y2_cl1, 1);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_100_120_X3Y3_cl1, 1);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_100_120_X4Y4_cl1, 1);

                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_100_120_X1Y1_cl2, 2);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_100_120_X2Y2_cl2, 2);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_100_120_X3Y3_cl2, 2);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_100_120_X4Y4_cl2, 2);
            }
            int limitxy = 1000;
            if (chi2X < limitxy && chi2Y < limitxy) {            
                if (angleX >= 60.0 && angleX <= 120.0) {
                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_X1Y1cl0, 0);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_X2Y2cl0, 0);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_X3Y3cl0, 0);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_X4Y4cl0, 0);
            /*
                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_X1Y1cl1, 1);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_X2Y2cl1, 1);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_X3Y3cl1, 1);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_X4Y4cl1, 1);

                fillStripPairMapAll(VClustersX1, VClustersY1, h_stripXY_X1Y1cl2, 2);
                fillStripPairMapAll(VClustersX2, VClustersY2, h_stripXY_X2Y2cl2, 2);
                fillStripPairMapAll(VClustersX3, VClustersY3, h_stripXY_X3Y3cl2, 2);
                fillStripPairMapAll(VClustersX4, VClustersY4, h_stripXY_X4Y4cl2, 2);
            */}
                if(chi2X<1000) {
        if(C1X.size()>1) h_stripschi2_1X->Fill(C1X.Get(0).Strip());
        if(C2X.size()>1) h_stripschi2_2X->Fill(C2X.Get(0).Strip());
        if(C3X.size()>1) h_stripschi2_3X->Fill(C3X.Get(0).Strip());
        if(C4X.size()>1) h_stripschi2_4X->Fill(C4X.Get(0).Strip());
    }
    if(chi2Y<1000) {
        if(C1Y.size()>1) h_stripschi2_1Y->Fill(C1Y.Get(0).Strip());
        if(C2Y.size()>1) h_stripschi2_2Y->Fill(C2Y.Get(0).Strip());
        if(C3Y.size()>1) h_stripschi2_3Y->Fill(C3Y.Get(0).Strip());
        if(C4Y.size()>1) h_stripschi2_4Y->Fill(C4Y.Get(0).Strip());
    }   
        }
        // ---> NEW LAYER 1X vs 2X MAP <---
                fillStripPairMapAll(VClustersX1, VClustersX2, h_stripXY_X1X2, 0);
                fillStripPairMapAll(VClustersX3, VClustersX4, h_stripXY_X3X4, 0);

            // --- Per-layer counts ---
const int nx1 = (int)VClustersX1.size();
const int nx2 = (int)VClustersX2.size();
const int nx3 = (int)VClustersX3.size();
const int nx4 = (int)VClustersX4.size();

const int ny1 = (int)VClustersY1.size();
const int ny2 = (int)VClustersY2.size();
const int ny3 = (int)VClustersY3.size();
const int ny4 = (int)VClustersY4.size();

const int nclusX[] = {nx1, nx2, nx3, nx4};
const int nclusY[] = {ny1, ny2, ny3, ny4};

// --- Per-layer AND logic (each layer must pass threshold) ---
const bool passX_gt0 = (nx1 > 0 && nx2 > 0 && nx3 > 0 && nx4 > 0);
const bool passX_gt1 = (nx1 > 1 && nx2 > 1 && nx3 > 1 && nx4 > 1);
const bool passX_gt2 = (nx1 > 2 && nx2 > 2 && nx3 > 2 && nx4 > 2);

const bool passY_gt0 = (ny1 > 0 && ny2 > 0 && ny3 > 0 && ny4 > 0);
const bool passY_gt1 = (ny1 > 1 && ny2 > 1 && ny3 > 1 && ny4 > 1);
const bool passY_gt2 = (ny1 > 2 && ny2 > 2 && ny3 > 2 && ny4 > 2);

// --- Fill X histograms (60-120 degrees) ---
if (angleX >= 60.0 && angleX <= 120.0) {
    // Total size (sum of 4 layers)
    if (passX_gt0) {
        h_cluSize_vs_AngleX_gt0->Fill(angleX, nclusX[0]);
        if (nclusX[0] > 0) for (auto *c : VClustersX1) h_cluCharge_vs_AngleX_gt0->Fill(angleX, c->Charge());
    }
    if (passX_gt0) {
        h_cluSize_vs_AngleX_gt0->Fill(angleX, nclusX[1]);    
        if (nclusX[1] > 0) for (auto *c : VClustersX2) h_cluCharge_vs_AngleX_gt0->Fill(angleX, c->Charge());
    }
    if (passX_gt0){
        h_cluSize_vs_AngleX_gt0->Fill(angleX, nclusX[2]);
        if (nclusX[2] > 0) for (auto *c : VClustersX3) h_cluCharge_vs_AngleX_gt0->Fill(angleX, c->Charge());
    }
    if (passX_gt0) {
        h_cluSize_vs_AngleX_gt0->Fill(angleX, nclusX[3]);
        if (nclusX[3] > 0) for (auto *c : VClustersX4) h_cluCharge_vs_AngleX_gt0->Fill(angleX, c->Charge());
    }
    if (passX_gt1) {
        h_cluSize_vs_AngleX_gt1->Fill(angleX, nclusX[0]);
        if (nclusX[0] > 1) for (auto *c : VClustersX1) h_cluCharge_vs_AngleX_gt1->Fill(angleX, c->Charge());
    }
    if (passX_gt1) {
        h_cluSize_vs_AngleX_gt1->Fill(angleX, nclusX[1]);
        if (nclusX[1] > 1) for (auto *c : VClustersX2) h_cluCharge_vs_AngleX_gt1->Fill(angleX, c->Charge());
    }
    if (passX_gt1) {
        h_cluSize_vs_AngleX_gt1->Fill(angleX, nclusX[2]);
        if (nclusX[2] > 1) for (auto *c : VClustersX3) h_cluCharge_vs_AngleX_gt1->Fill(angleX, c->Charge());
    }
    if (passX_gt1) {
        h_cluSize_vs_AngleX_gt1->Fill(angleX, nclusX[3]);
        if (nclusX[3] > 1) for (auto *c : VClustersX4) h_cluCharge_vs_AngleX_gt1->Fill(angleX, c->Charge());
    }
    if (passX_gt2) {
        h_cluSize_vs_AngleX_gt2->Fill(angleX, nclusX[0]);
        if (nclusX[0] > 2) for (auto *c : VClustersX1) h_cluCharge_vs_AngleX_gt2->Fill(angleX, c->Charge());
    }
        if (passX_gt2) {
        h_cluSize_vs_AngleX_gt2->Fill(angleX, nclusX[1]);   
        if (nclusX[1] > 2) for (auto *c : VClustersX2) h_cluCharge_vs_AngleX_gt2->Fill(angleX, c->Charge());
    }
    if (passX_gt2) {
        h_cluSize_vs_AngleX_gt2->Fill(angleX, nclusX[2]);    
        if (nclusX[2] > 2) for (auto *c : VClustersX3) h_cluCharge_vs_AngleX_gt2->Fill(angleX, c->Charge());
    }
    if (passX_gt2) {
        h_cluSize_vs_AngleX_gt2->Fill(angleX, nclusX[3]);
        if (nclusX[3] > 2) for (auto *c : VClustersX4) h_cluCharge_vs_AngleX_gt2->Fill(angleX, c->Charge());
    }
/*
    // Per-layer (individual checks)
    if (nclusX[0] > 0) h_cluSize_vs_AngleX1_gt0->Fill(angleX, nclusX[0]);
    if (nclusX[1] > 0) h_cluSize_vs_AngleX2_gt0->Fill(angleX, nclusX[1]);
    if (nclusX[2] > 0) h_cluSize_vs_AngleX3_gt0->Fill(angleX, nclusX[2]);
    if (nclusX[3] > 0) h_cluSize_vs_AngleX4_gt0->Fill(angleX, nclusX[3]);

    if (nclusX[0] > 1) h_cluSize_vs_AngleX1_gt1->Fill(angleX, nclusX[0]);
    if (nclusX[1] > 1) h_cluSize_vs_AngleX2_gt1->Fill(angleX, nclusX[1]);
    if (nclusX[2] > 1) h_cluSize_vs_AngleX3_gt1->Fill(angleX, nclusX[2]);
    if (nclusX[3] > 1) h_cluSize_vs_AngleX4_gt1->Fill(angleX, nclusX[3]);

    if (nclusX[0] > 2) h_cluSize_vs_AngleX1_gt2->Fill(angleX, nclusX[0]);
    if (nclusX[1] > 2) h_cluSize_vs_AngleX2_gt2->Fill(angleX, nclusX[1]);
    if (nclusX[2] > 2) h_cluSize_vs_AngleX3_gt2->Fill(angleX, nclusX[2]);
    if (nclusX[3] > 2) h_cluSize_vs_AngleX4_gt2->Fill(angleX, nclusX[3]);
*/
/*
    // Per-layer charge fills (gt0, gt1, gt2)
    if (nclusX[0] > 0) for (auto *c : VClustersX1) h_cluCharge_vs_AngleX1_gt0->Fill(angleX, c->Charge());
    if (nclusX[1] > 0) for (auto *c : VClustersX2) h_cluCharge_vs_AngleX2_gt0->Fill(angleX, c->Charge());
    if (nclusX[2] > 0) for (auto *c : VClustersX3) h_cluCharge_vs_AngleX3_gt0->Fill(angleX, c->Charge());
    if (nclusX[3] > 0) for (auto *c : VClustersX4) h_cluCharge_vs_AngleX4_gt0->Fill(angleX, c->Charge());

    if (nclusX[0] > 1) for (auto *c : VClustersX1) h_cluCharge_vs_AngleX1_gt1->Fill(angleX, c->Charge());
    if (nclusX[1] > 1) for (auto *c : VClustersX2) h_cluCharge_vs_AngleX2_gt1->Fill(angleX, c->Charge());
    if (nclusX[2] > 1) for (auto *c : VClustersX3) h_cluCharge_vs_AngleX3_gt1->Fill(angleX, c->Charge());
    if (nclusX[3] > 1) for (auto *c : VClustersX4) h_cluCharge_vs_AngleX4_gt1->Fill(angleX, c->Charge());

    if (nclusX[0] > 2) for (auto *c : VClustersX1) h_cluCharge_vs_AngleX1_gt2->Fill(angleX, c->Charge());
    if (nclusX[1] > 2) for (auto *c : VClustersX2) h_cluCharge_vs_AngleX2_gt2->Fill(angleX, c->Charge());
    if (nclusX[2] > 2) for (auto *c : VClustersX3) h_cluCharge_vs_AngleX3_gt2->Fill(angleX, c->Charge());
    if (nclusX[3] > 2) for (auto *c : VClustersX4) h_cluCharge_vs_AngleX4_gt2->Fill(angleX, c->Charge());
*/}

// Fill helper
auto fillSizeChargeByAngle = [&](double ang, MMCluster* c,
                                 TH2F* h60, TH2F* h90, TH2F* h120) {
    const double sz = c->size();
    const double ch = c->Charge();

    if (ang >= 60.0 && ang < 85.0)       h60->Fill(sz, ch);
    else if (ang >= 85.0 && ang < 95.0)  h90->Fill(sz, ch);
    else if (ang >= 95.0 && ang <= 120.0) h120->Fill(sz, ch);
};

/*
// X clusters 60-120 degrees (separate by layer)
for (auto *c : VClustersX1) fillSizeChargeByAngle(angleX, c, h_cluSize_vs_cluCharge_60X1, h_cluSize_vs_cluCharge_90X1, h_cluSize_vs_cluCharge_120X1);
for (auto *c : VClustersX2) fillSizeChargeByAngle(angleX, c, h_cluSize_vs_cluCharge_60X2, h_cluSize_vs_cluCharge_90X2, h_cluSize_vs_cluCharge_120X2);
for (auto *c : VClustersX3) fillSizeChargeByAngle(angleX, c, h_cluSize_vs_cluCharge_60X3, h_cluSize_vs_cluCharge_90X3, h_cluSize_vs_cluCharge_120X3);
for (auto *c : VClustersX4) fillSizeChargeByAngle(angleX, c, h_cluSize_vs_cluCharge_60X4, h_cluSize_vs_cluCharge_90X4, h_cluSize_vs_cluCharge_120X4);

// Y clusters 60-120 degrees (separate by layer)
for (auto *c : VClustersY1) fillSizeChargeByAngle(angleY, c, h_cluSize_vs_cluCharge_60Y1, h_cluSize_vs_cluCharge_90Y1, h_cluSize_vs_cluCharge_120Y1);
for (auto *c : VClustersY2) fillSizeChargeByAngle(angleY, c, h_cluSize_vs_cluCharge_60Y2, h_cluSize_vs_cluCharge_90Y2, h_cluSize_vs_cluCharge_120Y2);
for (auto *c : VClustersY3) fillSizeChargeByAngle(angleY, c, h_cluSize_vs_cluCharge_60Y3, h_cluSize_vs_cluCharge_90Y3, h_cluSize_vs_cluCharge_120Y3);
for (auto *c : VClustersY4) fillSizeChargeByAngle(angleY, c, h_cluSize_vs_cluCharge_60Y4, h_cluSize_vs_cluCharge_90Y4, h_cluSize_vs_cluCharge_120Y4);
*/

// --- Fill Y histograms (60-120 degrees) ---
if (angleY >= 60.0 && angleY <= 120.0) {
    // Total size (sum of 4 layers)
    if (passY_gt0) {
        h_cluSize_vs_AngleY_gt0->Fill(angleY, nclusY[0]);
        if (nclusY[0] > 0) for (auto *c : VClustersY1) h_cluCharge_vs_AngleY_gt0->Fill(angleY, c->Charge());
    }
    if (passY_gt0) {
        h_cluSize_vs_AngleY_gt0->Fill(angleY, nclusY[1]);
        if (nclusY[1] > 0) for (auto *c : VClustersY2) h_cluCharge_vs_AngleY_gt0->Fill(angleY, c->Charge());
    }
    if (passY_gt0) {
        h_cluSize_vs_AngleY_gt0->Fill(angleY, nclusY[2]);
        if (nclusY[2] > 0) for (auto *c : VClustersY3) h_cluCharge_vs_AngleY_gt0->Fill(angleY, c->Charge());
    }
    if (passY_gt0) {
        h_cluSize_vs_AngleY_gt0->Fill(angleY, nclusY[3]);    
        if (nclusY[3] > 0) for (auto *c : VClustersY4) h_cluCharge_vs_AngleY_gt0->Fill(angleY, c->Charge());
    }
    if (passY_gt1) {
        h_cluSize_vs_AngleY_gt1->Fill(angleY, nclusY[0]);
        if (nclusY[0] > 1) for (auto *c : VClustersY1) h_cluCharge_vs_AngleY_gt1->Fill(angleY, c->Charge());
    }
    if (passY_gt1) {
        h_cluSize_vs_AngleY_gt1->Fill(angleY, nclusY[1]);    
        if (nclusY[1] > 1) for (auto *c : VClustersY2) h_cluCharge_vs_AngleY_gt1->Fill(angleY, c->Charge());
    }
    if (passY_gt1){
        h_cluSize_vs_AngleY_gt1->Fill(angleY, nclusY[2]);        
        if (nclusY[2] > 1) for (auto *c : VClustersY3) h_cluCharge_vs_AngleY_gt1->Fill(angleY, c->Charge());
    }
    if (passY_gt1) {
        h_cluSize_vs_AngleY_gt1->Fill(angleY, nclusY[3]);    
        if (nclusY[3] > 1) for (auto *c : VClustersY4) h_cluCharge_vs_AngleY_gt1->Fill(angleY, c->Charge());
    }
    if (passY_gt2) {
        h_cluSize_vs_AngleY_gt2->Fill(angleY, nclusY[0]);
        if (nclusY[0] > 2) for (auto *c : VClustersY1) h_cluCharge_vs_AngleY_gt2->Fill(angleY, c->Charge());
    }
    if (passY_gt2) {
        h_cluSize_vs_AngleY_gt2->Fill(angleY, nclusY[1]);
        if (nclusY[1] > 2) for (auto *c : VClustersY2) h_cluCharge_vs_AngleY_gt2->Fill(angleY, c->Charge());
    }    
    if (passY_gt2) {
        h_cluSize_vs_AngleY_gt2->Fill(angleY, nclusY[2]);
        if (nclusY[2] > 2) for (auto *c : VClustersY3) h_cluCharge_vs_AngleY_gt2->Fill(angleY, c->Charge());
    } 
    if (passY_gt2) {
        h_cluSize_vs_AngleY_gt2->Fill(angleY, nclusY[3]);
        if (nclusY[3] > 2) for (auto *c : VClustersY4) h_cluCharge_vs_AngleY_gt2->Fill(angleY, c->Charge());
    }

    // Per-layer (individual checks)
    if (nclusY[0] > 0) h_cluSize_vs_AngleY1_gt0->Fill(angleY, nclusY[0]);
    if (nclusY[1] > 0) h_cluSize_vs_AngleY2_gt0->Fill(angleY, nclusY[1]);
    if (nclusY[2] > 0) h_cluSize_vs_AngleY3_gt0->Fill(angleY, nclusY[2]);
    if (nclusY[3] > 0) h_cluSize_vs_AngleY4_gt0->Fill(angleY, nclusY[3]);

    if (nclusY[0] > 1) h_cluSize_vs_AngleY1_gt1->Fill(angleY, nclusY[0]);
    if (nclusY[1] > 1) h_cluSize_vs_AngleY2_gt1->Fill(angleY, nclusY[1]);
    if (nclusY[2] > 1) h_cluSize_vs_AngleY3_gt1->Fill(angleY, nclusY[2]);
    if (nclusY[3] > 1) h_cluSize_vs_AngleY4_gt1->Fill(angleY, nclusY[3]);

    if (nclusY[0] > 2) h_cluSize_vs_AngleY1_gt2->Fill(angleY, nclusY[0]);
    if (nclusY[1] > 2) h_cluSize_vs_AngleY2_gt2->Fill(angleY, nclusY[1]);
    if (nclusY[2] > 2) h_cluSize_vs_AngleY3_gt2->Fill(angleY, nclusY[2]);
    if (nclusY[3] > 2) h_cluSize_vs_AngleY4_gt2->Fill(angleY, nclusY[3]);

     // Per-layer charge fills (gt0, gt1, gt2)
    if (nclusY[0] > 0) for (auto *c : VClustersY1) h_cluCharge_vs_AngleY1_gt0->Fill(angleY, c->Charge());
    if (nclusY[1] > 0) for (auto *c : VClustersY2) h_cluCharge_vs_AngleY2_gt0->Fill(angleY, c->Charge());
    if (nclusY[2] > 0) for (auto *c : VClustersY3) h_cluCharge_vs_AngleY3_gt0->Fill(angleY, c->Charge());
    if (nclusY[3] > 0) for (auto *c : VClustersY4) h_cluCharge_vs_AngleY4_gt0->Fill(angleY, c->Charge());

    if (nclusY[0] > 1) for (auto *c : VClustersY1) h_cluCharge_vs_AngleY1_gt1->Fill(angleY, c->Charge());
    if (nclusY[1] > 1) for (auto *c : VClustersY2) h_cluCharge_vs_AngleY2_gt1->Fill(angleY, c->Charge());
    if (nclusY[2] > 1) for (auto *c : VClustersY3) h_cluCharge_vs_AngleY3_gt1->Fill(angleY, c->Charge());
    if (nclusY[3] > 1) for (auto *c : VClustersY4) h_cluCharge_vs_AngleY4_gt1->Fill(angleY, c->Charge());

    if (nclusY[0] > 2) for (auto *c : VClustersY1) h_cluCharge_vs_AngleY1_gt2->Fill(angleY, c->Charge());
    if (nclusY[1] > 2) for (auto *c : VClustersY2) h_cluCharge_vs_AngleY2_gt2->Fill(angleY, c->Charge());
    if (nclusY[2] > 2) for (auto *c : VClustersY3) h_cluCharge_vs_AngleY3_gt2->Fill(angleY, c->Charge());
    if (nclusY[3] > 2) for (auto *c : VClustersY4) h_cluCharge_vs_AngleY4_gt2->Fill(angleY, c->Charge());
}



        }
       
        // ==========================================
            // --- ML DATA BRIDGE: EXPORT THE EVENT ---
            // ==========================================
            // Only export if this was a valid 4-layer track on the second pass
            if (useEmpirical && clpev1x>0 && clpev2x>0 && clpev3x>0 && clpev4x>0 && 
                                clpev1y>0 && clpev2y>0 && clpev3y>0 && clpev4y>0) {
                if (chi2 < 1000){
                    
                // Write Event ID
                ml_data << jentry << ",";
                
                // Write Layer 1 (X, Y, Z, Total Charge, Max Cluster Size)
                ml_data << xpos[0] << "," << ypos[0] << "," << zpos[0] << "," 
                        << VClustersX1.at(0)->Charge() + VClustersY1.at(0)->Charge() << "," 
                        << std::max(VClustersX1.at(0)->size(), VClustersY1.at(0)->size()) << ",";
                
                // Write Layer 2
                ml_data << xpos[1] << "," << ypos[1] << "," << zpos[1] << "," 
                        << VClustersX2.at(0)->Charge() + VClustersY2.at(0)->Charge() << "," 
                        << std::max(VClustersX2.at(0)->size(), VClustersY2.at(0)->size()) << ",";
                
                // Write Layer 3
                ml_data << xpos[2] << "," << ypos[2] << "," << zpos[2] << "," 
                        << VClustersX3.at(0)->Charge() + VClustersY3.at(0)->Charge() << "," 
                        << std::max(VClustersX3.at(0)->size(), VClustersY3.at(0)->size()) << ",";
                
                // Write Layer 4
                ml_data << xpos[3] << "," << ypos[3] << "," << zpos[3] << "," 
                        << VClustersX4.at(0)->Charge() + VClustersY4.at(0)->Charge() << "," 
                        << std::max(VClustersX4.at(0)->size(), VClustersY4.at(0)->size()) << ",";
                
                // Write Angles
                ml_data << angleX << "," << angleY << "\n";
            }
            
        }

        for(auto c : VClustersX1) delete c; VClustersX1.clear();
        for(auto c : VClustersX2) delete c; VClustersX2.clear();
        for(auto c : VClustersX3) delete c; VClustersX3.clear();
        for(auto c : VClustersX4) delete c; VClustersX4.clear();
        for(auto c : VClustersY1) delete c; VClustersY1.clear();
        for(auto c : VClustersY2) delete c; VClustersY2.clear();
        for(auto c : VClustersY3) delete c; VClustersY3.clear();           
        for(auto c : VClustersY4) delete c; VClustersY4.clear();
        
        HT1X.clear(); HT2X.clear(); HT3X.clear(); HT4X.clear();
        HT1Y.clear(); HT2Y.clear(); HT3Y.clear(); HT4Y.clear();
        C1X.Epmty(); C2X.Epmty(); C3X.Epmty(); C4X.Epmty(); C1Y.Epmty(); C2Y.Epmty(); C3Y.Epmty(); C4Y.Epmty();
        VClustersX1.clear(); VClustersX2.clear(); VClustersX3.clear(); VClustersX4.clear();
        VClustersY1.clear(); VClustersY2.clear(); VClustersY3.clear(); VClustersY4.clear();
        event_typeX=0;  event_typeY=0;
        clpev1x=0; clpev1y=0; clpev2x=0; clpev2y=0;
        clpev3x=0; clpev3y=0; clpev4x=0; clpev4y=0;
       //////
        countstrips=0;
        pos1X=0; pos2X=0; pos1Y=0; pos2Y=0;pos3X=0; pos4X=0; pos3Y=0; pos4Y=0;
        slope14=0; slope24=0; angle24=0; angleX=0; angleY=0;
        slope12=0; slope23=0; angle12=0; angle23=0;

    
        if(jentry%10000==0) 	std::cout <<  jentry <<endl;
        entries=jentry;
        
    }      //Loop over the events
    
    
    
    
    std::cout << "All events = " <<entries <<endl;

if (!useEmpirical) {
    // First pass: extract empirical sigmas from residual histograms
    h_res_vs_sizeX->FitSlicesY();
    h_res_vs_sizeY->FitSlicesY();
    UpdateEmpiricalSigmas();
    std::cout << "First pass completed. Empirical sigmas extracted.\n";
    std::cout << "Now call Loop(true) to run the second pass.\n";
} else {
    std::cout << "Second pass with empirical errors completed.\n";
}

    plot_Raw();
//    plot_Clusters();
    plot_Residuals();
// Close the ML export file safely
    if (useEmpirical) {
        ml_data.close();
        std::cout << "\n---> ML DATASET GENERATED: 'ml_real_data.csv' <---" << std::endl;
    }

std::cout << "\n==========================================" << endl;
std::cout << "      AUTOMATED LOOP ALIGNMENT RESULTS      " << endl;
std::cout << "==========================================" << endl;

    double best_sigma = 999.0;
    double best_shift = 0.0;
    double best_chi2 = 999.0;

    for(int i=0; i<=99; i++) {
        // Fit the sharp core of each shifted histogram
        scan_hist[i]->Fit("gaus", "Q", "", -8.0, 8.0);
        TF1 *fit = scan_hist[i]->GetFunction("gaus");
        if(fit) {
            double current_sigma = fit->GetParameter(2);
            double current_chi2 = fit->GetChisquare() / fit->GetNDF();
            
            // If this specific shift gave a tighter track (smaller sigma), save it!
            if(current_sigma < best_sigma) {
                best_sigma = current_sigma;
                best_chi2 = current_chi2;
                best_shift = -0.5 + (i * 0.05);
            }
        }
    }

    // ==========================================================
    // --- DRAW 2D RESIDUAL VS STRIP CANVASES ---
    // ==========================================================


    TCanvas *h_res_vs_sizeXY = new TCanvas("h_res_vs_sizeXY", "Res vs SizeXY", 1200, 800);
    h_res_vs_sizeXY->Divide(1,2);
    h_res_vs_sizeXY->cd(1); h_res_vs_sizeX->Draw("colz");
    h_res_vs_sizeXY->cd(2); h_res_vs_sizeY->Draw("colz");
    
    

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

    std::cout << "--> The absolute BEST corr_pos2X is : " << best_shift << " mm" << std::endl;
    std::cout << "--> This achieved a Sigma of        : " << best_sigma << " mm" << std::endl;
    std::cout << "--> This achieved a Chi2/ndf of     : " << best_chi2 << std::endl;
    std::cout << "==========================================\n" << std::endl;

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

// Set your color scale from 0 to 10 (or whatever max error you want)
p_master_graph->SetMinimum(0.0);
p_master_graph->SetMaximum(10.0);

p_master_graph->Draw("colz");
// Right side: The Profile (The mathematical average Chi2 vs width)
c_delta_angle->cd(2);
p_avgChi2_vs_deltaAngleX->SetMarkerStyle(20);
p_avgChi2_vs_deltaAngleX->Draw("P");

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
*/
//   TCanvas *c200 = new TCanvas("c200"," MAP ", 150,0,800,800);
//    TGraph* graph0 = new TGraph(mi1, mapX1, mapY1);

  //  graph0->Draw("A*");
    
 //   h_map_XY->Draw();
 //   cout<< "mi1 = "<<mi1 <<endl;
    
}
void mmanalys::Calibrate() {
    // First pass with analytic errors
    Loop(false);
    // Second pass with empirical errors
    Loop(true);
}
void residual() {
    mmanalys t;
    t.Calibrate();
}

