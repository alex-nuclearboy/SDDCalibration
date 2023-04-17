// XrayLines.h

#ifndef XRAYLINES_H
#define XRAYLINES_H

// "_RI" refers to Relative Intensity

/***** Energy of electron-hole pair creation in silicon *****/
double SiW = 3.71;

/********* Silicon *********/
// K-alpha 
static const double SiKa        = 1739. ;

/*********** Iron **********/
// K-alpha1 
static const double FeKa1       = 6403.84;
static const double FeKa1_RI    = 100 ;
// K-alpha2 
static const double FeKa2       = 6390.84;
static const double FeKa2_RI    = 50  ;
// K-beta1,3
static const double FeKb1_3     = 7058.0;
static const double FeKb1_3_RI  = 17.;
// Gamma 1-0
static const double FeG1_0      = 14412.95;
static const double FeG1_0_RI   = 87.67; // Internal conversion coefficients, percent 

/********* Chromium ********/
// K-alpha1
static const double CrKa1       = 5414.7;
static const double Crka1_RI    = 100;
// K-alpha2
static const double CrKa2       = 5405.5;
static const double CrKa2_RI    = 50;
// K-beta1,3
static const double CrKb13      = 5946.7;
static const double CrKb13_RI   = 15;

/********* Calcium *********/
// K-alpha1 
static const double CaKa1       = 3691.7;
static const double CaKa1_RI    = 100 ;
// K-alpha2 
static const double CaKa2       = 3688.1;
static const double CaKa2_RI    = 50  ;
// K-beta1,3
static const double CaKb1_3     = 4012.7;
static const double CaKb1_3_RI  = 13.;

/********** Copper *********/
// K-alpha1
static const double CuKa1       = 8047.78;
static const double CuKa1_W     = 2.11;
static const double CuKa1_RI    = 100 ;
// K-alpha2
static const double CuKa2       = 8027.83;
static const double CuKa2_W     = 2.17;
static const double CuKa2_RI    = 51  ;
// K-beta1
static const double CuKb1       = 8905.29;
static const double CuKb1_RI    = 17  ;

/*********** Zinc **********/
// K-alpha1 
static const double ZnKa1       = 8638.86;
static const double ZnKa1_RI    = 100 ;
// K-alpha2
static const double ZnKa2       = 8615.78;
static const double ZnKa2_RI    = 51  ;
// K-beta1
static const double ZnKb1       = 9572.0;
static const double ZnKb1_RI    = 17 ;
// L-alpha1
static const double ZnLa1       = 1011.7;
// L-alpha2
static const double ZnLa2       = 1011.7;

/********* Titanium ********/
// K-alpha1
static const double TiKa1       = 4510.84;
static const double TiKa1_W     = 1.16;
static const double TiKa1_RI    = 100;
// K-alpha2
static const double TiKa2       = 4504.86;
static const double TiKa2_W     = 1.18;
static const double TiKa2_RI    = 50 ;
// K-beta1
static const double TiKb1       = 4931.81;
static const double TiKb1_RI    = 15 ;

/******** Manganese ********/
// K-alpha1 
static const double MnKa1       = 5898.75;
static const double MnKa1_RI    = 100 ;
// K-alpha2 
static const double MnKa2       = 5887.65;
static const double MnKa2_RI    = 50  ;
// K-beta1 
static const double MnKb1       = 6490.45;
static const double MnKb1_RI    = 17  ;

/********** Nickel *********/ 
// K-alpha1 
static const double NiKa1       = 7478.15;
static const double NiKa1_RI    = 100;

static const double NiKa2       = 7460.89;
static const double NiKa2_RI    = 51;

static const double NiKb1_3     = 8264.7; 
static const double NiKb1_3_RI  = 17;

/********** Aurum  *********/
// L-alpha1
static const double AuLa1       = 9713.3;
static const double AuLa1_RI    = 100.;
// L-alpha2
static const double AuLa2       = 9628.0;
static const double AuLa2_RI    = 11.;
// L-beta 1
static const double AuLb1       = 11442.3;
static const double AuLb1_RI    = 67  ;
// L-beta 2
static const double AuLb2       = 11584.7;
static const double AuLb2_RI    = 23  ;

/********* Tungsten ********/
// L-alpha1 
static const double WLa1        = 8397.6;
static const double WLa1_RI     = 100.;
// L-alpha2 
static const double WLa2        = 8335.2;
static const double WLa2_RI     = 11.;
// L-beta1 
static const double WLb1        = 9672.4;
static const double WLb1_RI     = 67.;
// L-beta2 
static const double WLb2        = 9961.5;
static const double WLb2_RI     = 21.;
// L-gamma1
static const double WLg1        = 11285.9;
static const double WLg1_RI     = 13.;

/******** Zirconium ********/
// K-alpha1
static const double ZrKa1       = 15775.1;
static const double ZrKa1_RI    = 100 ;
// K-alpha2 
static const double ZrKa2       = 15690.9;
static const double ZrKa2_RI    = 51  ;
// K-beta1
static const double ZrKb1       = 17667.8;
static const double ZrKb1_RI    = 15  ;
// K-beta2
static const double ZrKb2       = 17970  ;
static const double ZrKb2_RI    = 3   ;
// K-beta3
static const double ZrKb3       = 17654  ;
static const double ZrKb3_RI    = 8   ;

/*********** Lead **********/
// L-alpha1
static const double PbLa1       = 10551.5;
static const double PbLa1_RI    = 100.;
// L-alpha2
static const double PbLa2       = 10449.5;
static const double PbLa2_RI    = 11.;
// L-beta1
static const double PbLb1       = 12613.7;
static const double PbLb1_RI    = 66.;
// L-beta2
static const double PbLb2       = 12622.6;
static const double PbLb2_RI    = 25.;

/********* Bromine *********/
// K-alpha1
static const double BrKa1       = 11924.2;
static const double BrKa1_RI    = 100.;
// K-alpha2
static const double BrKa2       = 11877.6;
static const double BrKa2_RI    = 52.;
// K-beta1
static const double BrKb1       = 13291.4;
static const double BrKb1_RI    = 14.;
// K-beta2
static const double BrKb2       = 13469.;
static const double BrKb2_RI    = 1.36;
// K-beta3
static const double BrKb3       = 13284.5;
static const double BrKb3_RI    = 7.;

/******** Strontium ********/
// K-alpha1
static const double SrKa1       = 14165.;
static const double SrKa1_RI    = 100.;
// K-alpha2
static const double SrKa2       = 14098.;
static const double SrKa2_RI    = 52.07;
// K-beta1
static const double SrKb1       = 15836.;
static const double SrKb1_RI    = 14.4;
// K-beta2
static const double SrKb2       = 16085.;
static const double SrKb2_RI    = 2.56;
// K-beta3
static const double SrKb3       = 15825.;
static const double SrKb3_RI    = 7.44;

/******** Cobalt 57 ********/
// Gamma 1, 0 Fe
static const double CoGamma10   = 14412.95;

/********** Silver *********/
// K-alpha1
static const double AgKa1       = 22162.9; 
static const double AgKa1_RI    = 100.;
// K-alpha2
static const double AgKa2       = 21990.3; 
static const double AgKa2_RI    = 53.;
// K-beta1
static const double AgKb1       = 24942.4;
static const double AgKb1_RI    = 16.;
// K-beta3
static const double AgKb3       = 24911.5;
static const double AgKb3_RI    = 8.;

/********* Antimony ********/   // In alloy with Lead
// K-alpha1
static const double SbKa1       = 26359.1;
static const double SbKa1_RI    = 100.;
// K-alpha2 
static const double SbKa2       = 26110.8;
static const double SbKa2_RI    = 53.;

/*********** Tin ***********/
// K-alpha1 
static const double SnKa1       = 25271.3;
static const double SnKa1_RI    = 100.;
// K-alpha2
static const double SnKa2       = 25044.0;
static const double SnKa2_RI    = 53.;

/********* Bismuth *********/
// M-alpha1
static const double BiMa1       = 2422.6; 
static const double BiMa1_RI    = 100.;
// L-lines
static const double BiLa2       = 10730.91; 
static const double BiLa2_RI    = 11.;
static const double BiLa1       = 10838.8;
static const double BiLa1_RI    = 100.;
static const double BiLb2       = 12979.9; 
static const double BiLb2_RI    = 25.;
static const double BiLb1       = 13023.5;
static const double BiLb1_RI    = 67.;
static const double BiLg1       = 15247.7;
static const double BiLg1_RI    = 14.;

/******** Palladium ********/   // In the electronics and soldering material
// K-alpha1
static const double PdKa1       = 21177.1;
static const double PdKa1_RI    = 100.;
// K-alpha2
static const double PdKa2       = 21020.1;
static const double PdKa2_RI    = 53.;
// K-beta1 
static const double PdKb1       = 23818.7;
static const double PdKb1_RI    = 16;
// K-beta2 
static const double PdKb2       = 24299.1;
static const double PdKb2_RI    = 4.;
// K-beta3 
static const double PdKb3       = 23791.1;
static const double PdKb3_RI    = 8.;

/****** Kaonic Carbon ******/
static const double KC54        = 10216.5 ;
static const double KC54_RI     = 100.;
static const double KC65        = 5544.9 ;
static const double KC65_RI     = 22.68;
static const double KC64        = 15759.4 ;
static const double KC75        = 8885.8 ;
static const double KC75_RI     = 2.62;
static const double KC86        = 5509.6 ;
static const double KC86_RI     = 3.09;

/***** Kaonic Nitrogen *****/
static const double KN54        = 13995.9 ;
static const double KN65        = 7595.4 ;
static const double KN76        = 4577.3 ;

/****** Kaonic Oxygen ******/
static const double KO65        = 9968.7 ;
static const double KO65_RI     = 100.;
static const double KO76        = 6006.8 ;
static const double KO76_RI     = 27.8;
static const double KO75        = 15973.3 ;
static const double KO75_RI     = 44.4;
static const double KO86        = 9902.7 ;
static const double KO86_RI     = 16.7;

/***** Kaonic Titanium *****/
static const double KTi109      = 14745.6 ;
static const double KTi1110     = 10910.1 ;
static const double KTi1211     = 8298.02;
static const double KTi1210     = 19208.1 ;
static const double KTi1312     = 6457.81;
static const double KTi1311     = 14755.8 ;
static const double KTi1413     = 5124.08;

/***** Kaonic Aluminum *****/
static const double KAl21       = 1634000;
static const double KAl32       = 303000;
static const double KAl43       = 106000;
static const double KAl54       = 49000;
static const double KAl65       = 26600;
static const double KAl76       = 16040.1 ;
static const double KAl87       = 10410.6 ;
static const double KAl98       = 7144.0 ;
static const double KAl109      = 5110.0;
static const double KAl1110     = 3781.0;
static const double KAl1211     = 2876.0;
static const double KAl1312     = 2238.0;
static const double KAl1413     = 1776.0;
static const double KAl1514     = 1433.0;

/****** Kaonic Helium ******/
static const double KHe32       = 6463.5 ;
static const double KHe42       = 8721.7 ;
static const double KHe52       = 9766.8 ;

#endif
