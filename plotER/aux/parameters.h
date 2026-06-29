


double BR_B0s_jpsiphi_mumuKK = 3.04*10e-5;
double BR_Bp_jpsiK_mumuK = 6.08*10e-5;

int nbinsmasshisto = 80;

double minhisto = 0;
double maxhisto = 999;
double minhisto_B=5.;
double maxhisto_B=5.8;
double minhisto_X=3.6;
double maxhisto_X=4.0;

const int N_pt_Bins_X = 4;
std::vector<double> ptbinsvec_X = {7.5, 12.5, 17.5, 22.5, 50};

///FOR TESTING TESTING
//const int N_pt_Bins_X = 1;
//std::vector<double> ptbinsvec_X = {15, 25};
///FOR TESTING TESTING

const int N_pt_Bins_B = 5;
std::vector<double> ptbinsvec_B = {5, 10, 15, 20, 30, 60};

const int N_y_Bins_X = 4;
std::vector<double> ybinsvec = {0.0, 0.8, 1.5, 2.0, 2.4};

const int N_mult_Bins_X = 4;
std::vector<double> nmbinsvec = {0, 15, 30, 50, 100};

const int N_cent_Bins_X = 3;
std::vector<double> centbinsvec = {0, 20, 40, 80};