#include <bits/stdc++.h>

using namespace std;

#define SCALE 1
#define INV_SCALE 1/SCALE
#define QTD_DIFFUSIONS 4
#define DIM_X 15
#define DIM_Y 5
#define DIM_Z 3
#define ALPHA 0.01
#define PI 3.14159265359
#define R 8.314472 // [N.m/(g.mol.K)] (Gas constant) Zaki Ahmad, in Principles of Corrosion Engineering and Corrosion Control, 2006
#define T 298.16 // [K] Ambient temperature
#define n_NA 1 // Sodium charge
#define n_CA 2 // Calcium charge
#define F 96485.3329 // [s.A/mol] Faraday's constant
#define low_threshold 0.0
#define high_threshold 0.01

double conc, tx_conc, diameter_cell = 5, deltaamp;

// #### GAP JUNCTIONS PROBABILITIES FOR ASTROCYTES #### //
double phl[50] = {3.33333333333e-01,9.51756626745e-02,2.71812917035e-02,7.76661124715e-03,2.22288659849e-03,6.39921679991e-04,1.87410416804e-04,5.81750667195e-05,2.17596143238e-05,1.1436923158e-05,7.88454209682e-06,7.43738619183e-06,7.37970057786e-06,7.29603316347e-06,7.27478942971e-06,7.26006289992e-06,7.26084787208e-06,7.26080132601e-06,7.26061054996e-06,7.26081361742e-06,7.26079620991e-06,7.26072365567e-06,7.26058079345e-06,7.26074725419e-06,7.26087576894e-06,7.26073008288e-06,7.26061194028e-06,7.26074336727e-06,7.26101528686e-06,7.26081974085e-06,7.26091667847e-06,7.26058059059e-06,7.26084014577e-06,7.2610063969e-06,7.26069065682e-06,7.26083092741e-06,7.26076153595e-06,7.26071756287e-06,7.26092535023e-06,7.26076421324e-06,7.26060026219e-06,7.26075209967e-06,7.26093367537e-06,7.26073986493e-06,7.26039032094e-06,7.26091299989e-06,7.26077756319e-06,7.26071491915e-06,7.2607710224e-06,7.26082337127e-06};
double plh[50] = {0.333333333333,0.706523083189,0.825908352326,0.86117728508,0.871356970354,0.874272895353,0.875107602476,0.875346523333,0.875414979813,0.875433066086,0.875439533589,0.875437420808,0.875439605921,0.875443525702,0.875437685437,0.875440602592,0.875441218644,0.875440938251,0.875438148928,0.875441400815,0.875441147657,0.875440025483,0.875437801975,0.875440355641,0.875442338151,0.875440081279,0.87543825197,0.875440285023,0.87544449347,0.875441466845,0.875442967173,0.875437765288,0.875441782611,0.875444355802,0.875439468868,0.875441639938,0.875440565916,0.875439885313,0.875443101387,0.875440607355,0.875438069767,0.875440419864,0.875443230241,0.875440230498,0.875434820356,0.875442910232,0.875440813981,0.875439844395,0.875440712745,0.875441522986};
double phh[50] = {0.333333333333,0.198301254137,0.14691035597,0.131056103673,0.126420143048,0.125087182967,0.124704987107,0.1245953016,0.124563260573,0.124555496991,0.124552581869,0.124555141806,0.124553014379,0.124549178265,0.124555039774,0.124552137345,0.124551520508,0.124551800947,0.124554590462,0.124551338371,0.124551591547,0.124552713793,0.124554937445,0.124552383612,0.124550400973,0.124552657991,0.124554487418,0.124552454233,0.124548245514,0.124551272335,0.12454977191,0.124554974131,0.124550956549,0.124548383192,0.124553270441,0.124551099232,0.124552173322,0.124552853969,0.124549637687,0.124552131881,0.124554669633,0.124552319384,0.124549508825,0.124552508762,0.124557919254,0.124549828855,0.124551925241,0.12455289489,0.124552026484,0.124551216191};

// Classe para representar todas as variáveis da célula (Ex.: astrócito)
class Cell {
public:
	map<string, double> parameters;
	int id;

	Cell() {
		/*
			0v1: taxa de liberação de Cálcio do Tx
			1Y:
			2vin: o fluxo de calcio que parte do espaço extracelular, por meio da membrana do astrócito, até o interior do citosol
			3VM2: o fluxo máximo de íons de cálcio fora da bomba (???)
			4C: Concentração de Cálcio no citosol
			5n: Coeficiente de Hill (2.02)
			6K2:
			10kf: constante que determina a liberação de cálcio do RE para o citosol
			22D: coneficiente de difusão
			23l: volume da célula
			24K: taxa máxima da ativação do receptor (Nakano, 2010; Eq. 3)
			25ka: taxa máxima da ativação do receptor (2.5) (Nakano, 2010; Eq. 3)
			26m: Coeficiente de Hill (2.2)
			27phh - 29plh: probabilidades das gap junctions
		*/

		double v1 = 100*SCALE; parameters["v1"] = v1;
		double Y = 0*SCALE; parameters["Y"] = Y;
		double vin = 0.05*SCALE; parameters["vin"] = vin;
		double VM2 = 15*SCALE; parameters["VM2"] = VM2;
		double C = 0.1*SCALE; parameters["C"] = C;
		double n = 2.02*SCALE; parameters["n"] = n;
		double K2 = 0.1*SCALE; parameters["K2"] = K2;
		double VM3 = 40*SCALE; parameters["VM3"] = VM3; // Por que nao 40, como diz no artigo?
		double ko = 0.5*SCALE; parameters["ko"] = ko;
		double ER = 1.5*SCALE; parameters["ER"] = ER;
		double kf = 0.5*SCALE; parameters["kf"] = kf;
		double kp = 0.3*SCALE; parameters["kp"] = kp;
		double kdeg = 0.08*SCALE; parameters["kdeg"] = kdeg;
		double vp = 0.05*SCALE; parameters["vp"] = vp;
		double kcaaa = 0.15*SCALE; parameters["kcaaa"] = kcaaa;
		double kcai = 0.15*SCALE; parameters["kcai"] = kcai;
		double kip3 = 0.1*SCALE; parameters["kip3"] = kip3;
		double IP3 = 0.1*SCALE; parameters["IP3"] = IP3;
		double q = 4*SCALE; parameters["q"] = q;
		double W = 0*SCALE; parameters["W"] = W;
		double A = 0*SCALE; parameters["A"] = A;
		double kia = 0.5*SCALE; parameters["kia"] = kia;
		double D = 350 * (4 * PI * diameter_cell / 3.0)/*122500*/; parameters["D"] = D;
		double l = PI *SCALE * pow(diameter_cell / 2, 2); parameters["l"] = l;
		double K = 0.0006 *SCALE; parameters["K"] = K;
		double ka = 2.5*SCALE; parameters["ka"] = ka;
		double m = 2.2*SCALE; parameters["m"] = m;
		parameters["phl"] = phl[0];
		parameters["plh"] = plh[0];
		parameters["phh"] = phh[0];

		double Na_i = SCALE*15000;  parameters["Na_i"] = Na_i; //Langer3, Chatton 2016 [uM] - Colocar uma relação
		double Na_o = SCALE*150000/(DIM_X * DIM_Y * DIM_Z);  parameters["Na_o"] = Na_o; //chatton 2016 [uM]
		double NaD = (4/3)*PI*diameter_cell*600; parameters["NaD"] = NaD; 
		double C_o = SCALE*2300/(DIM_X * DIM_Y * DIM_Z); parameters["C_o"] = C_o; //Kirischuk 1997 [uM]
		//double K_i = 0; parameters["K_i"] = K_i; // A definir
		//double K_o = 0; parameters["K_o"] = K_o; // A definir
		double V_NCX = -SCALE*80; parameters["V_NCX"] = V_NCX; //Verify this value later - Koenigsberger, 2004-->-40/Kirischuk 2012-->-80 [mV]
		double Vm = SCALE*R*T/(n_CA*F)*log(C_o/C)*1000; parameters["Vm"] = Vm; // [mV] VERIFICAR!!!! Nernst Eq./Eq (3) from  Koenigsberger, 2004
		double DDC_Na = -SCALE*115000; parameters["DDC_Na"] = DDC_Na; // Kirischuk 2012 [uM]
		double G_NCX = SCALE*0.00316; parameters["G_NCX"] = G_NCX; // Koenigsberger, 2004 [uM.m.V^(-1).S^(-1)]
		double C_NCX = SCALE*0.5; parameters["C_NCX"] = C_NCX; // Barros, 2015 [uM]
		// EQUAÇÕES PARALELAS
		double K_mCaAct = SCALE*0.394; parameters["K_mCaAct"] = K_mCaAct; // Weber, 2001 & Brazhe 2018 [uM]
		double n_Hill = SCALE*2; parameters["n_Hill"] = n_Hill; // Weber, 2001 & Brazhe 2018
		double J_max = SCALE*25; parameters["J_max"] = J_max; // Brazhe 2018 [A/F] 25-60
		double eta = SCALE*0.35; parameters["eta"] = eta; // Weber, 2001 & Brazhe 2018
		double k_sat = SCALE*0.27; parameters["k_sat"] = k_sat; // Weber, 2001 & Brazhe 2018
		double K_mNao = SCALE*87500; parameters["K_mNao"] = K_mNao; // Weber, 2001 & Brazhe 2018 [uM]
		double K_mNai = SCALE*12300; parameters["K_mNai"] = K_mNai; // Weber, 2001 & Brazhe 2018 [uM]
		double K_mCao = SCALE*1300; parameters["K_mCao"] = K_mCao; // Weber, 2001 & Brazhe 2018 [uM]
		double K_mCai = SCALE*3.6; parameters["K_mCai"] = K_mCai; // Weber, 2001 & Brazhe 2018 [uM]
		double H_Na = SCALE*2; parameters["H_Na"] = H_Na; // Matsuoka 1996 & Brazhe 2018 ou 1 ***********?????????????
		double H_Ca = SCALE*1.2; parameters["H_Ca"] = H_Ca; // Matsuoka 1996 & Brazhe 2018 ou 1.1 ***********?????????????
		double K_Na = SCALE*28000; parameters["K_Na"] = K_Na; // Matsuoka 1996 & Brazhe 2018 ou 1 [uM] ***********?????????????
		double K_Ca = SCALE*0.1; parameters["K_Ca"] = K_Ca; // Matsuoka 1996 & Brazhe 2018 ou 1.4 [uM] ***********?????????????
		double tau_0 = SCALE*0.1; parameters["tau_0"] = tau_0; // Matsuoka 1996 & Brazhe 2018 [s] ***********?????????????
		double H_tau = SCALE*0.5; parameters["H_tau"] = H_tau; // Matsuoka 1996 & Brazhe 2018 ***********?????????????
		double K_tau = SCALE*0.5; parameters["K_tau"] = K_tau; // Matsuoka 1996 & Brazhe 2018 ***********?????????????
		double tau_Ca = SCALE*0.001; parameters["tau_Ca"] = tau_Ca; // Brazhe 2018 ***********?????????????
		double C_rest = SCALE*0.073; parameters["C_rest"] = C_rest; // Kirischuk 2012 [uM] 

		// Modulation and Demodulation
		double C_variation = SCALE*0; parameters["C_variation"] = C_variation;
	}

	void setId(int ID) {
		id = ID;
		parameters["id"] = ID;
	}
};

class Network {
public:
	int NC; // Number of cells
	Cell tecido[DIM_Y][DIM_X][DIM_Z];
	list<int> *connect;
	vector<int> rx_id;

	Network() {
		NC = DIM_X * DIM_Y * DIM_Z;
		connect = new list<int>[NC];

		int id_cont = 0;

		for (int k = 0; k < DIM_Z; k++) {
			for (int i = 0; i < DIM_Y; i++) {
				for (int j = 0; j < DIM_X; j++) {
					Cell celula;
					celula.setId(id_cont);
					tecido[i][j][k] = celula;
					id_cont++;
				}
			}
		}
	}

	int getId(int x, int y, int z) {
		return tecido[y][x][z].parameters["id"];
	}

	double get(int x, int y, int z, string parameter) {
		return tecido[y][x][z].parameters[parameter];
	}

	double get(int id, string parameter) {
		int x = id % DIM_X;
		int y = (id / DIM_X) % DIM_Y;
		int z = id / (DIM_X * DIM_Y);
		return tecido[y][x][z].parameters[parameter];
	}

	void set(int x, int y, int z, string parameter, double value) {
		tecido[y][x][z].parameters[parameter] = value;
	}

	void set(int id, string parameter, double value) {
		int x = id % DIM_X;
		int y = (id / DIM_X) % DIM_Y;
		int z = id / (DIM_X * DIM_Y);
		tecido[y][x][z].parameters[parameter] = value;
	}

	void accumulate(int x, int y, int z, string parameter, double add) {
		tecido[y][x][z].parameters[parameter] += add;
	}
	
	void update_parameters(int x, int y, int z) {
		tecido[y][x][z].parameters["Vm"] = R*T/(n_CA*F)*log(tecido[y][x][z].parameters["C_o"]/tecido[y][x][z].parameters["C"])*1000; // ALTERAR!!!
	}

	void accumulate(int id, string parameter, double add) {
		int x = id % DIM_X;
		int y = (id / DIM_X) % DIM_Y;
		int z = id / (DIM_X * DIM_Y);
		tecido[y][x][z].parameters[parameter] += add;
	}

	void changeSignal(int x, int y, int z, string parameter) {
		tecido[y][x][z].parameters[parameter] *= -1;
	}

	void printTissue() {
		cout << endl;
		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				cout << fixed << setprecision(2) << get(i, j, trunc(DIM_Z / 2), "C") << ":";
				cout << fixed << setprecision(2) << get(i, j, trunc(DIM_Z / 2), "C_o") << "  ";

			}
			cout << endl;
		}
		cout << endl << endl;
		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				cout << fixed << setprecision(2) << get(i, j, trunc(DIM_Z / 2), "Na_i") << ":";
				cout << fixed << setprecision(2) << get(i, j, trunc(DIM_Z / 2), "Na_o") << "  ";

			}
			cout << endl;
		}
	}

	void writeFileHeader(ofstream &file){
		for (int i = 0; i < DIM_Y; i++) {
		  for (int j = 0; j < DIM_X; j++) {
			if(i == DIM_Y-1 && j == DIM_X-1)
			  file << i << "_" << j;
			else
			  file <<  i << "_" << j << ",";
		  }
		}
		file << endl;
	}

	void printTissuef(ofstream &file, int time_int){
		file << time_int << ",";
		for (int i = 0; i < DIM_Y; i++) {
		  for (int j = 0; j < DIM_X; j++) {
			if(i == DIM_Y-1 && j == DIM_X-1)
			  file << get(j, i, trunc(DIM_Z / 2), "C");
			else
			  file <<  get(j, i, trunc(DIM_Z / 2), "C") << ",";
		  }
		}
		file << endl;
	}

	void regularDegree() {
		for (int x = 0; x < DIM_X; x++) {
			for (int y = 0; y < DIM_Y; y++) {
				for (int z = 0; z < DIM_Z; z++) {
					// x + 1
					if (y + 1 < DIM_Y) {
						connect[getId(x, y, z)].push_back(getId(x, y + 1, z));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
					// x - 1
					if (y - 1 >= 0) {
						connect[getId(x, y, z)].push_back(getId(x, y - 1, z));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
					// y + 1
					if (x + 1 < DIM_X) {
						connect[getId(x, y, z)].push_back(getId(x + 1, y, z));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
					// y - 1
					if (x - 1 >= 0) {
						connect[getId(x, y, z)].push_back(getId(x - 1, y, z));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
					// z + 1
					if (z + 1 < DIM_Z) {
						connect[getId(x, y, z)].push_back(getId(x, y, z + 1));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
					// z - 1
					if (z - 1 >= 0) {
						connect[getId(x, y, z)].push_back(getId(x, y, z - 1));
					} else {
						connect[getId(x, y, z)].push_back(-1);
					}
				}
			}
		}

		// IMPRIMINDO AS CONEXÕES DE CADA CÉLULA
		list<int>::iterator it;

		for (int v = 0; v < NC; v++) {
			cout << "Célula " << v << ": ";
			for (it = connect[v].begin(); it != connect[v].end(); it++){
				cout << *it << " ";
			}
			cout << endl;
		}
	}

	void linkRadius(int radius) {
		for (int x = 0; x < DIM_X; x++) {
			for (int y = 0; y < DIM_Y; y++) {
				for (int z = 0; z < DIM_Z; z++) {
					for (int r = 1; r <= radius; r++) {
						// x + r
						if (y + r < DIM_Y) {
							connect[getId(x, y, z)].push_back(getId(x, y + r, z));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
						// x - r
						if (y - r >= 0) {
							connect[getId(x, y, z)].push_back(getId(x, y - r, z));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
						// y + r
						if (x + r < DIM_X) {
							connect[getId(x, y, z)].push_back(getId(x + r, y, z));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
						// y - r
						if (x - r >= 0) {
							connect[getId(x, y, z)].push_back(getId(x - r, y, z));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
						// z + r
						if (z + r < DIM_Z) {
							connect[getId(x, y, z)].push_back(getId(x, y, z + r));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
						// z - r
						if (z - r >= 0) {
							connect[getId(x, y, z)].push_back(getId(x, y, z - r));
						} else {
							connect[getId(x, y, z)].push_back(-1);
						}
					}
				}
			}
		}

		// IMPRIMINDO AS CONEXÕES DE CADA CÉLULA
		list<int>::iterator it;

		for (int v = 0; v < NC; v++) {
			cout << "Célula " << v << ": ";
			for (it = connect[v].begin(); it != connect[v].end(); it++){
				cout << *it << " ";
			}
			cout << endl;
		}
	}

	vector<int> getConnections(int id) {
		vector<int> connections;
		list<int>::iterator it;

		for (it = connect[id].begin(); it != connect[id].end(); it++){
			connections.push_back(*it);
		}

		return connections;
	}

	int numberConnections() {
		return connect[0].size();
	}

	// Modulation and demodulation
	int mod_demod(int x, int y, int z){

		if (tecido[y][x][z].parameters["C_variation"]>=high_threshold){
			return 1;
		}else if (tecido[y][x][z].parameters["C_variation"]>=low_threshold && tecido[y][x][z].parameters["C_variation"]<high_threshold){
			return 0;
		}else{
			return -1;
		}
	}
	int mod_demod(int id){
		int x = id % DIM_X;
		int y = (id / DIM_X) % DIM_Y;
		int z = id / (DIM_X * DIM_Y);

		if (tecido[y][x][z].parameters["C_variation"]>=high_threshold){
			return 1;
		}else if (tecido[y][x][z].parameters["C_variation"]>=low_threshold && tecido[y][x][z].parameters["C_variation"]<high_threshold){
			return 0;
		}else{
			return -1;
		}
	}

	double conditional_accumulate(vector<int>::iterator first, vector<int>::iterator last, int value){
		double sum = 0;
		for (; first != last; first++){
			if(*first == value) sum += *first;
			// cout << "bit: " << *first << endl;
		}
		return sum;
	}

	void rx_geometry(string rx_geometry, int Rx_number, int destination, int tx_x, int tx_y, int tx_z){
		if (rx_geometry == "VL"){ // Up to DIM_X
			for (int i = 0; i < Rx_number; i++){
				rx_id.push_back(getId(tx_x+destination+i, tx_y+1, tx_z));
			}
		}else if (rx_geometry == "HL"){ // Up to DIM_Y
			for (int i = 0; i < Rx_number; i++){
				rx_id.push_back(getId(tx_x+destination, tx_y-1+i, tx_z));
			}
		}else if (rx_geometry == "X"){ // Up to five cells
			for (int i = 0; i < Rx_number; i++){
				if (i<3){
					rx_id.push_back(getId(tx_x+destination+i, tx_y+1, tx_z));
				}else if (i==3){
					rx_id.push_back(getId(tx_x+destination+1, tx_y+1, tx_z+1));
				}else if (i==4){
					rx_id.push_back(getId(tx_x+destination+1, tx_y+1, tx_z-1));
				}
			}
		}else if (rx_geometry == "D"){ // Just for five cells
			rx_id.push_back(getId(tx_x+destination, tx_y+1, tx_z));
			rx_id.push_back(getId(tx_x+destination, tx_y, tx_z));
			rx_id.push_back(getId(tx_x+destination, tx_y+1, tx_z+1));
			rx_id.push_back(getId(tx_x+destination, tx_y+2, tx_z));
			rx_id.push_back(getId(tx_x+destination, tx_y+1, tx_z-1));
		}
	}

};

class Gillespie {
	Network *tecido;

public:
	Gillespie(Network *rede) {
		tecido = rede;

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				cout << "(" << j << "," << i << ") " << tecido->getId(i, j, trunc(DIM_Z / 2)) << " " << tecido->get(i, j, trunc(DIM_Z / 2), "C") << " ";
			}
			cout << endl;
		}

		//cout << sigma0(0) << endl;
	}

	// Reaction 1
	double sigma0(int id) {
		return tecido->get(id, "vin");
	}

	// Reaction 2
	double sigma1(int id) {
		return 4 * tecido->get(id, "VM3") * ((pow(tecido->get(id, "kcaaa"), tecido->get(id, "n")) * pow(tecido->get(id, "C"), tecido->get(id, "n")) ) / ((pow(tecido->get(id, "C"), tecido->get(id, "n")) + pow(tecido->get(id, "kcaaa"), tecido->get(id, "n"))) * (pow(tecido->get(id, "C"), tecido->get(id, "n")) + pow(tecido->get(id, "kcai"), tecido->get(id, "n"))) ) * (pow(tecido->get(id, "IP3"), tecido->get(id, "m")) / (pow(tecido->get(id, "kip3"), tecido->get(id, "m")) + pow(tecido->get(id, "IP3"), tecido->get(id, "m")) ) )) * (tecido->get(id, "ER") - tecido->get(id, "C"));
	}

	// Reaction 3
	double sigma2(int id) {
		return tecido->get(id, "VM2") * ( pow(tecido->get(id, "C"), 2) / (pow(tecido->get(id, "K2"), 2) + pow(tecido->get(id, "C"), 2)) );
	}

	// Reaction 4
	double kf_Ea(int id) {
		return tecido->get(id, "kf") * tecido->get(id, "ER");
	}

	// Reaction 5
	double kf_Ca(int id) {
		return tecido->get(id, "kf") * tecido->get(id, "C");
	}

	// Reaction 6
	double ko_Ca(int id) {
		return tecido->get(id, "ko") * tecido->get(id, "C");
	}

	// Reaction 7
	double sigma3(int id) {
		return tecido->get(id, "vp") * (pow(tecido->get(id, "C"), 2) / (pow(tecido->get(id, "C"), 2) + pow(tecido->get(id, "kp"), 2)) );
	}

	// Reaction 8
	double kd_Ia(int id) {
		return tecido->get(id, "kdeg") * tecido->get(id, "IP3");
	}

	// Reaction 9
	vector<double> diffusions(int id) {
		int nConnections = tecido->numberConnections();
		vector<double> diffusions(nConnections * 3);

		vector<int> connections(nConnections);
		connections = tecido->getConnections(id);
		double value;

		for (int i = 0; i < connections.size(); i++){
			for (int gj = 0; gj < 3; gj++) {
				if (connections[i] != -1)
					value = diffusionEquation(id, connections[i], gj);
				else
					value = 0;

				diffusions[(3 * i) + gj] = value;
			}
		}

		return diffusions;
	};
	double diffusionEquation(int id1, int id2, int gap_junction) {
		double vol_cell = (4 / 3) * (PI * pow((diameter_cell / 2), 3));
		double diff;

		if (tecido->get(id1, "C") <= tecido->get(id2, "C"))
			return 0;

		if (gap_junction == 0) {
			diff = (tecido->get(id1, "D") / vol_cell) * (tecido->get(id1, "C") - tecido->get(id2, "C")) * tecido->get(id1, "phh");
		} else if (gap_junction == 1) {
			diff = (tecido->get(id1, "D") / vol_cell) * (tecido->get(id1, "C") - tecido->get(id2, "C")) * tecido->get(id1, "phl");
		} else {
			diff = (tecido->get(id1, "D") / vol_cell) * (tecido->get(id1, "C") - tecido->get(id2, "C")) * tecido->get(id1, "plh");
		}

		return diff;
	}

	// Reaction 10
	vector<double> Na_i_diffusions(int id) {
		int nConnections = tecido->numberConnections();
		vector<double> diffusions(nConnections * 3);

		vector<int> connections(nConnections);
		connections = tecido->getConnections(id);
		double value;

		for (int i = 0; i < connections.size(); i++){
			for (int gj = 0; gj < 3; gj++) {
				if (connections[i] != -1)
					value = Na_i_diffusionEquation(id, connections[i], gj);
				else
					value = 0;

				diffusions[(3 * i) + gj] = value;
			}
		}

		return diffusions;
	};
	double Na_i_diffusionEquation(int id1, int id2, int gap_junction) {
		double vol_cell = (4 / 3) * (PI * pow((diameter_cell / 2), 3));
		double diff;

		if (tecido->get(id1, "Na_i") <= tecido->get(id2, "Na_i"))
			return 0;

		if (gap_junction == 0) {
			diff = (tecido->get(id1, "NaD") / vol_cell) * (tecido->get(id1, "Na_i") - tecido->get(id2, "Na_i")) * tecido->get(id1, "phh");
		} else if (gap_junction == 1) {
			diff = (tecido->get(id1, "NaD") / vol_cell) * (tecido->get(id1, "Na_i") - tecido->get(id2, "Na_i")) * tecido->get(id1, "phl");
		} else {
			diff = (tecido->get(id1, "NaD") / vol_cell) * (tecido->get(id1, "Na_i") - tecido->get(id2, "Na_i")) * tecido->get(id1, "plh");
		}

		return diff;
	}

	// Reaction 11
	vector<double> Ca_o_diffusions(int id) {
		int nConnections = tecido->numberConnections();
		vector<double> diffusions(nConnections * 3);

		vector<int> connections(nConnections);
		connections = tecido->getConnections(id);
		double value;

		for (int i = 0; i < connections.size(); i++){
			for (int gj = 0; gj < 3; gj++) {
				if (connections[i] != -1)
					value = Ca_o_diffusionEquation(id, connections[i], gj);
				else
					value = 0;

				diffusions[(3 * i) + gj] = value;
			}
		}

		return diffusions;
	};
	double Ca_o_diffusionEquation(int id1, int id2, int gap_junction) {
		double vol_cell = (4 / 3) * (PI * pow((diameter_cell / 2), 3));
		double diff;

		if (tecido->get(id1, "C_o") <= tecido->get(id2, "C_o"))
			return 0;

		if (gap_junction == 0) {
			diff = (tecido->get(id1, "D") / vol_cell) * (tecido->get(id1, "C_o") - tecido->get(id2, "C_o")) * tecido->get(id1, "phh");
		} else if (gap_junction == 1) {
			diff = (tecido->get(id1, "D") / vol_cell) * (tecido->get(id1, "C_o") - tecido->get(id2, "C_o")) * tecido->get(id1, "phl");
		} else {
			diff = (tecido->get(id1, "D") / vol_cell) * (tecido->get(id1, "C_o") - tecido->get(id2, "C_o")) * tecido->get(id1, "plh");
		}

		return diff;
	}

	// Reaction 12
	vector<double> Na_o_diffusions(int id) {
		int nConnections = tecido->numberConnections();
		vector<double> diffusions(nConnections * 3);

		vector<int> connections(nConnections);
		connections = tecido->getConnections(id);
		double value;

		for (int i = 0; i < connections.size(); i++){
			for (int gj = 0; gj < 3; gj++) {
				if (connections[i] != -1)
					value = Na_o_diffusionEquation(id, connections[i], gj);
				else
					value = 0;

				diffusions[(3 * i) + gj] = value;
			}
		}

		return diffusions;
	};
	double Na_o_diffusionEquation(int id1, int id2, int gap_junction) {
		double vol_cell = (4 / 3) * (PI * pow((diameter_cell / 2), 3));
		double diff;

		if (tecido->get(id1, "Na_o") <= tecido->get(id2, "Na_o"))
			return 0;

		if (gap_junction == 0) {
			diff = (tecido->get(id1, "NaD") / vol_cell) * (tecido->get(id1, "Na_o") - tecido->get(id2, "Na_o")) * tecido->get(id1, "phh");
		} else if (gap_junction == 1) {
			diff = (tecido->get(id1, "NaD") / vol_cell) * (tecido->get(id1, "Na_o") - tecido->get(id2, "Na_o")) * tecido->get(id1, "phl");
		} else {
			diff = (tecido->get(id1, "NaD") / vol_cell) * (tecido->get(id1, "Na_o") - tecido->get(id2, "Na_o")) * tecido->get(id1, "plh");
		}

		return diff;
	}

	//Reaction 13 - NCX
	double Na_Ca_exchanger(int id, int *NCX_mode) {
		double dCa_dt, J_NCX, Allo, Delta_E_num, Delta_E_denom;
		//cout << "Entrei no NCX " << endl;
		//(tecido->get(id, "Erev") >= tecido->get(id, "Vm")) &&
		// Modo reverso: Potencial de reversão >= Potencial de repouso; DDC de Na (intra-extra) >= Concentracao de Ref  
		if (tecido->get(id, "Na_i") >= 15000 && tecido->get(id, "C") <= 0.5) { // Kirischuk 2012 [uM]
			// Eq. (6) from Koenigsberger, 2004
			*NCX_mode = 1;
			J_NCX = tecido->get(id, "G_NCX")*tecido->get(id, "C")*(tecido->get(id, "Vm") - tecido->get(id, "V_NCX"))/(tecido->get(id, "C") + tecido->get(id, "C_NCX"));
			/*
			Allo = 1/(1 + pow(tecido->get(id, "K_mCaAct")/tecido->get(id, "C"), tecido->get(id, "n_Hill")));
			Delta_E_num = tecido->get(id, "J_max")*( pow(tecido->get(id, "Na_i"),3)*tecido->get(id, "C_o")*exp(tecido->get(id, "eta")*(tecido->get(id, "V_NCX")/1000)*F/(R*T)) - pow(tecido->get(id, "Na_o"),3)*tecido->get(id, "C")*exp((tecido->get(id, "eta")-1)*(tecido->get(id, "V_NCX")/1000)*F/(R*T)));
			Delta_E_denom = ( tecido->get(id, "K_mCao")*pow(tecido->get(id, "Na_i"),3) + pow(tecido->get(id, "K_mNao"),3)*tecido->get(id, "C") + pow(tecido->get(id, "K_mNai"),3)*tecido->get(id, "C_o")*(1+tecido->get(id, "C")/tecido->get(id, "K_mCai")) + tecido->get(id, "K_mCai")*pow(tecido->get(id, "Na_o"),3)*(1+pow(tecido->get(id, "Na_i")/tecido->get(id, "K_mCai"),3)) + pow(tecido->get(id, "Na_i"),3)*tecido->get(id, "C_o") + pow(tecido->get(id, "Na_o"),3)*tecido->get(id, "C") )*(1 + tecido->get(id, "k_sat")*exp((tecido->get(id, "eta")-1)*(tecido->get(id, "V_NCX")/1000)*F/(R*T)));
			
			J_NCX = Allo*Delta_E_num/Delta_E_denom;
			//dCa_dt = J_NCX - (tecido->get(id, "C_rest")-tecido->get(id, "C"))/tecido->get(id, "tau_Ca");

			*///cout << setprecision(16) << "Reverse Mode -- J_NCX = " << J_NCX << endl;
			//tecido->get(id, "vin");
			//cout << "Reverse Mode = " << tecido->get(id, "Vm") << endl;
			return J_NCX;

		}
		//(tecido->get(id, "Erev") < tecido->get(id, "Vm")) && 
		// Modo direto: Potencial de reversão < Potencial de repouso; DDC de Na (intra-extra) < Concentracao de Ref
		else if (tecido->get(id, "Na_i") < 15000 && tecido->get(id, "C") > 0.5) { // Kirischuk 2012 [uM]
			// Eq. (6) from Koenigsberger, 2004
			*NCX_mode = 2;
			J_NCX = tecido->get(id, "G_NCX")*tecido->get(id, "C")*(tecido->get(id, "Vm") - tecido->get(id, "V_NCX"))/(tecido->get(id, "C") + tecido->get(id, "C_NCX"));
			/*
			Allo = 1/(1 + pow(tecido->get(id, "K_mCaAct")/tecido->get(id, "C"), tecido->get(id, "n_Hill")));
			Delta_E_num = tecido->get(id, "J_max")*( pow(tecido->get(id, "Na_i"),3)*tecido->get(id, "C_o")*exp(tecido->get(id, "eta")*(tecido->get(id, "V_NCX")/1000)*F/(R*T)) - pow(tecido->get(id, "Na_o"),3)*tecido->get(id, "C")*exp((tecido->get(id, "eta")-1)*(tecido->get(id, "V_NCX")/1000)*F/(R*T)));
			Delta_E_denom = ( tecido->get(id, "K_mCao")*pow(tecido->get(id, "Na_i"),3) + pow(tecido->get(id, "K_mNao"),3)*tecido->get(id, "C") + pow(tecido->get(id, "K_mNai"),3)*tecido->get(id, "C_o")*(1+tecido->get(id, "C")/tecido->get(id, "K_mCai")) + tecido->get(id, "K_mCai")*pow(tecido->get(id, "Na_o"),3)*(1+pow(tecido->get(id, "Na_i")/tecido->get(id, "K_mCai"),3)) + pow(tecido->get(id, "Na_i"),3)*tecido->get(id, "C_o") + pow(tecido->get(id, "Na_o"),3)*tecido->get(id, "C") )*(1 + tecido->get(id, "k_sat")*exp((tecido->get(id, "eta")-1)*(tecido->get(id, "V_NCX")/1000)*F/(R*T)));
			
			J_NCX = Allo*Delta_E_num/Delta_E_denom;
			//dCa_dt = J_NCX + (tecido->get(id, "C_rest")-tecido->get(id, "C"))/tecido->get(id, "tau_Ca");

			*///cout << setprecision(16) << "Forward Mode -- J_NCX = " << J_NCX << endl;
			//tecido->get(id, "ko") * tecido->get(id, "C");
			//cout << "Forward Mode = " << tecido->get(id, "Vm") << endl;
			return J_NCX;
		}
		else{
			*NCX_mode = 0;
			return 0;
		}

	}

	//Reaction 14 - Na+/K+ Pump

	// REACTIONS - FIM

	vector<double> calciumReactions() {
		int nConnections = tecido->numberConnections();
		int NC = DIM_X * DIM_Y * DIM_Z; // Total number of cells
		int num_reactions = 9; // Number of reactions (7 Intracellular + 1 Intercelular)
		double max_reaction = 0, reaction_choice, alfa_0 = 0, reaction_value;
		vector<double> retorno(5);
		vector<double> diffusion(nConnections * 3);
		double reactions[DIM_Y][DIM_X][DIM_Z][num_reactions + nConnections * 3 - 1];

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {

				for (int k = 0; k < DIM_Z; k++) {
					// << Begin Reactions
					for (int r = 0; r < num_reactions; r++) {
						if (r == 0) {
							reaction_value = sigma0(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 1) {
							reaction_value = sigma1(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 2) {
							reaction_value = sigma2(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 3) {
							reaction_value = kf_Ea(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 4) {
							reaction_value = kf_Ca(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 5) {
							reaction_value = ko_Ca(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 6) {
							reaction_value = sigma3(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 7) {
							reaction_value = kd_Ia(tecido->getId(i, j, k));
							reactions[j][i][k][r] = reaction_value;
						} else if (r == 8) {
							diffusion = diffusions(tecido->getId(i, j, k));
							for (int d = 0; d < diffusion.size(); d++) {
								reaction_value = diffusion[d];
								reactions[j][i][k][r + d] = reaction_value;

								/*if (reaction_value >= max_reaction) {
									max_reaction = reaction_value;
									cell[0] = j; cell[1] = i; cell[2] = k;
									if (d >= 0 && d <= 2)
										reaction_choice = 8;
									else if (d >= 3 && d <= 5)
										reaction_choice = 9;
									else if (d >= 6 && d <= 8)
										reaction_choice = 10;
									else if (d >= 9 && d <= 11)
										reaction_choice = 11;
									else if (d >= 12 && d <= 14)
										reaction_choice = 12;
									else
										reaction_choice = 13;
								}*/

								alfa_0 += reaction_value;
							}
						}
						
						if (r != 8) {
							/*if (reaction_value >= max_reaction) {
								max_reaction = reaction_value;
								cell[0] = j; cell[1] = i; cell[2] = k;
								reaction_choice = r;
							}*/
							
							alfa_0 += reaction_value;
						}
						//cout << reaction_value << endl;
					}
					// End Reactions >>
				}
			}
		}


		// for (int i = 0; i < DIM_X; i++) {
		// 	for (int j = 0; j < DIM_Y; j++) {
		// 		for (int k = 0; k < DIM_Z; k++) {
		// 			cout << tecido->getId(i, j, k) << ": ";
		// 			for (int r = 0; r < num_reactions; r++) {
		// 				cout << reactions[i][j][k][r] << " ";
		// 			}
		// 			cout << endl;
		// 		}
		// 	}
		// }


		// Gerando dois números aleatórios: r1 e r2
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator (seed);
		std::uniform_real_distribution<double> distribution (0.0,1.0);

		double r1 = distribution(generator);
		double r2 = distribution(generator);

		// Calculando o tempo tau
		//cout << setprecision(8) << "alfa0 " << alfa_0 << endl;
		double tau = (1 / alfa_0) * log(1 / r1);

		// Definindo a reação que será executada
		double sum_upper = 0, sum_down = 0;
		bool flag = false;

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				for (int k = 0; k < DIM_Z; k++) {
					for (int r = 0; r < (num_reactions + nConnections * 3 - 1); r++) {
						sum_upper += reactions[j][i][k][r];

						if (sum_upper >= alfa_0 * r2) {
							//cout << i << " " << j << " " << k << " " << r << endl;

							flag = false;
							sum_down = 0;
							for (int x = 0; x < DIM_X && flag == false; x++) {
								for (int y = 0; y < DIM_Y && flag == false; y++) {
									for (int z = 0; z < DIM_Z && flag == false; z++) {
										for (int n = 0; n < (num_reactions + nConnections * 3 - 1) && flag == false; n++) {
											if (x == i && y == j && z == k && n == r) {
												//cout << x << " " << y << " " << z << " " << n << endl;
												//cout << sum_down << " " << alfa_0 * r2 << " " << sum_upper << endl;
												flag = true;
											} else {
												sum_down += reactions[y][x][z][n];
											}
										}
									}
								}
							}

							if (sum_down < alfa_0 * r2) {
								//cout << sum_upper << " " << alfa_0 * r2 << " " << sum_down << endl;
								retorno[0] = r + 1; // [1, 26]
								retorno[1] = i;
								retorno[2] = j;
								retorno[3] = k;
								retorno[4] = tau;
								//cout << setprecision(5) << tau << endl;

								return retorno;
							}
						}
					}
				}
			}
		}
	}

	vector<double> sodiumInter() {
		int nConnections = tecido->numberConnections();
		int NC = DIM_X * DIM_Y * DIM_Z; // Total number of cells
		int num_reactions = 1; 
		double max_reaction = 0, reaction_choice, alfa_0 = 0, reaction_value;
		vector<double> retorno(5);
		vector<double> Na_i_diffusion(nConnections * 3);
		double reactions[DIM_Y][DIM_X][DIM_Z][num_reactions + (nConnections * 3 - 1)]; //0 a 43

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {

				for (int k = 0; k < DIM_Z; k++) {
					// << Begin Reactions
					Na_i_diffusion = Na_i_diffusions(tecido->getId(i, j, k));
					for (int d = 0; d < Na_i_diffusion.size(); d++) {
						reaction_value = Na_i_diffusion[d];
						reactions[j][i][k][d] = reaction_value;
						//printf("alfa_sodium = %f \n", reaction_value);
						alfa_0 += reaction_value;
					}
					// End Reactions >>
					
				}
			}
		}


		// for (int i = 0; i < DIM_X; i++) {
		// 	for (int j = 0; j < DIM_Y; j++) {
		// 		for (int k = 0; k < DIM_Z; k++) {
		// 			cout << tecido->getId(i, j, k) << ": ";
		// 			for (int r = 0; r < num_reactions; r++) {
		// 				cout << reactions[i][j][k][r] << " ";
		// 			}
		// 			cout << endl;
		// 		}
		// 	}
		// }


		// Gerando dois números aleatórios: r1 e r2
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator (seed);
		std::uniform_real_distribution<double> distribution (0.0,1.0);

		double r1 = distribution(generator);
		double r2 = distribution(generator);

		// Calculando o tempo tau
		double tau = (1 / alfa_0) * log(1 / r1);

		// Definindo a reação que será executada
		double sum_upper = 0, sum_down = 0;
		bool flag = false;

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				for (int k = 0; k < DIM_Z; k++) {
					for (int r = 0; r < (num_reactions + (nConnections * 3 - 1)); r++) {
						// cout << "R: " << r <<endl;
						sum_upper += reactions[j][i][k][r];

						if (sum_upper >= alfa_0 * r2) {
							//cout << i << " " << j << " " << k << " " << r << endl;

							flag = false;
							sum_down = 0;
							for (int x = 0; x < DIM_X && flag == false; x++) {
								for (int y = 0; y < DIM_Y && flag == false; y++) {
									for (int z = 0; z < DIM_Z && flag == false; z++) {
										for (int n = 0; n < (num_reactions +  (nConnections * 3 - 1)) && flag == false; n++) { //0 a 43 -> tem que ir
											// cout << "N: " << r <<endl;

											if (x == i && y == j && z == k && n == r) {
												// cout << x << " " << y << " " << z << " " << n << endl;
												// cout << sum_down << " " << alfa_0 * r2 << " " << sum_upper << endl;
												flag = true;
											} else {
												sum_down += reactions[y][x][z][n];
											}
										}
									}
								}
							}

							if (sum_down < alfa_0 * r2) {
								//cout << sum_upper << " " << alfa_0 * r2 << " " << sum_down << endl;
								retorno[0] = r + 1; // [1, 26]
								retorno[1] = i;
								retorno[2] = j;
								retorno[3] = k;
								retorno[4] = tau;
								// cout << "Reaction: " << r+1 <<endl;
								// cout << fixed << setprecision(10) << "Tau: " << tau <<endl;
								return retorno;
							}
						}
					}
				}
			}
		}
	}

	vector<double> sodiumExtra() {
		int nConnections = tecido->numberConnections();
		int NC = DIM_X * DIM_Y * DIM_Z; // Total number of cells
		int num_reactions = 1; 
		double max_reaction = 0, reaction_choice, alfa_0 = 0, reaction_value;
		vector<double> retorno(5);
		vector<double> Na_o_diffusion(nConnections * 3);
		double reactions[DIM_Y][DIM_X][DIM_Z][num_reactions + (nConnections * 3 - 1)]; //0 a 43

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {

				for (int k = 0; k < DIM_Z; k++) {
					// << Begin Reactions
					Na_o_diffusion = Na_o_diffusions(tecido->getId(i, j, k));
					for (int d = 0; d < Na_o_diffusion.size(); d++) {
						reaction_value = Na_o_diffusion[d];
						reactions[j][i][k][d] = reaction_value;

						alfa_0 += reaction_value;
					}
					// End Reactions >>

				}
			}
		}


		// for (int i = 0; i < DIM_X; i++) {
		// 	for (int j = 0; j < DIM_Y; j++) {
		// 		for (int k = 0; k < DIM_Z; k++) {
		// 			cout << tecido->getId(i, j, k) << ": ";
		// 			for (int r = 0; r < num_reactions; r++) {
		// 				cout << reactions[i][j][k][r] << " ";
		// 			}
		// 			cout << endl;
		// 		}
		// 	}
		// }


		// Gerando dois números aleatórios: r1 e r2
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator (seed);
		std::uniform_real_distribution<double> distribution (0.0,1.0);

		double r1 = distribution(generator);
		double r2 = distribution(generator);

		// Calculando o tempo tau
		double tau = (1 / alfa_0) * log(1 / r1);

		// Definindo a reação que será executada
		double sum_upper = 0, sum_down = 0;
		bool flag = false;

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				for (int k = 0; k < DIM_Z; k++) {
					for (int r = 0; r < (num_reactions + (nConnections * 3 - 1)); r++) {
						// cout << "R: " << r <<endl;
						sum_upper += reactions[j][i][k][r];

						if (sum_upper >= alfa_0 * r2) {
							//cout << i << " " << j << " " << k << " " << r << endl;

							flag = false;
							sum_down = 0;
							for (int x = 0; x < DIM_X && flag == false; x++) {
								for (int y = 0; y < DIM_Y && flag == false; y++) {
									for (int z = 0; z < DIM_Z && flag == false; z++) {
										for (int n = 0; n < (num_reactions +  (nConnections * 3 - 1)) && flag == false; n++) { //0 a 43 -> tem que ir
											// cout << "N: " << r <<endl;

											if (x == i && y == j && z == k && n == r) {
												// cout << x << " " << y << " " << z << " " << n << endl;
												// cout << sum_down << " " << alfa_0 * r2 << " " << sum_upper << endl;
												flag = true;
											} else {
												sum_down += reactions[y][x][z][n];
											}
										}
									}
								}
							}

							if (sum_down < alfa_0 * r2) {
								//cout << sum_upper << " " << alfa_0 * r2 << " " << sum_down << endl;
								retorno[0] = r + 1; // [1, 26]
								retorno[1] = i;
								retorno[2] = j;
								retorno[3] = k;
								retorno[4] = tau;
								// cout << "Reaction: " << r+1 <<endl;
								// cout << fixed << setprecision(10) << "Tau: " << tau <<endl;
								return retorno;
							}
						}
					}
				}
			}
		}
	}

	vector<double> calciumExtra() {
		int nConnections = tecido->numberConnections();
		int NC = DIM_X * DIM_Y * DIM_Z; // Total number of cells
		int num_reactions = 1; 
		double max_reaction = 0, reaction_choice, alfa_0 = 0, reaction_value;
		vector<double> retorno(5);
		vector<double> Ca_o_diffusion(nConnections * 3);
		double reactions[DIM_Y][DIM_X][DIM_Z][num_reactions + (nConnections * 3 - 1)]; //0 a 43

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				for (int k = 0; k < DIM_Z; k++) {
					// << Begin Reactions

					Ca_o_diffusion = Ca_o_diffusions(tecido->getId(i, j, k));
					for (int d = 0; d < Ca_o_diffusion.size(); d++) {
						reaction_value = Ca_o_diffusion[d];
						reactions[j][i][k][d] = reaction_value;

						alfa_0 += reaction_value;
					}
					// End Reactions >>

				}
			}
		}


		// for (int i = 0; i < DIM_X; i++) {
		// 	for (int j = 0; j < DIM_Y; j++) {
		// 		for (int k = 0; k < DIM_Z; k++) {
		// 			cout << tecido->getId(i, j, k) << ": ";
		// 			for (int r = 0; r < num_reactions; r++) {
		// 				cout << reactions[i][j][k][r] << " ";
		// 			}
		// 			cout << endl;
		// 		}
		// 	}
		// }


		// Gerando dois números aleatórios: r1 e r2
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator (seed);
		std::uniform_real_distribution<double> distribution (0.0,1.0);

		double r1 = distribution(generator);
		double r2 = distribution(generator);

		// Calculando o tempo tau
		double tau = (1 / alfa_0) * log(1 / r1);

		// Definindo a reação que será executada
		double sum_upper = 0, sum_down = 0;
		bool flag = false;

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				for (int k = 0; k < DIM_Z; k++) {
					for (int r = 0; r < (num_reactions + (nConnections * 3 - 1)); r++) {
						// cout << "R: " << r <<endl;
						sum_upper += reactions[j][i][k][r];

						if (sum_upper >= alfa_0 * r2) {
							//cout << i << " " << j << " " << k << " " << r << endl;

							flag = false;
							sum_down = 0;
							for (int x = 0; x < DIM_X && flag == false; x++) {
								for (int y = 0; y < DIM_Y && flag == false; y++) {
									for (int z = 0; z < DIM_Z && flag == false; z++) {
										for (int n = 0; n < (num_reactions +  (nConnections * 3 - 1)) && flag == false; n++) { //0 a 43 -> tem que ir
											// cout << "N: " << r <<endl;

											if (x == i && y == j && z == k && n == r) {
												// cout << x << " " << y << " " << z << " " << n << endl;
												// cout << sum_down << " " << alfa_0 * r2 << " " << sum_upper << endl;
												flag = true;
											} else {
												sum_down += reactions[y][x][z][n];
											}
										}
									}
								}
							}

							if (sum_down < alfa_0 * r2) {
								//cout << sum_upper << " " << alfa_0 * r2 << " " << sum_down << endl;
								retorno[0] = r + 1; // [1, 26]
								retorno[1] = i;
								retorno[2] = j;
								retorno[3] = k;
								retorno[4] = tau;
								// cout << "Reaction: " << r+1 <<endl;
								// cout << fixed << setprecision(10) << "Tau: " << tau <<endl;
								return retorno;
							}
						}
					}
				}
			}
		}
	}

	vector<double> NCX_reaction(vector<int> &NCX_mode_vector) {
		int NCX_mode, num_reactions = 40; // number of reactions
		int NCX_mode_array[DIM_Y][DIM_X][DIM_Z][num_reactions];
		double max_reaction = 0, reaction_choice, alfa_0 = 0, reaction_value;
		vector<double> retorno(5);
		double reactions[DIM_Y][DIM_X][DIM_Z];

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				for (int k = 0; k < DIM_Z; k++) {
					for (int r = 0; r < num_reactions; r++) {
						tecido->update_parameters(i, j , k);
						
						reaction_value = Na_Ca_exchanger(tecido->getId(i, j, k), &NCX_mode);
						reactions[j][i][k] += reaction_value;
						NCX_mode_array[j][i][k][r] = NCX_mode;
						
						alfa_0 += reaction_value;
						// cout << reaction_value << endl;
					}
				}
			}
		}

		// for (int i = 0; i < DIM_X; i++) {
		// 	for (int j = 0; j < DIM_Y; j++) {
		// 		for (int k = 0; k < DIM_Z; k++) {
		// 			cout << tecido->getId(i, j, k) << ": ";
		// 				cout << reactions[i][j][k][r] << " ";
		// 			cout << endl;
		// 		}
		// 	}
		// }

		// Gerando dois números aleatórios: r1 e r2
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator (seed);
		std::uniform_real_distribution<double> distribution (0.0,1.0);

		double r1 = distribution(generator);
		double r2 = distribution(generator);

		// Calculando o tempo tau
		//cout << setprecision(8) << "alfa0 " << alfa_0 << endl;
		double tau = (1 / alfa_0) * log(1 / r1);

		// Definindo a reação que será executada
		double sum_upper = 0, sum_down = 0;
		bool flag = false;

		for (int i = 0; i < DIM_X; i++) {
			for (int j = 0; j < DIM_Y; j++) {
				for (int k = 0; k < DIM_Z; k++) {
					sum_upper += reactions[j][i][k];

					if (sum_upper >= alfa_0 * r2) {
						//cout << i << " " << j << " " << k << " " << r << endl;

						flag = false;
						sum_down = 0;
						for (int x = 0; x < DIM_X && flag == false; x++) {
							for (int y = 0; y < DIM_Y && flag == false; y++) {
								for (int z = 0; z < DIM_Z && flag == false; z++) {
									if (x == i && y == j && z == k) {
										//cout << x << " " << y << " " << z << " " << n << endl;
										//cout << sum_down << " " << alfa_0 * r2 << " " << sum_upper << endl;
										flag = true;
									} else {
										sum_down += reactions[y][x][z];
									}
								}
							}
						}
						if (sum_down < alfa_0 * r2) {
							//cout << sum_upper << " " << alfa_0 * r2 << " " << sum_down << endl;
							for (int n=0; n<num_reactions; n++){
								NCX_mode_vector.push_back(NCX_mode_array[j][i][k][n]);
							}
							retorno[0] = 0;
							retorno[1] = i;
							retorno[2] = j;
							retorno[3] = k;
							retorno[4] = tau;
							//cout << tau << endl;

							return retorno;
						}
					}
				}
			}
		}
	}
};

void simulation(int destination, double frequency, string topology, double time_slot, ofstream& file_results, int Rx_number, string rx_geometry) {
	Network tecido;
	int tx_x = trunc(DIM_X / 15), tx_x_2 = trunc(DIM_X / 15); // Rx depends on that position
	int tx_y = trunc(DIM_Y / 5), tx_y_2 = trunc(3 * DIM_Y / 5); // Rx depends on that position
	int tx_z = trunc(DIM_Z / 3), tx_z_2 = trunc(DIM_Z / 3); // Rx depends on that position
	int radius = 0;

	// SETTING THE Rx GEOMETRY
	tecido.rx_geometry(rx_geometry, Rx_number, destination, tx_x, tx_y, tx_z);

	// SETTING THE TOPOLOGY OF THE TISSUE
	if (topology == "RD"){
		tecido.regularDegree();
		radius = 1;
	}
	else if (topology == "LR2"){ 
		tecido.linkRadius(2);
		radius = 2;
	}
	else if (topology == "LR3"){
		tecido.linkRadius(3);
		radius = 3;
	}

	// SETTING THE VALUES
	tecido.set(tx_x, tx_y, tx_z, "C", 0.5*SCALE); // Tx 1
	tecido.set(tx_x_2, tx_y_2, tx_z_2, "C", 0.5*SCALE); // Tx 2
	// tecido.set(tx_x, tx_y, tx_z, "Na_i", 20000*SCALE); // Quais as referências?
	// tecido.set(tx_x, tx_y, tx_z, "C_o", 20*SCALE); // Quais as referências?
	// tecido.set(tx_x, tx_y, tx_z, "Na_o", 1000*SCALE); // Quais as referências?

	// Print tissue
	tecido.printTissue();

	// Opening file containing the data of calcium concentration that will be plotted
	// ofstream exportfile;
	// exportfile.open("temp/data.txt");
	// tecido.writeFileHeader(exportfile);
	// tecido.printTissuef(exportfile, 0); // Writes the tissue's initial state to the file
	
	// Opening file for sotring concentration data at time simulation
	// ofstream cdatafile;
	// cdatafile.open("temp/cdata.csv");
	// cdatafile << "time,c_in,\n"; 

	// Inicializando o Algoritmo de Gillespie
	Gillespie gillespie(&tecido);

	int nConnections = tecido.numberConnections(), num_reactions = 40;
	cout << "Connections per cell: " << nConnections << endl;
	vector<double> choice(5);
	vector<int> connections(nConnections), qtd_reactions(9 + QTD_DIFFUSIONS * (nConnections * 3)), NCX_mode_vector, Rx_states, Tx_states;
	double simulation_time = 200, current_time = 0, current_time_calcium = 0, current_time_sodium_inter = 0, current_time_NCX = 0;
	double tau_max = 100000, tau_calcium=0, tau_sodium_inter=0, tau_NCX=0, E_signal = 0, E_noise = 0, c_in=0, c_out=0, current_time_mod_demod = 0;
	int reaction, int_time = 0, x_c, y_c, z_c, bit, bit2, bit3, time_slots_number = destination-1;
	bool diffusion_error = false, tau_flag = false;
	vector<int>::iterator first;
	vector<double> C_tx, C_rx;

	while (simulation_time > current_time) {
		
		// tau_flag = false;
		// current_time_calcium = 0; //current_time_sodium_inter = 0; current_time_sodium_extra = 0; current_time_calcium_extra = 0; current_time_NCX = 0;
		// tau_max = 100000;
		
		// do
		// {	
			// Updating time and concentrations
			if (trunc(current_time) != int_time) {
				int_time = trunc(current_time);

				cout << "Time: " << int_time << endl;

				// Print Reactions
				cout << "Reactions: ";
				for (int i = 0; i < qtd_reactions.size(); i++) {
					cout << qtd_reactions[i] << " ";
				}
				cout << endl;

				// End Print Reactions

				// Print Tissue
				tecido.printTissue();
				cout << endl;	
				// End Print Tissue

				// Write Calcium concentration of tissue on file
				// tecido.printTissuef(exportfile, int_time);
				// End Print Tissue

				// Update Gap Junctions
				if (int_time < 50) {
					for (int i = 0; i < DIM_X; i++) {
						for (int j = 0; j < DIM_Y; j++) {
							for (int k = 0; k < DIM_Z; k++) {
								tecido.set(i, j, k, "phl", phl[int_time]);
								tecido.set(i, j, k, "plh", plh[int_time]);
								tecido.set(i, j, k, "phh", phh[int_time]);
							}
						}
					}
				}

				// Calcium and sodium oscillations
				if (int_time % ((int) (1 / frequency)) == 0){
					tecido.accumulate(tx_x, tx_y, tx_z, "C", 0.5*SCALE);
					tecido.accumulate(tx_x_2, tx_y_2, tx_z_2, "C", 0.5*SCALE);
					//tecido.accumulate(tx_x, tx_y, tx_z, "Na_i", 20000*SCALE);
				}
			}

			// CHOICE OF THE REACTIONS
			diffusion_error = false;

			// INTRACELULAR CALCIUM REACTIONS
			// if (current_time_calcium<tau_max) {	
				choice = gillespie.calciumReactions();
				reaction = choice[0];
				tau_max = choice[4];
				current_time_mod_demod += tau_max*1000;
				
				x_c = choice[1];
				y_c = choice[2];
				z_c = choice[3];

				//cout << setprecision(5) << "Tau_calcium = " << tau_calcium << endl;
				qtd_reactions[reaction - 1]++;

				if (reaction == 1) {
					tecido.accumulate(choice[1], choice[2], choice[3], "C", ALPHA);
				} else if (reaction == 2) {
					tecido.accumulate(choice[1], choice[2], choice[3], "C", ALPHA);
					tecido.accumulate(choice[1], choice[2], choice[3], "ER", -ALPHA);

					if (tecido.get(choice[1], choice[2], choice[3], "ER") < 0) {
						tecido.changeSignal(choice[1], choice[2], choice[3], "ER");
					}
				} else if (reaction == 3) {
					tecido.accumulate(choice[1], choice[2], choice[3], "C", -ALPHA);
					tecido.accumulate(choice[1], choice[2], choice[3], "ER", ALPHA);

					if (tecido.get(choice[1], choice[2], choice[3], "C") < 0) {
						tecido.changeSignal(choice[1], choice[2], choice[3], "C");
					}
				} else if (reaction == 4) {
					tecido.accumulate(choice[1], choice[2], choice[3], "C", ALPHA);
					tecido.accumulate(choice[1], choice[2], choice[3], "ER", -ALPHA);

					if (tecido.get(choice[1], choice[2], choice[3], "ER") < 0) {
						tecido.changeSignal(choice[1], choice[2], choice[3], "ER");
					}
				} else if (reaction == 5) {
					tecido.accumulate(choice[1], choice[2], choice[3], "C", -ALPHA);
					tecido.accumulate(choice[1], choice[2], choice[3], "ER", ALPHA);

					if (tecido.get(choice[1], choice[2], choice[3], "C") < 0) {
						tecido.changeSignal(choice[1], choice[2], choice[3], "C");
					}
				} else if (reaction == 6) {
					tecido.accumulate(choice[1], choice[2], choice[3], "C", -ALPHA);

					if (tecido.get(choice[1], choice[2], choice[3], "C") < 0) {
						tecido.changeSignal(choice[1], choice[2], choice[3], "C");
					}
				} else if (reaction == 7) {
					tecido.accumulate(choice[1], choice[2], choice[3], "IP3", ALPHA);
				} else if (reaction == 8) {
					tecido.accumulate(choice[1], choice[2], choice[3], "IP3", -ALPHA);

					if (tecido.get(choice[1], choice[2], choice[3], "IP3") < 0) {
						tecido.changeSignal(choice[1], choice[2], choice[3], "IP3");
					}
				}

				/* DIFFUSION REACTIONS OF CALCIUM */
				else {
					for (int conn = 0; conn < nConnections; conn++) {
						if (reaction >= 9 + (conn * 3) && reaction <= 11 + (conn * 3)) {
							connections = tecido.getConnections(tecido.getId(choice[1], choice[2], choice[3]));

							if (connections[conn] != -1 && tecido.get(choice[1], choice[2], choice[3], "C") > tecido.get(connections[conn], "C")) {
								tecido.accumulate(choice[1], choice[2], choice[3], "C", -ALPHA);
								tecido.accumulate(connections[conn], "C", ALPHA);

								// Modulation and Demodulation - Begin
								if ( (x_c == tx_x && y_c == tx_y && z_c == tx_z) || (x_c == tx_x_2 && y_c == tx_y_2 && z_c == tx_z_2) ){
									tecido.accumulate(x_c, y_c, z_c, "C_variation", ALPHA);
								}
								first = tecido.rx_id.begin();
								for (; first != tecido.rx_id.end(); first++){
									if (connections[conn] == *first) tecido.accumulate(*first, "C_variation", ALPHA);
								}
								// Modulation and Demodulation - End

							} else {
								//cout << tecido.getId(choice[1], choice[2], choice[3]) << " " << connections[conn] << endl;
								diffusion_error = true;
							}
						}
						// E_signal = E_signal + pow(0.1,2);
					}
				}
						// /* DIFFUSION ERROR => NO INCREMENT TIME */
				if (diffusion_error){
					current_time -= (tau_max * 1000);
					current_time_mod_demod -= (tau_max * 1000);
				}
				// if (reaction<9) E_signal = E_signal + pow(0.1,2);

				// current_time_calcium += tau_calcium;					
			// }

			// Modulation and Demodulation - Begin
			if (current_time_mod_demod >= time_slot){
				// Transmitter
				bit = tecido.mod_demod(tx_x, tx_y, tx_z); // Tx1
				bit2 = tecido.mod_demod(tx_x_2, tx_y_2, tx_z_2); // Tx2
				if (bit == 1 || bit2 == 1){
					Tx_states.push_back(1);
				}else{
					Tx_states.push_back(bit);
				}
				// cout << "Tx Bit: " << bit << endl;
				tecido.set(tx_x, tx_y, tx_z, "C_variation", 0);
				tecido.set(tx_x_2, tx_y_2, tx_z_2, "C_variation", 0);

				// Receiver
				first = tecido.rx_id.begin();
				bit2 = 0;
				bit3 = 0;
				for (; first != tecido.rx_id.end(); first++){ // Rx
					bit = tecido.mod_demod(*first);
					if (bit == 1) bit2++;
					else if (bit == 0) bit3++;
					// cout << "Rx Bit: " << bit << endl;
					tecido.set(*first, "C_variation", 0);
				}
				if (bit2 != 0) Rx_states.push_back(1);
				else if (bit3 > trunc(Rx_number/2)) Rx_states.push_back(0);
				else Rx_states.push_back(-1);

				current_time_mod_demod = 0;
			}
			// Modulation and Demodulation - End

			// Saving concentration signal into csv file 
			// if (x_c != tx_x && y_c != tx_y && z_c != tx_z) {
			// 		c_in = tecido.get(x_c, y_c, z_c, "C");
			// 		cdatafile << int_time << "," << c_in << ",\n";
			// }else{
			// 	cdatafile << int_time << "," << c_in << ",\n";
			// }

			//INTERCELLULAR SODIUM REACTIONS
			// if (current_time_sodium_inter<tau_max) {
			// 	choice = gillespie.sodiumInter();
			// 	tau_sodium_inter = choice[4];
			// 	reaction = choice[0];
			// 	//cout << setprecision(5) << "Tau_sodium_inter = " << tau_sodium_inter << endl;
			// 	for (int conn = 0; conn < nConnections; conn++) {
			// 		if (reaction >= 1 + (conn * 3) && reaction <= 3 + (conn * 3)) {
			// 			connections = tecido.getConnections(tecido.getId(choice[1], choice[2], choice[3]));

			// 			if (connections[conn] != -1 && tecido.get(choice[1], choice[2], choice[3], "Na_i") > tecido.get(connections[conn], "Na_i")) {
			// 				tecido.accumulate(choice[1], choice[2], choice[3], "Na_i", -ALPHA);
			// 				tecido.accumulate(connections[conn], "Na_i", ALPHA);
			// 			} else {
			// 				//cout << tecido.getId(choice[1], choice[2], choice[3]) << " " << connections[conn] << endl;
			// 				diffusion_error = true;
			// 			}
			// 		}
			// 	}
			// 	current_time_sodium_inter += tau_sodium_inter;
			// 	/* DIFFUSION ERROR => NO INCREMENT TIME */
			// 	if (diffusion_error) current_time_sodium_inter -= (tau_sodium_inter);
			// }

			// NCX - Reaction
			// if (current_time_NCX < tau_max) {
			// 	choice = gillespie.NCX_reaction(NCX_mode_vector);
			// 	tau_NCX = choice[4];

			// 	for (int i = 0; i < num_reactions; i++){
					
			// 		if (NCX_mode_vector[i] == 1) { // Modo Reverso
			// 			if (!i) cout << "MODO REVERSO" << endl;
			// 			tecido.accumulate(choice[1], choice[2], choice[3], "Na_o", 3*ALPHA);
			// 			tecido.accumulate(choice[1], choice[2], choice[3], "Na_i", -3*ALPHA);
			// 			tecido.accumulate(choice[1], choice[2], choice[3], "C", ALPHA);
			// 			tecido.accumulate(choice[1], choice[2], choice[3], "C_o", -ALPHA);
			// 		} else if (NCX_mode_vector[i] == 2) { // Modo Direto
			// 			if (!i) cout << "MODO DIRETO" << endl;
			// 			tecido.accumulate(choice[1], choice[2], choice[3], "Na_i", 3*ALPHA);
			// 			tecido.accumulate(choice[1], choice[2], choice[3], "Na_o", -3*ALPHA);
			// 			tecido.accumulate(choice[1], choice[2], choice[3], "C_o", ALPHA);
			// 			tecido.accumulate(choice[1], choice[2], choice[3], "C", -ALPHA);
			// 		}
			// 		E_noise = E_noise + pow(0.3,2);
			// 	}
			// 	current_time_NCX += tau_NCX;
			// }

			//cout << setprecision(5) << "C = " << tecido.get(choice[1], choice[2], choice[3], "C") << endl;
			// if (!tau_flag){
			// 	tau_max = max(tau_calcium, max(tau_sodium_inter, tau_NCX));
			// 	tau_flag=true;
			current_time += tau_max*1000;
			// }
			// if (tau_NCX == 0 || tau_calcium == 0){
			// 	cout << "Não houve reacao!" << endl;
			// 	break;
			// }

			/* STORAGE OF CALCIUM CONCENTRATION */
			C_tx.push_back(tecido.get(tx_x,tx_y,tx_z, "C")+tecido.get(tx_x_2,tx_y_2,tx_z_2, "C"));
			
			first = tecido.rx_id.begin();
			c_out = 0;
			for (; first != tecido.rx_id.end(); first++) c_out = c_out + tecido.get(*first, "C");
			C_rx.push_back(c_out);

			//  || current_time_sodium_inter<tau_max

		// } while (current_time_calcium<tau_max || current_time_NCX<tau_max);
	}
	// cdatafile.close();

	/* ### CALCULATING GAIN ### */

	double acc_c_tx = accumulate(C_tx.begin(), C_tx.end(), 0.0);
	double acc_c_rx = accumulate(C_rx.begin(), C_rx.end(), 0.0);

	// cout << "E_signal = " << E_signal << "; E_noise = " << E_noise << endl;

	// double calc_SNR = 10 * log10(E_signal/E_noise);
	double calc_gain = 10 * log10((acc_c_rx / C_rx.size()) / ((acc_c_tx) / C_tx.size()));

	/* ### END GAIN ### */

	/* ### CALCULATING SISO CHANEL CAPACITY AND BER - BEGIN ### */

	int bit_number = 2;
	double channel_capacity = 0, px[bit_number] = {}, py[bit_number] = {}, pyx_joint[bit_number][bit_number] = {}, BER = 0;
	vector<double> I_xy;

	px[1] = tecido.conditional_accumulate(Tx_states.begin(), Tx_states.end(), 1)/Tx_states.size();
	px[0] = 1 - px[1];
	
	py[1] = tecido.conditional_accumulate(Rx_states.begin(), Rx_states.end(), 1)/Rx_states.size();
	if (py[1] == 0) py[1] = 0.00000000001;
	py[0] = 1 - py[1];

	for (int i = 0; i < Tx_states.size()-time_slots_number; i++){
		// cout << Rx_states[i+time_slots_number] << endl;

		if (Rx_states[i+time_slots_number] != -1 && Tx_states[i] != -1) pyx_joint[Rx_states[i+time_slots_number]][Tx_states[i]]++;
		if (Rx_states[i+time_slots_number] != Tx_states[i]) BER++;
	}
	pyx_joint[0][0] = Tx_states.size()-pyx_joint[1][0];
	pyx_joint[0][1] = Tx_states.size()-pyx_joint[1][1];
	BER = BER/Tx_states.size();

	for (int y = 0; y < bit_number; y++){ // Number of y1 given x0; Number of y1 given x1; Number of y0 given x0; Number of y0 given x1
		for (int x = 0; x < bit_number; x++){
			pyx_joint[y][x] = pyx_joint[y][x]/Tx_states.size();
			// cout << "Pyx_joint: " << pyx_joint[y][x] << endl;
			// cout << "Py: " << py[y] << endl;
			if (pyx_joint[y][x] == 0) pyx_joint[y][x] = 0.00000000001;

			I_xy.push_back(px[x] * pyx_joint[y][x] * log2(pyx_joint[y][x]/py[y])); // Mutual Information
		}
	}

	// first = I_xy.begin();
	// for (; first != I_xy.end(); first++){
	// 	cout << "I_xy: " << *first << endl;
	// }

	channel_capacity = *max_element(I_xy.begin(), I_xy.end());

	/* ### CALCULATING SISO CHANEL CAPACITY - END ### */

	file_results << topology << "," << destination << "," << frequency <<  "," << calc_gain << "," << channel_capacity << "," << BER << ",\n";

};

/* MAIN */
int main(){

	int simulation_number = 1, Rx_number = 5; // 0 < Rx_number < x && Rx_number <= y
	double time_slot = 0.1; //0.1,,0.5,0.8,1 s
	vector<double> frequencies{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1}; //{0.6};[Hz]
	string topology = "RD";
	string rx_geometry = "HL"; // HL = Horizontal Line; VL = Vertical Line; X = Cross (<= 5 cells); D = Diamond (== 5 cells)
	// int destination = 1;

	ofstream file_results;
	file_results.open("results/results.csv");
	file_results << "Topology,Range,Freq (Hz),Gain (dB),Channel Capacity (bits),BER,\n"; 

	for (int j = 0; j < simulation_number; j++)
	{
		for (double frequency : frequencies)
		{
			for (int destination = 1; destination < 7; destination++)
			{
				simulation(destination, frequency, topology, time_slot, file_results, Rx_number, rx_geometry);
			}
		}
	}
	file_results.close();

	return 0;
};
