#pragma once
#ifndef mes
#define mes

class Node
{
	int id;
	double x, y, t;
	int status;
public:
	Node();
	Node(double x, double y);
	int get_id() { return id; };
	double get_x() { return x; };
	double get_y() { return y; };
	double get_t() { return t; };
	double get_status() { return status; };
	void set_t(double t) { this->t = t; };

};

class Surface
{
	Node *n;
public:
	Surface();
	Surface(Node n1, Node n2);
	Node *get_n() { return n; };
};

class Element
{
	friend class Fourier;
	int nr;
	Node *ID;
	Surface *p;
	int *global_id;
	int bound_surface_number;				//ilosc powierzchni brzegowych
	int *local_surface_number;				//lokalne numery powierzchni
	double **local_H;
	double *local_P;
	int bc_counter;
public:
	Element();
	Element(double x, double y);
	Node *get_n() { return ID; };
	int *get_licz_pow() { return local_surface_number; };
	int get_lp() { return bound_surface_number; };
};

class Global_Data
{
	friend class Grid;
	friend class Fourier;
	double H, B, nH, nB, t0, t, delta_tau, tau, alpha, cw, k, ro;
	int nh, ne;
public:
	Global_Data();
	double get_tau() { return tau; };
	double get_deltaTau() { return delta_tau; };
	double get_nh() { return nh; };
	double get_nB() { return nB; };
	double get_nH() { return nH; };
	double get_t0() { return t0; };
	double get_H() { return H; };
	double get_B() { return B; };
	double get_ne() { return ne; };
};

class Grid
{
	friend class Fourier;
public:
	Node * nd;
	Element *el;
	Grid(Global_Data data);
	Node *get_nd() { return nd; };
};

class Local_Node
{
public:
	double ksi, eta;
	Local_Node(double ksi, double eta);
	Local_Node();
	double get_ksi() { return ksi; };
	double get_eta() { return eta; };
};

class Local_Surface
{
public:
	Local_Node node[2];// [2][2];
	double N[2][4];
	void set_n(Local_Node n1, Local_Node n2);
};

class Local_Element
{
	double **dN_ksi;
	double **dN_eta;
	double **N;
	Local_Surface *pl;
	Local_Node *pc;

	double N1_ksi(double eta) { return (-(1.0 / 4.0) * (1 - eta)); };
	double N1_eta(double ksi) { return (-(1.0 / 4.0) * (1 - ksi)); };

	double N2_ksi(double eta) { return ((1.0 / 4.0) * (1 - eta)); };
	double N2_eta(double ksi) { return (-(1.0 / 4.0) * (1 + ksi)); };

	double N3_ksi(double eta) { return ((1.0 / 4.0) * (1 + eta)); };
	double N3_eta(double ksi) { return ((1.0 / 4.0) * (1 + ksi)); };

	double N4_ksi(double eta) { return (-(1.0 / 4.0) * (1 + eta)); };
	double N4_eta(double ksi) { return ((1.0 / 4.0) * (1 - ksi)); };

	double N1(double ksi, double eta) { return 0.25 * (1.0 - ksi) * (1.0 - eta); };
	double N2(double ksi, double eta) { return 0.25 * (1.0 + ksi) * (1.0 - eta); };
	double N3(double ksi, double eta) { return 0.25 * (1.0 + ksi) * (1.0 + eta); };
	double N4(double ksi, double eta) { return 0.25 * (1.0 - ksi) * (1.0 + eta); };

public:
	Local_Element();
	double **get_dN_ksi() { return dN_ksi; };
	double **get_dN_eta() { return dN_eta; };
	double **get_N() { return N; };
	Local_Surface *get_pl() { return pl; };
};

class Jacobian
{
	double **J;
	double det_J;
	double **J_odw;

public:
	Jacobian();
	void set_J(Local_Element el, double *x, double *y, int pc);
	double **get_J() { return J; }
	double get_det_J() { return det_J; };
	double **get_J_odw() { return J_odw; };

};

class Fourier {
	double **global_H;
	double *global_P;
	double dNdx[4];
	double dNdy[4];
	double C;
	Local_Element el_l;
public:
	Fourier(Global_Data gd);
	void compute(Global_Data gd, Grid g);
	double **get_globalH() { return global_H; };
	double *get_globalP() { return global_P; };
};

double *Gauss(int n, double **H, double *P);
double Min(double *tab, int n);
double Max(double *tab, int n);
#endif // !1