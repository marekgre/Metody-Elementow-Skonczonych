#include "mes.h"
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

Node::Node()
{
	Global_Data global_data;
	t = global_data.get_t0();
}

Node::Node(double xx, double yy)
{
	Global_Data global_data;
	double B = global_data.get_B(); 
	double H = global_data.get_H();
	this->x = xx;
	this->y = yy;
	this->t = global_data.get_t0();

	//if (this->y == 0 || this->y >= H) status = 1;                   //  dla
	//else status = 0;                                                //  2 

	if ((this->y == 0 || this->y >= H) || (this->x == 0 || this->x >= B)) status = 1;   // dla 4
	else status = 0;                                                                    // dla 4
}

Surface::Surface()
{

}

Surface::Surface(Node n1, Node n2)
{
	n = new Node[2];
	n[0] = n1;
	n[1] = n2;
}

Element::Element()
{
	bc_counter = 0;
	local_H = new double *[4];
	local_P = new double[4];
	for (int i = 0; i < 4; i++)
		local_H[i] = new double[4];
}

Element::Element(double x, double y)
{
	local_H = new double *[4];
	local_P = new double[4];
	for (int i = 0; i < 4; i++)
		local_H[i] = new double[4];
	Global_Data global_data;
	ID = new Node[4];

	global_id = new int[4];
	int i = x;
	int j = y;

	global_id[0] = global_data.get_nH() * i + j;
	global_id[1] = global_data.get_nH() * (i + 1) + j;
	global_id[2] = global_data.get_nH() * (i + 1) + (j + 1);
	global_id[3] = global_data.get_nH() * i + (j + 1);
	
	x = x*(global_data.get_B() / (global_data.get_nB() - 1));
	y = y*(global_data.get_H() / (global_data.get_nH() - 1));

	ID[0] = Node(x, y);
	ID[1] = Node(x + (global_data.get_B() / (global_data.get_nB() - 1)), y);
	ID[2] = Node(x + (global_data.get_B() / (global_data.get_nB() - 1)), y + (global_data.get_H() / (global_data.get_nH() - 1)));
	ID[3] = Node(x, y + (global_data.get_H() / (global_data.get_nH() - 1)));

	p = new Surface[4];
	p[0] = Surface(ID[3], ID[0]);
	p[1] = Surface(ID[0], ID[1]);
	p[2] = Surface(ID[1], ID[2]);
	p[3] = Surface(ID[2], ID[3]);

	bound_surface_number = 0;
	for (int i = 0; i < 4; i++)
	{
		if (p[i].get_n()[0].get_status() == 1 && p[i].get_n()[1].get_status() == 1)
			bound_surface_number++;
	}
	local_surface_number = new int[bound_surface_number];
	int counter = 0;
	for (int i = 0; i < 4; i++)
	{
		if (p[i].get_n()[0].get_status() == 1 && p[i].get_n()[1].get_status() == 1)
			local_surface_number[counter++] = i;
	}
}

Grid::Grid(Global_Data data)
{
	Global_Data global_data;
	el = new Element[global_data.get_ne()];
	int k = 0;
	for (int i = 0; i < (global_data.get_nB()-1); i++)
	{
		for (int j = 0; j < (global_data.get_nH() - 1); j++)
			el[k++] = Element(i, j);
	}
	nd = new Node[global_data.get_nh()];
	k = 0;
	for (int i = 0; i < global_data.get_nB(); i++)
	{
		for (int j = 0; j < global_data.get_nH(); j++)
			nd[k++] = Node(i*(global_data.get_B() / (global_data.get_nB() - 1)), j*(global_data.get_H() / (global_data.get_nH() - 1)));
	}
}

Global_Data::Global_Data()
{
	ifstream file;
	file.open("File.txt");
	file >> H;
	file >> B;
	file >> nH;
	file >> nB;
	nh = nB * nH;
	ne = (nB - 1) * (nH - 1);
	file >> t0;
	file >> t;
	file >> delta_tau;
	file >> tau;
	file >> alpha;
	file >> cw;
	file >> k;
	file >> ro;
	file.close();
	double asr = (k / (cw*ro));
	delta_tau = (pow((B / nB), 2) / (0.5*asr));
}

Jacobian::Jacobian()
{
	J = new double *[2];
	J_odw = new double *[2];
	for (int i = 0; i < 2; i++)
	{
		J[i] = new double[2];
		J_odw[i] = new double[2];
	}
}

void Jacobian::set_J(Local_Element le, double *x, double *y, int pc)
{
	J[0][0] = (le.get_dN_ksi()[pc][0] * x[0] + le.get_dN_ksi()[pc][1] * x[1] + le.get_dN_ksi()[pc][2] * x[2] + le.get_dN_ksi()[pc][3] * x[3]);
	J[1][0] = (le.get_dN_ksi()[pc][0] * y[0] + le.get_dN_ksi()[pc][1] * y[1] + le.get_dN_ksi()[pc][2] * y[2] + le.get_dN_ksi()[pc][3] * y[3]);
	J[0][1] = (le.get_dN_eta()[pc][0] * x[0] + le.get_dN_eta()[pc][1] * x[1] + le.get_dN_eta()[pc][2] * x[2] + le.get_dN_eta()[pc][3] * x[3]);
	J[1][1] = (le.get_dN_eta()[pc][0] * y[0] + le.get_dN_eta()[pc][1] * y[1] + le.get_dN_eta()[pc][2] * y[2] + le.get_dN_eta()[pc][3] * y[3]);

	J_odw[0][0] = J[1][1];
	J_odw[1][0] = -1.0*J[1][0];
	J_odw[0][1] = -1.0*J[0][1];
	J_odw[1][1] = J[0][0];

	det_J = J[0][0] * J[1][1] - J[1][0] * J[0][1];
}

Local_Node::Local_Node(double ksi, double eta)
{
	this->ksi = ksi;
	this->eta = eta;
}

Local_Node::Local_Node()
{

}

void Local_Surface::set_n(Local_Node n1, Local_Node n2)
{
	node[0] = n1;
	node[1] = n2;
}

Local_Element::Local_Element()
{
	dN_ksi = new double*[4];
	dN_eta = new double*[4];
	N = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		dN_ksi[i] = new double[4];
		dN_eta[i] = new double[4];
		N[i] = new double[4];
	}
	pc = new Local_Node[4];
	pc[0].ksi = -1.0 / sqrt(3.0);
	pc[0].eta = -1.0 / sqrt(3.0);

	pc[1].ksi = 1.0 / sqrt(3.0);
	pc[1].eta = -1.0 / sqrt(3.0); 
	
	pc[2].ksi = -1.0 / sqrt(3.0);
	pc[2].eta = 1.0 / sqrt(3.0);
	
	pc[3].ksi = 1.0 / sqrt(3.0);
	pc[3].eta = 1.0 / sqrt(3.0);

	for (int i = 0; i < 4; i++)
	{
		dN_ksi[i][0] = N1_ksi(pc[i].get_ksi());
		dN_ksi[i][1] = N2_ksi(pc[i].get_ksi());
		dN_ksi[i][2] = N3_ksi(pc[i].get_ksi());
		dN_ksi[i][3] = N4_ksi(pc[i].get_ksi());

		dN_eta[i][0] = N1_eta(pc[i].get_eta());
		dN_eta[i][1] = N2_eta(pc[i].get_eta());
		dN_eta[i][2] = N3_eta(pc[i].get_eta());
		dN_eta[i][3] = N4_eta(pc[i].get_eta());

		N[i][0] = N1(pc[i].get_ksi(), pc[i].get_eta());
		N[i][1] = N2(pc[i].get_ksi(), pc[i].get_eta());
		N[i][2] = N3(pc[i].get_ksi(), pc[i].get_eta());
		N[i][3] = N4(pc[i].get_ksi(), pc[i].get_eta());
	}
	pl = new Local_Surface[4];
	pl[0].set_n(Local_Node(-1.0, 1.0 / sqrt(3.0)),Local_Node(-1.0, -1.0 / sqrt(3.0)));
	pl[1].set_n(Local_Node(-1.0 / sqrt(3.0), -1.0), Local_Node(1.0 / sqrt(3.0), -1.0));
	pl[2].set_n(Local_Node(1.0, -1.0 / sqrt(3.0)), Local_Node(1.0, 1.0 / sqrt(3.0)));
	pl[3].set_n(Local_Node(1.0 / sqrt(3.0), 1.0), Local_Node(-1.0 / sqrt(3.0), 1.0));

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			pl[i].N[j][0] = N1(pl[i].node[j].get_ksi(), pl[i].node[j].get_eta());
			pl[i].N[j][1] = N2(pl[i].node[j].get_ksi(), pl[i].node[j].get_eta());
			pl[i].N[j][2] = N3(pl[i].node[j].get_ksi(), pl[i].node[j].get_eta());
			pl[i].N[j][3] = N4(pl[i].node[j].get_ksi(), pl[i].node[j].get_eta());
		}
	}
}

Fourier::Fourier(Global_Data gd)
{
	global_H = new double *[gd.nh];
	global_P = new double[gd.nh];
	for (int i = 0; i < gd.nh; i++)
		global_H[i] = new double[gd.nh];
}

void Fourier::compute(Global_Data global_data, Grid grid)
{	
	double temp_interpol;
	int id;
	double x[4], y[4], t_el[4], detJ = 0;

	for (int i = 0; i < global_data.nh; i++)
	{
		for (int j = 0; j < global_data.nh; j++)
			global_H[i][j] = 0;
		global_P[i] = 0;
	}
	for (int i = 0; i < global_data.ne; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
				grid.el[i].local_H[j][k] = 0;
			grid.el[i].local_P[j] = 0;
		}
		for (int j = 0; j < 4; j++)
		{
			id = grid.el[i].global_id[j];
			x[j] = grid.el[i].ID[j].get_x();
			y[j] = grid.el[i].ID[j].get_y();
			t_el[j] = grid.get_nd()[id].get_t();
		}
		for (int pc = 0; pc < 4; pc++)
		{
			Jacobian J;
			J.set_J(el_l, x, y, pc);
			temp_interpol = 0;
			for (int j = 0; j < 4; j++)
			{
				dNdx[j] = 1.0 / J.get_det_J()*(J.get_J_odw()[0][0] * el_l.get_dN_ksi()[pc][j]+ J.get_J_odw()[0][1] * el_l.get_dN_eta()[pc][j]);
				dNdy[j] = 1.0 / J.get_det_J()*(J.get_J_odw()[1][0] * el_l.get_dN_ksi()[pc][j] + J.get_J_odw()[1][1] * el_l.get_dN_eta()[pc][j]);
				temp_interpol += t_el[j] * el_l.get_N()[pc][j];
			}
			detJ = abs(J.get_det_J());
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					C = global_data.cw*global_data.ro*el_l.get_N()[pc][j] * el_l.get_N()[pc][k] * detJ;
					grid.el[i].local_H[j][k] += global_data.k*(dNdx[j]*dNdx[k]+dNdy[j]*dNdy[k])*detJ+C/global_data.get_deltaTau();					
					grid.el[i].local_P[j] += C / global_data.get_deltaTau()*temp_interpol;
				}
			}
		}
		for (int j = 0; j < grid.el[i].get_lp(); j++)
		{
			id = grid.el[i].get_licz_pow()[j];
			switch (id) {
				case 0:
					detJ = sqrt(pow(grid.el[i].ID[3].get_x() - grid.el[i].ID[0].get_x(), 2) + pow(grid.el[i].ID[3].get_y() - grid.el[i].ID[0].get_y(), 2))/2.0;
					break;
				case 1:
					detJ = sqrt(pow(grid.el[i].ID[0].get_x() - grid.el[i].ID[1].get_x(), 2) + pow(grid.el[i].ID[0].get_y() - grid.el[i].ID[1].get_y(), 2))/2.0;
					break;
				case 2:
					detJ = sqrt(pow(grid.el[i].ID[1].get_x() - grid.el[i].ID[2].get_x(), 2) + pow(grid.el[i].ID[1].get_y() - grid.el[i].ID[2].get_y(), 2))/2.0;
					break;
				case 3:
					detJ = sqrt(pow(grid.el[i].ID[2].get_x() - grid.el[i].ID[3].get_x(), 2) + pow(grid.el[i].ID[2].get_y() - grid.el[i].ID[3].get_y(), 2))/2.0;
					break;
			}
			for (int k = 0; k < 2; k++)
			{
				for (int m = 0; m < 4; m++)
				{
					for (int n = 0; n < 4; n++)
						grid.el[i].local_H[m][n] += global_data.alpha * el_l.get_pl()[id].N[k][m] * el_l.get_pl()[id].N[k][n] * detJ;
					grid.el[i].local_P[m] += global_data.alpha* global_data.t *detJ * el_l.get_pl()[id].N[k][m];
				}
			}
		}
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
				global_H[grid.el[i].global_id[j]][grid.el[i].global_id[k]] += grid.el[i].local_H[j][k];
			global_P[grid.el[i].global_id[j]] += grid.el[i].local_P[j];
		}
	}
}

double *Gauss(int n, double **H, double *P) {
	double m, s, e = pow(10, -12);
	double *result = new double[n];
	for (int i = 0; i < n; i++)
		result[i] = 0.0;
	double **tab = new double*[n];
	for (int i = 0; i < n; i++)
		tab[i] = new double[n + 1];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			tab[j][i] = H[j][i];
	}
	for (int i = 0; i < n; i++)
		tab[i][n] = P[i];
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i+1; j < n; j++)
		{
			if (abs(tab[i][i])<e)
			{
				cout << ("dzielnik rowny 0") << endl;
				break;
			}
			m = -1.0*tab[j][i] / tab[i][i];
			for (int k = i+1; k <= n; k++)
			{
				tab[j][k] += m*tab[i][k];
			}
		}
	}
	for (int i = n-1; i >= 0; i--)
	{
		s = tab[i][n];
		for (int j = n - 1; j >= i+1; j--) {
			s -= tab[i][j] * result[j];	
		}
		if (abs(tab[i][i]) < e)
		{
			cout << "dzielnik rowny 0" << endl;
			break;
		}
		result[i] = s / tab[i][i];
	}
	return result;
}

double Min(double *tab, int n)
{
	double min = tab[0];
	for (int i = 0; i < n; i++)
	{
		if (tab[i] < min)
			min = tab[i];
	}
	return min;
}

double Max(double *tab, int n)
{
	double max = tab[0];
	for (int i = 0; i < n; i++)
	{
		if (tab[i] > max)
			max = tab[i];
	}
	return max;
}