/**************************************************************************
 * Calling code for the simulation of a neuronal network with
 * channelrhodopsin-2 channels
 **************************************************************************
 *
 * Copyright (C) Jannik Luboeinski 2015-2020
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

// Mandatory compiler options: -std=c++0x
// Mandatory linker options: -lgsl -lgslcblas [-static]
// Compiler option to define Linux platform (default): -D LINUX
// Compiler option to define Windows platform: -D WINDOWS

//#define SEEK_I_CONST // if defined, I_const will be varied and the I_const value that leads to a mean frequency of 5 Hz (in the absence of light) will be determined
//#define GAUSS_PROFILES // uses a Gaussian connection profile instead of a uniform probability distribution for network connections

#include "NetworkSimulation.cpp"

vector<double> stimplot {2.5, 5.0, 7.5, 10.}; // the irradiances for which firing rate map plots shall be created
vector<int> rasterplot {0, 1799, 3599, 3600, 4499}; // the neurons which shall be included in the spike raster plot

int Nl = 60; // number of neurons in one row of the excitatory population (total number: Nl^2)
int Nl_inh = 30; // number of neurons in one row of the inhibitory population (total number: Nl_inh^2)
const double dt = 0.1; // ms, duration of one time step
double dE = 2.5; // mW/mm^2, size of one irradiance step
double t_max = 20000.0; // ms, total duration of the simulation
double t_offset = 0.0; // ms, offset for the beginning of the stimulus
double E_max = 20.; // mW/mm^2, maximum irradiance value
double t_pulse = 4.0; // ms, duration of one light pulse
double raster_bin = 10.0; // ms, time bin for spike number and raster plot

/* default, for Nl=60, Nl_inh=30 (similar in Amit&Brunel 1997a) */
double J_ee = 0.110; // pC, coupling strength for excitatory inputs to excitatory neurons
double J_ei = 0.190; // pC, coupling strength for excitatory inputs to inhibitory neurons
double J_ie = 0.340; // pC, coupling strength for inhibitory inputs to excitatory neurons
double J_ii = 0.540; // pC, coupling strength for inhibitory inputs to inhibitory neurons

// WARNING: cal_table is restricted to certain pc values!
double pc_step = 0.001; // connection probability, step size
double pc_start = 0.005; // connection probability, initial value
double pc_end = 0.015; // connection probability, final value

// WARNING: cal_table is restricted to certain tau_syn values!
double tau_syn_step = 1.00; // ms, synaptic time constant, step size
double tau_syn_start = 5.0; // ms, synaptic time constant, initial value
double tau_syn_end = 5.0; // ms, synaptic time constant, final value

double frequency_step = 2.0; // Hz, stimulus frequency, step size
double frequency_start = 50.0; // Hz, stimulus frequency, initial value
double frequency_end = 50.0; // Hz, stimulus frequency, final value

// WARNING: cal_table is restricted to certain J values!
double J_step = 0.2; // factor to rescale the synaptic coupling strengths, step size
double J_start = 1.0; // factor to rescale the synaptic coupling strengths, initial value
double J_end = 1.0; // factor to rescale the synaptic coupling strengths, final value

double sigma_light = 8; // standard deviation of spatial Gaussian light stimulus

// constant currents in case of an uncoupled, non-illuminated network (in nA)
// sigma_WN       I_const
// 0.003          0.974262
// 0.010          0.914576
// 0.020          0.825474
// 0.040          0.642673
// 0.095          0.116017

double I_const = 0.914576; // nA, constant current applied to LIF neurons (initially for uncoupled neuron, may be re-set using the tables below
                           // if SEEK_I_CONST is not set)
bool use_cal_table = true; // specifies whether or not the calibration table for I_const is used (else, the I_const value upon declaration is used)
                           // for certain parameter settings, this adjusts the current to achieve a mean firing rate of 5 Hz without stimulation

const double cal_table [3][15][10] = { // transfer matrix for I_const in mixed network
                                       // -> values that were computed for the respective networks to reach 5 Hz on average
                                       //    firing rates are slightly above 5 Hz (computation a bit different as described, thus values are bit higher)
                                       //    (values for Ne=3600, Ni=900, sigma_WN=0.010 nA)
// tau_syn = 1 ms
/*   J	=	0.2		 0.4		  0.6			0.8		 1.0		  1.2			1.4		 1.6		  1.8			2.0  */
			 {{0., 		 0., 		  0., 		0., 		 0.873978, 0., 		0., 		 0., 		  0., 		0.}, 		 // pc = 0.001
			  {0., 		 0., 		  0., 		0., 		 0.873224, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.002
			  {0., 		 0., 		  0., 		0., 		 0.872537, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.003
			  {0., 		 0., 		  0., 		0., 		 0.871913, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.004
			  {0., 		 0., 		  0., 		0., 		 0.871349, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.005
			  {0., 		 0., 		  0., 		0., 		 0.870840, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.006
			  {0., 		 0., 		  0., 		0., 		 0.870384, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.007
			  {0., 		 0., 		  0., 		0., 		 0.869978, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.008
			  {0., 		 0., 		  0., 		0., 		 0.869617, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.009
			  {0.873359, 0.872415, 0.871418, 0.870383, 0.869300, 0.868138, 0.866852, 0.865390, 0.863707, 0.861771},// 		0.010
			  {0., 		 0., 		  0., 		0., 		 0.869023, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.011
			  {0., 		 0., 		  0., 		0., 		 0.868783, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.012
			  {0., 		 0., 		  0., 		0., 		 0.868578, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.013
			  {0., 		 0., 		  0., 		0., 		 0.868406, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.014
			  {0., 		 0., 		  0., 		0., 		 0.868264, 0., 		0., 		 0., 		  0., 		0.}},		 // 		0.015 */

// tau_syn = 5 ms
/*   J	=	0.2		 0.4		  0.6			0.8		 1.0		  1.2			1.4		 1.6		  1.8			2.0  */
			 {{0., 		 0., 		  0., 		0., 		 0.914488, 0., 		0., 		 0., 		  0., 		0.}, 		 // pc = 0.001
			  {0., 		 0., 		  0., 		0., 		 0.913883, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.002
			  {0., 		 0., 		  0., 		0., 		 0.913334, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.003
			  {0., 		 0., 		  0., 		0., 		 0.912838, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.004
			  {0., 		 0., 		  0., 		0., 		 0.912391, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.005
			  {0., 		 0., 		  0., 		0., 		 0.911990, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.006
			  {0., 		 0., 		  0., 		0., 		 0.911632, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.007
			  {0., 		 0., 		  0., 		0., 		 0.911314, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.008
			  {0., 		 0., 		  0., 		0., 		 0.911032, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.009
			  {0.913777, 0.912985, 0.912221, 0.911491, 0.910786, 0.910081, 0.909343, 0.908538, 0.907634, 0.906604},// 		0.010
			  {0., 		 0., 		  0., 		0., 		 0.910571, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.011
			  {0., 		 0., 		  0., 		0., 		 0.910387, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.012
			  {0., 		 0., 		  0., 		0., 		 0.910230, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.013
			  {0., 		 0., 		  0., 		0., 		 0.910099, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.014
			  {0., 		 0., 		  0., 		0., 		 0.909992, 0., 		0., 		 0., 		  0., 		0.}},		 // 		0.015 */

// tau_syn = 10 ms
/*   J	=	0.2		 0.4		  0.6			0.8		 1.0		  1.2			1.4		 1.6		  1.8			2.0  */
			 {{0., 		 0., 		  0., 		0., 		 0.944842, 0., 		0., 		 0., 		  0., 		0.}, 		 // pc = 0.001
			  {0., 		 0., 		  0., 		0., 		 0.944348, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.002
			  {0., 		 0., 		  0., 		0., 		 0.943900, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.003
			  {0., 		 0., 		  0., 		0., 		 0.943496, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.004
			  {0., 		 0., 		  0., 		0., 		 0.943133, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.005
			  {0., 		 0., 		  0., 		0., 		 0.942807, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.006
			  {0., 		 0., 		  0., 		0., 		 0.942517, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.007
			  {0., 		 0., 		  0., 		0., 		 0.942259, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.008
			  {0., 		 0., 		  0., 		0., 		 0.942031, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.009
			  {0.944064, 0.943384, 0.942789, 0.942277, 0.941832, 0.941432, 0.941051, 0.940660, 0.940235, 0.939756},// 		0.010
			  {0., 		 0., 		  0., 		0., 		 0.941659, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.011
			  {0., 		 0., 		  0., 		0., 		 0.941510, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.012
			  {0., 		 0., 		  0., 		0., 		 0.941383, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.013
			  {0., 		 0., 		  0., 		0., 		 0.941278, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.014
			  {0., 		 0., 		  0., 		0., 		 0.941193, 0., 		0., 		 0., 		  0., 		0.}},		 // 		0.015 */

};

const double cal_table_exc [3][15][10] = { // transfer matrix for I_const in purely excitatory network
                                           // -> values that were computed for the respective networks to reach 5 Hz on average
                                           //    (values for Ne=3600, Ni=0, sigma_WN=0.010 nA)
// tau_syn = 1 ms <==== not updated!!
/*   J	=	0.2		 0.4		  0.6			0.8		 1.0		  1.2			1.4		 1.6		  1.8			2.0  */
			 {{0., 		 0., 		  0., 		0., 		 0.873978, 0., 		0., 		 0., 		  0., 		0.},		 // pc = 0.001
			  {0., 		 0., 		  0., 		0., 		 0.873224, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.002
			  {0., 		 0., 		  0., 		0., 		 0.872537, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.003
			  {0., 		 0., 		  0., 		0., 		 0.871913, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.004
			  {0., 		 0., 		  0., 		0., 		 0.871349, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.005
			  {0., 		 0., 		  0., 		0., 		 0.870840, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.006
			  {0., 		 0., 		  0., 		0., 		 0.870384, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.007
			  {0., 		 0., 		  0., 		0., 		 0.869978, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.008
			  {0., 		 0., 		  0., 		0., 		 0.869617, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.009
			  {0.873359, 0.872415, 0.871418, 0.870383, 0.869300, 0.868138, 0.866852, 0.865390, 0.863707, 0.861771},// 		0.010
			  {0., 		 0., 		  0., 		0., 		 0.869023, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.011
			  {0., 		 0., 		  0., 		0., 		 0.868783, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.012
			  {0., 		 0., 		  0., 		0., 		 0.868578, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.013
			  {0., 		 0., 		  0., 		0., 		 0.868406, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.014
			  {0., 		 0., 		  0., 		0., 		 0.868264, 0., 		0., 		 0., 		  0., 		0.}},		 // 		0.015 */

// tau_syn = 5 ms
/*   J	=	0.2		 0.4		  0.6			0.8		 1.0		  1.2			1.4		 1.6		  1.8			2.0  */
			 {{0.,	    0.,  	  0.,  		0.,  		 0.908598, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // pc = 0.001
			  {0.,  		 0.,  	  0.,  		0.,  		 0.906541, 0.,  		0.,  		 0.,  	  0.,  		0.}, 		 // 		0.002
			  {0.,  		 0.,  	  0.,  		0.,  		 0.904357, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // 		0.003
			  {0.,  	    0.,  	  0.,  		0.,  		 0.902363, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // 		0.004
			  {0.,  		 0.,  	  0.,  		0.,  		 0.900081, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // 		0.005
			  {0.,  		 0.,  	  0.,  		0.,  		 0.898079, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // 		0.006
			  {0.,  		 0.,  	  0.,  		0.,  		 0.895898, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // 		0.007
			  {0.,  		 0.,  	  0.,  		0.,  		 0.894024, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // 		0.008
			  {0.,  		 0.,  	  0.,  		0.,  		 0.891900, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // 		0.009
			  {0.906640, 0.902484, 0.898382,	0.894104, 0.889633, 0.885228,	0.880764, 0.875952, 0.870570,	0.866064},// 		0.010
			  {0.,  		 0.,  	  0.,  		0.,  		 0.887764, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // 		0.011
			  {0.,  		 0.,  	  0.,  		0.,  		 0.885433, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // 		0.012
			  {0.,  		 0.,  	  0.,  		0.,  		 0.883325, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // 		0.013
			  {0.,  		 0.,  	  0.,  		0.,  		 0.881259, 0.,  		0.,  		 0.,  	  0.,  		0.},		 // 		0.014
			  {0.,  		 0.,  	  0.,  		0.,  		 0.878886, 0.,  		0.,  		 0.,  	  0.,  		0.}},		 // 		0.015 */

// tau_syn = 10 ms <==== not updated!!
/*   J	=	0.2		 0.4		  0.6			0.8		 1.0		  1.2			1.4		 1.6		  1.8			2.0  */
			 {{0., 		 0., 		  0., 		0., 		 0.944842, 0., 		0., 		 0., 		  0., 		0.}, 		 // pc = 0.001
			  {0., 		 0., 		  0., 		0., 		 0.944348, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.002
			  {0., 		 0., 		  0., 		0., 		 0.943900, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.003
			  {0., 		 0., 		  0., 		0., 		 0.943496, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.004
			  {0., 		 0., 		  0., 		0., 		 0.943133, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.005
			  {0., 		 0., 		  0., 		0., 		 0.942807, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.006
			  {0., 		 0., 		  0., 		0., 		 0.942517, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.007
			  {0., 		 0., 		  0., 		0., 		 0.942259, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.008
			  {0., 		 0., 		  0., 		0., 		 0.942031, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.009
			  {0.944064, 0.943384, 0.942789, 0.942277, 0.941832, 0.941432, 0.941051, 0.940660, 0.940235, 0.939756},// 		0.010
			  {0., 		 0., 		  0., 		0., 		 0.941659, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.011
			  {0., 		 0., 		  0., 		0., 		 0.941510, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.012
			  {0., 		 0., 		  0., 		0., 		 0.941383, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.013
			  {0., 		 0., 		  0., 		0., 		 0.941278, 0., 		0., 		 0., 		  0., 		0.}, 		 // 		0.014
			  {0., 		 0., 		  0., 		0., 		 0.941193, 0., 		0., 		 0., 		  0., 		0.}},		 // 		0.015 */

};

/*** main ***/
/* - return: -1 if an error occurred */
int main(int argc, char** argv)
{
	string path = getCurrentDir() + string("/Network_") + dateStr("", true); // path to working directory
#ifdef SEEK_I_CONST
	path += "_seek";
#endif
	double frequency; // current frequency value
	double pc; // current connection prob. value
	double tau_syn; // current synaptic time constant
	double J; // current scale factor for synaptic coupling strength

#ifdef WINDOWS
	system(""); // Windows needs this to start using colors in the console...
#endif
	cout << getSeparator('=') << endl
	     << "\x1b[36mSimulation Tool for Neural Networks with Light-Sensitive Channels\x1b[0m" << endl
	     << getSeparator('=') << endl
	     << "Code and documentation: https://jlubo.net/nn-lightchannels-sim" << endl
	     << "Contact: Jannik Luboeinski <jannik.lubo[at]gmx.de>" << endl
	     << getSeparator('-') << endl;

	// Change parameters that specify values for frequency, pc, J, I_const, Nl, Nl_inh and tau_syn
	if (argc > 1) // process comand-line arguments if given
	{
		for(int i=1; i<argc; i++)
		{
			char* pt;
			double read = 0.0;
			if( (pt = strstr(argv[i], "=")) != NULL )
			{
				pt++;
				read = atof(pt);
			}

			if ( strstr(argv[i], "-frequency=") == argv[i] || strstr(argv[i], "-f=") == argv[i] )
			{
				frequency_start = read;
				frequency_end = read;
			}
			else if (strstr(argv[i], "-f_start=") == argv[i])
				frequency_start = read;
			else if (strstr(argv[i], "-f_end=") == argv[i])
				frequency_end = read;
			else if (strstr(argv[i], "-pc=") == argv[i])
			{
				pc_start = read;
				pc_end = read;
			}
			else if (strstr(argv[i], "-pc_start=") == argv[i])
				pc_start = read;
			else if (strstr(argv[i], "-pc_end=") == argv[i])
				pc_end = read;
			else if (strstr(argv[i], "-tau_syn=") == argv[i])
			{
				tau_syn_start = read;
				tau_syn_end = read;
			}
			else if (strstr(argv[i], "-tau_syn_start=") == argv[i] || strstr(argv[i], "-tsstart=") == argv[i])
				tau_syn_start = read;
			else if (strstr(argv[i], "-tau_syn_end=") == argv[i] || strstr(argv[i], "-tsend=") == argv[i])
				tau_syn_end = read;
			else if (strstr(argv[i], "-J=") == argv[i])
			{
				J_start = read;
				J_end = read;
			}
			else if (strstr(argv[i], "-J_start=") == argv[i])
				J_start = read;
			else if (strstr(argv[i], "-J_end=") == argv[i])
				J_end = read;
			else if (strstr(argv[i], "-J_ee=") == argv[i])
				J_ee = read;
			else if (strstr(argv[i], "-I_const=") == argv[i])
			{
				I_const = read;
				use_cal_table = false;
			}
			else if (strstr(argv[i], "-Nl=") == argv[i])
				Nl = read;
			else if (strstr(argv[i], "-Nl_inh=") == argv[i])
				Nl_inh = read;
			else if (strstr(argv[i], "-sigma_light=") == argv[i])
				sigma_light = read;
		}
	}
	else // demand arguments from user
	{
		cout << "\x1b[33mPlease enter your desired parameter values.\x1b[0m If you do not specify anything, the default value is used." << endl << endl;

		cout << "Minimum light intensity (in mW/mm^2, fixed): 0" << endl;
		cout << "Maximum light intensity (in mW/mm^2, default is " << E_max << " mW/mm^2): ";
		readInp(E_max);

		cout << "Light intensity step size (in mW/mm^2, default is " << dE << " mW/mm^2): ";
		readInp(dE, 0., E_max);

		cout << "Duration of one light pulse (in ms, default is " << t_pulse << " ms): ";
		readInp(t_pulse);

		cout << "Stimulus frequency to start with (in Hz, default is " << frequency_start << " Hz): ";
		readInp(frequency_start);

		cout << "Stimulus frequency to end with (in Hz, default is " << frequency_end << " Hz): ";
		readInp(frequency_end);

		if (frequency_end-frequency_start > EPSILON)
		{
			cout << "Step size for the stimulus frequency (in Hz, default is " << frequency_step << " Hz): ";
			readInp(frequency_step, 0., frequency_end-frequency_start);
		}

		cout << "Prob. of connectivity to start with (default is " << pc_start << "): ";
		readInp(pc_start, 0., 1.);

		cout << "Prob. of connectivity to end with (default is " << pc_end << "): ";
		readInp(pc_end, 0., 1.);

		if (pc_end-pc_start > EPSILON)
		{
			cout << "Step size for the prob. of connectivity (default is " << pc_step << "): ";
			readInp(pc_step, 0., pc_end-pc_start);
		}

		cout << "Coupling strength scaling to start with (default is " << J_start << "): ";
		readInp(J_start);

		cout << "Coupling strength scaling to end with (default is " << J_end << "): ";
		readInp(J_end);

		if (J_end-J_start > EPSILON)
		{
			cout << "Step size for the coupling strength scaling (default is " << J_step << "): ";
			readInp(J_step, 0., J_end-J_start);
		}

		cout << "Synaptic time constant to start with (in ms, default is " << tau_syn_start << " ms): ";
		readInp(tau_syn_start);

		cout << "Synaptic time constant to end with (in ms, default is " << tau_syn_end << " ms): ";
		readInp(tau_syn_end);

		if (tau_syn_end-tau_syn_start > EPSILON)
		{
			cout << "Step size for the synaptic time constant (default is " << tau_syn_step << "): ";
			readInp(tau_syn_step, 0., tau_syn_end-tau_syn_start);
		}

		cout << "Mean external input current (in nA, default is " << I_const << " nA/automatic): ";

		if (readInp(I_const, std::numeric_limits<double>::lowest()))
		{
			use_cal_table = false;
			cout << "\x1b[33mNote:\x1b[0m Automatic adjustment of the mean input current has been disabled." << endl;
		}
		else
		{
			string inp;
			cout << "Use automatic adjustment of the mean input current? This only works with certain parameter settings (see NetworkBatch.cpp). Enter 'y' for yes or anything else for no: ";
			getline(cin, inp);
			if (inp == "y")
			{
				cout << "\x1b[33mNote:\x1b[0m Automatic adjustment of the mean input current will be used." << endl;
			}
			else
			{
				use_cal_table = false;
				cout << "\x1b[33mNote:\x1b[0m Automatic adjustment of the mean input current has been disabled." << endl;
			}
		}

		cout << "Number of neurons in one row, excitatory (default is " << Nl << "): ";
		readInp(Nl);

		cout << "Number of neurons in one row, inhibitory (default is " << Nl_inh << "): ";
		readInp(Nl_inh);

		cout << "Width of the Gaussian light stimulus in grid units (default is " << sigma_light << "): ";
		readInp(sigma_light);

		cout << "Simulation offset time (in ms, default is " << t_offset << " ms): ";
		readInp(t_offset);

		cout << "Total duration of one simulation (in ms, default is " << t_max << " ms): ";
		readInp(t_max);

		cout << getSeparator() << endl;
	}

	// Create working directory and files
	mkDir(path); // create working directory
	if (chdir(path.c_str()) == -1) { // try to change directory
		showChDirErrMessage();
		return -1;
	}

	writePalViridis(); // create Viridis color palette file
	ofstream colorPlotData(dateStr("_cplotdata.txt"));

	// Create gnuplot script files
	if (frequency_end - frequency_start > EPSILON)
	{
		createOverviewColorPlot(dE, 0.0, E_max, "mean_fr_freq", "f / Hz", frequency_start, frequency_end, frequency_step,
								  		"'%.1f'", "⟨{/Symbol n}⟩ / Hz", 5, 1, 6);
		createOverviewColorPlot(dE, dE, E_max, "gauss_sigma_freq", "f / Hz", frequency_start, frequency_end, frequency_step,
									  "'%.1f'", "{/Symbol s}_{FR}", 5, 1, 7);
		createOverviewColorPlot(dE, 0.0, E_max, "mean_fr_inh_freq", "f / Hz", frequency_start, frequency_end, frequency_step,
								  		"'%.1f'", "⟨{/Symbol n}⟩ / Hz", 5, 1, 6);
		createOverviewColorPlot(dE, dE, E_max, "gauss_sigma_inh_freq", "f / Hz", frequency_start, frequency_end, frequency_step,
									  "'%.1f'", "{/Symbol s}_{FR}", 5, 1, 7);
	}
	if (pc_end - pc_start > EPSILON)
	{
		createOverviewColorPlot(dE, 0.0, E_max, "mean_fr_pc", "p_c / %%", 100*pc_start, 100*pc_end, 100*pc_step,
									  "'%.1f'", "⟨{/Symbol n}⟩ / Hz", 5, 2, 6);
		createOverviewColorPlot(dE, dE, E_max, "gauss_sigma_pc", "p_c / %%", 100*pc_start, 100*pc_end, 100*pc_step,
									  "'%.1f'", "{/Symbol s}_{FR}", 5, 2, 7);
		createOverviewColorPlot(dE, 0.0, E_max, "mean_fr_inh_pc", "p_c / %%", 100*pc_start, 100*pc_end, 100*pc_step,
									  "'%.1f'", "⟨{/Symbol n}⟩ / Hz", 5, 2, 8);
		createOverviewColorPlot(dE, dE, E_max, "gauss_sigma_inh_pc", "p_c / %%", 100*pc_start, 100*pc_end, 100*pc_step,
									  "'%.1f'", "{/Symbol s}_{FR}", 5, 2, 9);
	}
	if (tau_syn_end - tau_syn_start > EPSILON)
	{
		createOverviewColorPlot(dE, 0.0, E_max, "mean_fr_tau_syn", "{/Symbol t}_{syn} / ms", tau_syn_start, tau_syn_end, tau_syn_step,
									  "'%.0f'", "⟨{/Symbol n}⟩ / Hz", 5, 3, 6);
		createOverviewColorPlot(dE, dE, E_max, "gauss_sigma_tau_syn", "{/Symbol t}_{syn} / ms", tau_syn_start, tau_syn_end, tau_syn_step,
									  "'%.0f'", "{/Symbol s}_{FR}", 5, 3, 7);
		createOverviewColorPlot(dE, 0.0, E_max, "mean_fr_inh_tau_syn", "{/Symbol t}_{syn} / ms", tau_syn_start, tau_syn_end, tau_syn_step,
									  "'%.0f'", "⟨{/Symbol n}⟩ / Hz", 5, 3, 8);
		createOverviewColorPlot(dE, dE, E_max, "gauss_sigma_inh_tau_syn", "{/Symbol t}_{syn} / ms", tau_syn_start, tau_syn_end, tau_syn_step,
									  "'%.0f'", "{/Symbol s}_{FR}", 5, 3, 9);
	}
	if (J_end - J_start > EPSILON)
	{
		createOverviewColorPlot(dE, 0.0, E_max, "mean_fr_J", "J_{ee} / fC", J_start, J_end, J_step,
									  "'%.1f'", "⟨{/Symbol n}⟩_{exc} / Hz", 5, 4, 6);
		createOverviewColorPlot(dE, dE, E_max, "gauss_sigma_J", "J_{ee} / fC", J_start, J_end, J_step,
									  "'%.1f'", "{/Symbol s}_{FR}", 5, 4, 7);
		createOverviewColorPlot(dE, 0.0, E_max, "mean_fr_inh_J", "J_{ee} / fC", J_start, J_end, J_step,
									  "'%.1f'", "⟨{/Symbol n}⟩_{inh} / Hz", 5, 4, 8);
		createOverviewColorPlot(dE, dE, E_max, "gauss_sigma_inh_J", "J_{ee} / fC", J_start, J_end, J_step,
									  "'%.1f'", "{/Symbol s}_{FR}^{inh}", 5, 4, 9);
	}

	// Seeking I_const
#ifdef SEEK_I_CONST//----------------------------------------------------------------------
	ofstream seekICData(dateStr("_I_const_nu.txt"));
	ofstream seekICResult(dateStr("_I_const.txt"));
	double seekICVar; // firing rate value determined in current simulation
	double seekICVar_old = -1.0; // firing rate value determined in previous simulation, initially has to be set to -1.0

	double I_const_start = I_const; // starting with the adjusted I_const
	double I_const_step = 0.001; // current increment step size
	frequency_start = frequency_end; // in case of "seeking I_const", there is no sense in variation of frequency (because no light is being considered)
	E_max = 0.0; // irradiance has to be switched off (otherwise, simulation takes unnecessary long and seekICVar won't contain the right value!)
	vector<double>().swap(stimplot); // empty vector for plot irradiances
#endif//-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

	// Create NetworkSimulation object and set computational parameters
	NetworkSimulation sn = NetworkSimulation(Nl, Nl_inh, dt, dE, raster_bin, t_offset, t_max, E_max, stimplot, rasterplot, &colorPlotData);
#ifdef SEEK_I_CONST//----------------------------------------------------------------------
	sn.setSeekICVar(&seekICVar); // pass address of seekICVar to NetworkSimulation object
#endif//-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

	// Loop over parameters and run simulation
	frequency = frequency_start;
	pc = pc_start;
	tau_syn = tau_syn_start;
	J = J_start;

	while(lrint(frequency*1000.0)
		<= lrint(frequency_end*1000.0)) // "lrint" condition because otherwise, loop sometimes exits before last round
	{
		while(lrint(pc*1000.0)
			<= lrint(pc_end*1000.0))
		{
			while(lrint(tau_syn*1000.0)
				<= lrint(tau_syn_end*1000.0))
			{
				while(lrint(J*1000.0)
					<= lrint(J_end*1000.0))
				{
					// if not fixed, draw I_const value from table of calibrated values
					if (use_cal_table)
					{
						int ctts;
						if (abs(tau_syn - 1.0) < EPSILON)
							ctts = 0;
						else if (abs(tau_syn - 5.0) < EPSILON)
							ctts = 1;
						else if (abs(tau_syn - 10.0) < EPSILON)
							ctts = 2;
						else
							ctts = -1;
						if (Nl_inh > 0)
							I_const = cal_table [ctts] [int(round(pc*1000.0))-1] [int(round(J/J_step))-1]; // values for mixed network
						else
							I_const = cal_table_exc [ctts] [int(round(pc*1000.0))-1] [int(round(J/J_step))-1]; // values for excitatory network
					}

#ifdef SEEK_I_CONST//----------------------------------------------------------------------
					while(true) // is cancelled by "break" once 5 Hz has been found
					{
#endif//-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
					sn.setParams(t_pulse, frequency, I_const, tau_syn, J, J_ee, J_ei, J_ie, J_ii, pc, sigma_light);
					sn.simulate(path,
#ifdef SEEK_I_CONST
						(frequency == frequency_start && pc == pc_start && tau_syn == tau_syn_start && J == J_start && I_const == I_const_start) ?
#else
						(frequency == frequency_start && pc == pc_start && tau_syn == tau_syn_start && J == J_start) ?
#endif
						true : false); // run the actual simulation

#ifdef SEEK_I_CONST//----------------------------------------------------------------------
						// after the simulation, seekICVar contains the mean firing rate
						seekICData << fixed << I_const << "\t\t\t" << pc << "\t\t\t" << tau_syn << "\t\t\t" << J << "\t\t\t";
						if (seekICVar > 0) // zeros should not be used
							seekICData << seekICVar << endl;
						else
							seekICData << "nan" << endl;

						// if 5 Hz has been crossed, compute I_const using linear interpolation
						if ( (seekICVar_old > 5.0 && seekICVar <= 5.0) ||  									// crossing from above 5.0 Hz
							  (seekICVar_old < 5.0 && seekICVar >= 5.0 && seekICVar_old >= 0.0) )		// crossing from below 5.0 Hz
						{
							double m = (seekICVar_old-seekICVar)/(I_const+I_const_step-I_const);
							double b = seekICVar_old - m*(I_const+I_const_step);
							seekICResult << fixed << pc << "\t\t\t" << tau_syn << "\t\t\t" << J << "\t\t\t" << (5.0 - b) / m << endl;
							break;
						}

						// initial value of I_const caused a firing rate lower than 5 Hz -> increase I_const from now on
						if (seekICVar < 5.0 && seekICVar_old < 0.0) // only if firing rate was below 5.0 and if this was the initial step
							I_const_step *= (-1.0);

						seekICVar_old = seekICVar;
						I_const -= I_const_step;
					}
					seekICVar_old = -1.0; // reset seekICVar_old
					I_const_step = (I_const_step < 0 ? I_const_step * (-1.0) : I_const_step); // reset I_const_step
					//I_const = I_const_start; // reset I_const // do not reset - for next round should start with previous value
					seekICData << endl;
#endif//-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

					J += J_step;
				}
				J = J_start; //reset J
				tau_syn += tau_syn_step;
			}
			tau_syn = tau_syn_start; // reset tau_syn
			pc += pc_step;
		}
		pc = pc_start; // reset pc
		frequency += frequency_step;
	}

	// Finish color plots (only create the necessary ones)
	colorPlotData.close(); // close color plot data file
	if (chdir(path.c_str()) == -1) { // try to re-change directory
		showChDirErrMessage();
		return -1;
	}

	if (frequency_start != frequency_end)
	{
		system("gnuplot mean_fr_freq_cplot.gpl");
		system("gnuplot gauss_sigma_freq_cplot.gpl");
	}
	if (pc_start != pc_end)
	{
		system("gnuplot mean_fr_pc_cplot.gpl");
		system("gnuplot gauss_sigma_pc_cplot.gpl");
	}
	if (tau_syn_start != tau_syn_end)
	{
		system("gnuplot mean_fr_tau_syn_cplot.gpl");
		system("gnuplot gauss_sigma_tau_syn_cplot.gpl");
	}
	if (J_start != J_end)
	{
		system("gnuplot mean_fr_J_cplot.gpl");
		system("gnuplot gauss_sigma_J_cplot.gpl");
	}

	// Finish seeking I_const
	#ifdef SEEK_I_CONST
	seekICData.close();
	seekICResult.close();
	#endif

	return 0;
}
