/**************************************************************************
 * Calling code for the simulation of a single neuron with
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
// Mandatory linker options: -lgsl -lgslcblas
// Compiler option to define Linux platform (default): -D LINUX
// Compiler option to define Windows platform: -D WINDOWS

#include "NeuronSimulation.cpp"

const vector<double> stimplot {5.}; // mW/mm^2 - all light intensities for which plots shall be created
const double plot_start = 1000.0; //ms
const double plot_end = 2000.0; //ms

const double dt = 0.01; // ms - be aware that constants given in time steps cause different behavior if dt is changed!
double dE = 0.5; // mW/mm^2
const int dE_average = 1; // number of bins
double t_max = 20000.0; // ms
double t_offset = 200.0; // ms
double E_max = 10.; // mW/mm^2
double t_pulse = 4.0; // ms

double frequency_step = 5.0; // Hz
double frequency_start = 5.0; // Hz
double frequency_end = 60.0; // Hz

double I_const = 0.914576; // nA - corresponding to sigma_WN = 0.010 nA

/*** main ***/
/* - return: -1 if an error occurred */
int main()
{
	string path = getCurrentDir() + string("/SingleNeuron_") + dateStr("", true); // path to working directory
	string cplotdf = dateStr("_cplotdata.txt"); // color plot data file (with first time stamp in its name)
	double frequency = frequency_start; // current frequency value

#ifdef WINDOWS
	system(""); // Windows needs this to start using colors in the console...
#endif
	cout << getSeparator('=') << endl
	     << "\x1b[36mSimulation Tool for a Single Neuron with Light-Sensitive Channels\x1b[0m" << endl
	     << getSeparator('=') << endl
	     << "Code and documentation: https://jlubo.net/nn-lightchannels-sim" << endl
	     << "Contact: Jannik Luboeinski <jannik.lubo[at]gmx.de>" << endl
	     << getSeparator('-') << endl;

	// demand arguments from user
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

	cout << "Mean input current (in nA, default is " << I_const << " nA): ";
	readInp(I_const, std::numeric_limits<double>::lowest());

	cout << "Simulation offset time (in ms, default is " << t_offset << " ms): ";
	readInp(t_offset);

	cout << "Total duration of one simulation (in ms, default is " << t_max << " ms): ";
	readInp(t_max);

	mkDir(path); // create working directory
	if (chdir(path.c_str()) == -1) { // Try to change directory
		showChDirErrMessage();
		return -1;
	}

	createOverviewColorPlot(dE, dE, E_max, "fr_steady-height", "f / Hz", frequency_start, frequency_end, frequency_step, "'%.0f'",
									"{/Symbol n}_{s,max} / Hz", 2, 1, 3);
	createOverviewColorPlot(dE, dE, E_max, "fr_FWHM", "f / Hz", frequency_start, frequency_end, frequency_step, "'%.0f'",
									"t^{FR}_{FWHM} / ms", 2, 1, 4);
	createOverviewColorPlot(dE, dE, E_max, "fr_peak", "f / Hz", frequency_start, frequency_end, frequency_step, "'%.0f'",
									"{/Symbol n}_{peak} / Hz", 2, 1, 5);
	createOverviewColorPlot(dE, dE, E_max, "fr_steady-base", "f / Hz", frequency_start, frequency_end, frequency_step, "'%.0f'",
									"{/Symbol n}_{s,base} / Hz", 2, 1, 6);
	createOverviewColorPlot(dE, dE, E_max, "fr_tau_adapt", "f / Hz", frequency_start, frequency_end, frequency_step, "'%.0f'",
									"{/Symbol t}^{FR}_{adapt} / ms", 2, 1, 7);

	createOverviewColorPlot(dE, dE, E_max, "O_steady-height", "f / Hz", frequency_start, frequency_end, frequency_step, "'%.0f'",
									"O_{s,max}", 2, 1, 8);
	createOverviewColorPlot(dE, dE, E_max, "O_FWHM", "f / Hz", frequency_start, frequency_end, frequency_step, "'%.0f'",
									"t^{O}_{FWHM} / ms", 2, 1, 9);
	createOverviewColorPlot(dE, dE, E_max, "O_peak", "f / Hz", frequency_start, frequency_end, frequency_step, "'%.0f'",
									"O_{peak}", 2, 1, 10);
	createOverviewColorPlot(dE, dE, E_max, "O_steady-base", "f / Hz", frequency_start, frequency_end, frequency_step, "'%.0f'",
									"O_{s,base}", 2, 1, 11);
	createOverviewColorPlot(dE, dE, E_max, "O_tau_adapt", "f / Hz", frequency_start, frequency_end, frequency_step, "'%.0f'",
									"{/Symbol t}^{O}_{adapt} / ms", 2, 1, 12);

	writePalViridis(); // create Viridis color palette file
	//writePalRainbow(); // create Rainbow color palette file

	NeuronSimulation sn = NeuronSimulation(dt, dE, dE_average, t_offset, t_max, E_max, stimplot, cplotdf, plot_start, plot_end);	// create object and set 																																												// computational parameters

	while(lrint(frequency*1000.0)
		<= lrint(frequency_end*1000.0)) // "lrint" condition because otherwise, loop sometimes exits before last round
	{
		sn.setParams(t_max, E_max, t_pulse, frequency, I_const);
		sn.simulate(path, (frequency == frequency_start) ? true : false);
		frequency += frequency_step;
	}

	if (chdir(path.c_str()) == -1) { // Try to re-change directory
		showChDirErrMessage();
		return -1;
	}

	system("gnuplot fr_steady-height_cplot.gpl");
	system("gnuplot fr_FWHM_cplot.gpl");
	system("gnuplot fr_peak_cplot.gpl");
	system("gnuplot fr_steady-base_cplot.gpl");
	system("gnuplot fr_tau_adapt_cplot.gpl");

	system("gnuplot O_steady-height_cplot.gpl");
	system("gnuplot O_FWHM_cplot.gpl");
	system("gnuplot O_peak_cplot.gpl");
	system("gnuplot O_steady-base_cplot.gpl");
	system("gnuplot O_tau_adapt_cplot.gpl");

	return 0;
}
