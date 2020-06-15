/**************************************************************************
 * Simulation of a neuronal network with Channelrhodopsin-2 channels
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

#include <iostream>
#include <fstream>
#include <unistd.h> // for chdir()

using namespace std;

//#define STIMULUS_COMPARISON // performs a comparison between pulsed light stimulus and pulsed current stimulus // INCOMPLETE AND NOT TESTED
//#define CONNECTION_PLOT // creates a gnuplot file for a plot showing the detailed interneuronal connections within the excitatory network
//#define CONNECTION_PLOT_CREATION // actually creates the plot of interneuronal connections within the excitatory network
#define SPIKES_PLOT // creates a plot of the total number of spikes per raster_bin and a spike raster plot (only for minimum stimulus amplitude!)

#include "Tools.hpp"
#include "Plots.hpp"
#include "Network.hpp"

/*** NetworkSimulation class ***
 * simulates a network of neurons, has instance of Network class */
class NetworkSimulation {

private:

/*** Simulation parameters ***/
const double dt; // ms, one time step for numerical integration
const double dE; // mW/mm^2, irradiance accuracy (in case of STIMULUS_COMPARISON: also dimension nA)
const int Nl; // number of neurons in one line (row or column) of the excitatory population
const int Nl_inh; // number of neurons in one line (row or column) of the inhibitory population
const double center; // the center position of the excitatory neuron distribution
const double center_inh; // the center position of the inhibitory neuron distribution
const int offset; // timesteps of the temporal offset up to the stimulus onset - is not included in t_max; during this time, no synaptic transmission takes place
double pc; // connection probability for unidirectional neuron connections
double tau_syn; // ms, synaptic time constant
double t_max;  // ms, total duration of simulation
double E_max; // mW/mm^2, irradiance to go up to (in case of STIMULUS_COMPARISON: also dimension nA)
double t_pulse; // ms, duration of one stimulus pulse
double frequency; // Hz, frequency for stimulus pulses
double J; // synaptic coupling strength (scale factor for all four different coupling strengths)
Network net;  // the network
#ifdef STIMULUS_COMPARISON
Network net2; // instance of another network
#endif

/*** Output parameters ***/
const vector<double> stimplot; // stimulus amplitudes for which data and plots will be created (only multiples of dE with precision up to .2 are processed)
const vector<int> rasterplot; // neurons for which spikes will be shown in the spike raster plot (if SPIKES_PLOT is defined)
#ifdef SPIKES_PLOT
const double raster_bin; // ms, binning of neuron activities for spike number plot and spike raster plot
#endif
ofstream *contour; // pointer to a data file for creating contour plot (includes freq., intensity, tau_syn, pc, J, mean firing rate and spatial Gaussian std. dev.)
#ifdef SEEK_I_CONST
double *seekic; // pointer to a variable to communicate with NetworkBatch class while seeking I_const
#endif
string purpose; // a small text describing the purpose of this simulation

/*** considerAmplitude ***
 * Checks if given stimulus amplitude must be considered (based on array stimplot) *
 * - E: discretized irradiance (in units of dE)
 * - return: true, if given irradiance must be considered */
bool considerAmplitude(int E) const
{
	for (int i=0; i<stimplot.size(); i++)
	{
		if (int(round(stimplot[i]/dE)) == E)
			return true;
	}
	return false;
}

/*** saveParams ***
 * Saves the crucial parameters in a file *
 * - str: string containing additional information, like the elapsed time */
void saveParams(const char* str)
{
	ofstream f (dateStr("_PARAMS.txt"));
	f << dateStr("") << endl; // time stamp
	f << endl;

	// Parameters
	f << "Simulation parameters:" << endl;
	f << "dt = " << dt << " ms" << endl;
	f << "dE = " << dE << " mW/mm^2" << endl;
	f << "t_max = " << t_max << " ms (" << int(ceil(t_max / dt)) << " steps)" << endl;
	f << "t_offset = " << offset*dt << " ms (" << offset << " steps)" << endl;
	f << "E_max = " << E_max << " mW/mm^2 (" << int(ceil(E_max / dE)) + 1 << " steps)" << endl;
	f << "t_pulse = " << t_pulse << " ms" << endl;
	f << "frequency = " << frequency << " Hz" << endl;
	f << "J = " << J << endl;
	net.saveNetworkParams(&f);
	f << endl;

	// Additional information
	f << "Constant stimulation beyond " << 1000.0/t_pulse << " Hz" << endl;
	f << "Purpose: " << purpose << endl;
	f << str << endl;
	f.close();
}

public:

/*** simulate ***
 * Runs the network simulation *
 * - working_dir: working directory *
 * - first_sim: tells if this is the first simulation run by the batch code */
int simulate(string working_dir, bool first_sim)
{
// ==============================================================================================================================
	// Initialize

	// Start time measurement
	timeMeasure(true);

	// Demand the purpose
	//if (first_sim)
	//{
	//	cout << "Purpose: ";
	//	getline(cin, purpose);
	//}

	// Constants
	const int m = int(ceil(E_max / dE)) + 1; // number of irradiance steps (including zero)
	const int n = int(ceil(t_max / dt)); // number of time steps
	const int pulse_length = int(ceil(t_pulse / dt)); // number of time steps of one pulse
	const int period_steps = int(round((1000/frequency)/dt)); // number of time steps of one period
	const string path = working_dir + "/f=" + dtos(frequency, 1) + ",pc=" + dtos(pc, 3) + ",tau_syn=" + dtos(tau_syn, 0)
									+ ",J=" + dtos(J, 1); // path to working directory
	const string separator = getSeparator(); // string of characters for a separator in command line

	// Try to change directory
	mkDir(path);
	if (chdir(path.c_str()) == -1) {
		showChDirErrMessage();
		cout << separator << endl;
		return -1;
	}

	// Output with general information
	// if this is the first simulation of this batch (first_sim = true), use time stamp from NetworkBatch.cpp, else, set a new time stamp
	cout << "\x1b[33mSimulating network with Ne = " << Nl*Nl << ", Ni = " << Nl_inh*Nl_inh
		  << " for " << t_max << " ms (" << dateStr("", !first_sim) << ")\x1b[0m" << endl;
	cout << "Stimulus: rectangle, t_pulse = " << t_pulse << " ms, \x1b[34mfrequency = " << frequency << " Hz\x1b[0m" << endl;
	cout << "Other parameters: \x1b[32mpc = " << pc << "\x1b[0m, tau_syn = "
#ifdef DELTA_SYNAPSES
		  << 0
#else
		  << tau_syn
#endif
		  << " ms, J = " << J << ", \x1b[35mI_const = " << net.getConstCurrent(1,1) << " nA\x1b[0m" << endl;

	// Declarations and initializations
	double max_firing_rate = 0.0; // contains the maximum firing rate over all exc. neurons and stimulus amplitudes
	double min_firing_rate = numeric_limits<double>::max(); // contains the minimum firing rate over all exc. neurons and stimulus amplitudes
	double max_firing_rate_inh = 0.0; // contains the maximum firing rate over all inh. neurons and stimulus amplitudes
	double min_firing_rate_inh = numeric_limits<double>::max(); // contains the minimum firing rate over all inh. neurons and stimulus amplitudes
	ofstream gpl_mf("mean_fr.gpl"); // gnuplot file for creating a plot of mean firing rate over irradiance/current amplitude in a PDF file
	ofstream gpl_sf("gauss_sigma_fr.gpl"); // gnuplot file for creating a plot of the firing rate std. dev. over irradiance/current amplitude in a PDF file
	ofstream txt_msf(dateStr("_fr.txt")); // data file for creating plots over stimulus amplitude
	ofstream txt_stim; // data for creating map plot of the stimulus for specific stimulus amplitude
	ofstream* txt_fr = new ofstream[stimplot.size()]; // data for creating exc. plots for specific stimulus amplitudes over time
	ofstream* txt_fr_inh = new ofstream[stimplot.size()]; // data for creating inh. plots for specific stimulus amplitudes over time
	ofstream gpl_light_stim; // gnuplot file for creating irradiance map plot for last specific stimulus amplitude (only one is necessary, for all the other
									 // amplitudes the only difference is the color code range, which always goes up to the particular amplitude)
#ifdef STIMULUS_COMPARISON
	ofstream gpl_curr_stim; // gnuplot file for creating current stimulus map plot for last specific stimulus amplitude (only one is necessary, for all the other
									// amplitudes the only difference is the color code range, which always goes up to the particular amplitude)
#endif
	ofstream* gpl_lfr = new ofstream[stimplot.size()]; // gnuplot files for creating exc. firing rate map plots for specific light stimulus amplitudes
	ofstream* gpl_lfr_inh = new ofstream[stimplot.size()]; // gnuplot files for creating inh. firing rate map plots for specific light stimulus amplitudes
	ofstream* gpl_fit_light = new ofstream[stimplot.size()]; // gnuplot files for Gauss sigma fit with light stimulus for specific light stimulus amplitudes
#ifdef STIMULUS_COMPARISON
	ofstream* gpl_cfr = new ofstream[stimplot.size()]; // gnuplot files for creating exc. firing rate map plots for specific current stimulus amplitudes
	ofstream* gpl_cfr_inh = new ofstream[stimplot.size()]; // gnuplot files for creating inh. firing rate map plots for specific current stimulus amplitudes
	ofstream* gpl_fit_current = new ofstream[stimplot.size()]; // gnuplot files for Gauss sigma fit with current stimulus for specific current stimulus amplitudes
#endif
#ifdef CONNECTION_PLOT
	ofstream gpl_cp("connections_arrows.gpl"); // gnuplot file for creating a plot of the interneuronal connections
#endif
	ofstream gpl_cn("inc_connection_map.gpl"); // gnuplot file for creating a color plot of the numbers of incoming interneuronal connections
	ofstream txt_cn(dateStr("_inc_connection.txt")); // data file for creating a color plot of the numbers of incoming interneuronal connections
	ofstream gplscript("gpl"); // shell script for calling all the gnuplot scripts (to re-generate the plots, e.g., for another gnuplot version)
#ifdef SPIKES_PLOT
	ofstream gpl_spikes("spike_num.gpl"); // gnuplot file for creating a plot of the number of spikes over time
	ofstream gpl_srp("spike_raster_plot.gpl"); // gnuplot file for creating a spike raster plot in a PDF file
	ofstream txt_spikes(dateStr("_spikes.txt")); // data file containing total number of spikes and activities of specific neurons over time
#endif
	ofstream logf("gslfit.log"); // log file containing information about data processing

	writePalViridis(); // create Viridis color palette file

	Stimulus lst = Stimulus(period_steps); // light stimulus
#ifdef STIMULUS_COMPARISON
	Stimulus cst = Stimulus(period_steps); // current stimulus
#endif

	double p = 0.0; // percentage of process completeness
	int current_E = 0; // specifies current position (index) in stimplot[] array

// ==============================================================================================================================
	// Prepare files for writing

	if (stimplot.size() > 0)
	{
		// Prepare file for irradiance map plot for specific stimulus amplitude
		string filename;
		filename = "_irradiance_" + dtos(stimplot.back(), 2) + ".txt"; // choose last specific stimulus amplitude (stimplot.back()) for this plot
		txt_stim.open(dateStr(filename));
		filename = "irradiance_" + dtos(stimplot.back(), 2) + "_map.gpl";
		gpl_light_stim.open(filename);
#ifdef STIMULUS_COMPARISON
		filename = "currentstim_" + dtos(stimplot.back(), 2) + "_map.gpl";
		gpl_curr_stim.open(filename);
#endif
		if (!txt_stim.is_open() || !gpl_light_stim.is_open()
#ifdef STIMULUS_COMPARISON
			 || !gpl_curr_stim.is_open()
#endif
												 )
		{
			cout << "Unable to open file!" << endl << separator << endl;
			return -1;
		}

		// Prepare files for specific stimulus plots (specified in "stimplot" vector)
		for (int i=0; i<stimplot.size(); i++)
		{
			filename = "_fr_exc_" + dtos(stimplot[i], 2) + ".txt";
			txt_fr[i].open(dateStr(filename)); // open data file
			filename = "fr_exc_" + dtos(stimplot[i], 2) + "_map_light.gpl";
			gpl_lfr[i].open(filename); // open gnuplot file for firing rate map (exc. population) with light stimulus

			filename = "_fr_inh_" + dtos(stimplot[i], 2) + ".txt";
			txt_fr_inh[i].open(dateStr(filename)); // open data file
			filename = "fr_inh_" + dtos(stimplot[i], 2) + "_map_light.gpl";
			gpl_lfr_inh[i].open(filename); // open gnuplot file for firing rate map (inh. population) with light stimulus

			filename = "gauss_sigma_fit_" + dtos(stimplot[i], 2) + "_light.gpl";
			gpl_fit_light[i].open(filename); // open gnuplot file for Gauss sigma fit with light stimulus
#ifdef STIMULUS_COMPARISON
			filename = "fr_exc_" + dtos(stimplot[i], 2) + "_map_current.gpl";
			gpl_cfr[i].open(filename); // open gnuplot file for firing rate map (exc. population) with current stimulus

			filename = "fr_inh_" + dtos(stimplot[i], 2) + "_map_current.gpl";
			gpl_cfr_inh[i].open(filename); // open gnuplot file for firing rate map (inh. population) with current stimulus

			filename = "gauss_sigma_fit_" + dtos(stimplot[i], 2) + "_current.gpl";
			gpl_fit_current[i].open(filename); // open gnuplot file for Gauss sigma fit with current stimulus
#endif
			if ( !txt_fr[i].is_open() || !gpl_lfr[i].is_open() || !txt_fr_inh[i].is_open() || !gpl_lfr_inh[i].is_open() || !gpl_fit_light[i].is_open()
#ifdef STIMULUS_COMPARISON
					|| !gpl_cfr[i].is_open() || !gpl_cfr_inh[i].is_open() || !gpl_fit_current[i].is_open()
#endif
				)
			{
				cout << "Unable to open file!" << endl << separator << endl;
				return -1;
			}
		}
	}

	// Check if files have been opened properly
  	if ( !gplscript.is_open() || !gpl_mf.is_open() || !txt_msf.is_open() || ((contour != NULL) && (!contour->is_open()))	)
	{
		cout << "Unable to open file!" << endl << separator << endl;
		return -1;
	}

	if (
#ifdef CONNECTION_PLOT
		!gpl_cp.is_open() ||
#endif
#ifdef SPIKES_PLOT
		!gpl_spikes.is_open() ||
		!gpl_srp.is_open() ||
		!txt_spikes.is_open() ||
#endif
		!gpl_cn.is_open() || !txt_cn.is_open())
	{
		cout << "Unable to open file!" << endl << separator << endl;
		return -1;
	}


// ==============================================================================================================================
	// Output script for connection plot

	int total_c_count_exc = 0; // total number of exc. connections
	int total_c_count_inh = 0; // total number of inh. connections

#ifdef CONNECTION_PLOT
	int color_n = 0;
	char color[10];
	gpl_cp << "set term pdf enhanced font \"Sans, 12\" color lw 1" << endl;
	gpl_cp << "set output '" << dateStr("_connections.pdf'") << endl << endl;
	gpl_cp << "set xrange [1" << ":" << Nl << "]" << endl;
	gpl_cp << "set yrange [1" << ":" << Nl << "]" << endl;
	gpl_cp << "set xlabel \"Neuron\"" << endl;
	gpl_cp << "set ylabel \"Neuron\"" << endl;
#endif

	// loops over k and l replace one loop over m in consecutive numbering
	for (int k=1; k<=Nl; k++) // row/column neuron numbering
	{
		for (int l=1; l<=Nl; l++) // row/column neuron numbering
		{
			int m = cNN(k,l); // achieve consecutive number
			const int c_count_exc = net.getNumberIncoming(TYPE_EXC,k,l); // number of incoming excitatory connections to neuron m
			const int c_count_inh = net.getNumberIncoming(TYPE_INH,k,l); // number of incoming inhibitory connections to neuron m

#ifdef CONNECTION_PLOT
			// draw connections between excitatory neurons
			for (int n=0; n<Nl*Nl; n++) // consecutive neuron numbering
			{
				switch (color_n)
				{
					case 0:
						strcpy(color, "red");
						break;
					case 1:
						strcpy(color, "blue");
						break;
					case 2:
						strcpy(color, "green");
						break;
					case 3:
						strcpy(color, "purple");
						break;
					default:
						strcpy(color, "orange");
						color_n = -1;
						break;
				}

				if (net.areConnected(n, m)) // check if there is connection n->m
				{
					gpl_cp << "set arrow from " << net.rowNN(n) << "," << net.colNN(n) << " to " << k << "," << l
							 << " head size 0.2,20 filled back lc rgb '" << color << "'" << endl;

					if (net.areConnected(m, n) && (m > n)) // check if there is also a connection m->n (if yes, draw a yellow dashed second arrow),
																		// but not to do this twice, only consider lower triangular matrix of connection matrix (m > n condition)
					{
						gpl_cp << "set arrow from " << k << "," << l << " to " << net.rowNN(n) << "," << net.colNN(n)
								 << " head size 0.2,20 filled back lc rgb 'yellow' dt 3" << endl;
					}
					color_n++;
				}
			}

			// draw neurons
			if (c_count_exc > 0) // if neurons possesses incoming connections, draw filled circle
				gpl_cp << "set object circle at first " << k << "," << l
					 	 << " radius char 0.15 front fillcolor rgb 'black' fillstyle solid noborder" << endl;
			else  // if not, draw empty circle
				gpl_cp << "set object circle at first " << k << "," << l
						  << " radius char 0.15 front fillstyle empty border lc rgb 'black' lw 1" << endl;
#endif //CONNECTION_PLOT

			total_c_count_exc += c_count_exc;
			total_c_count_inh += c_count_inh;

			txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << c_count_exc+c_count_inh << endl;
		} // end of for(l)
		txt_cn << endl; // empty line for color plot
	}

#ifdef CONNECTION_PLOT
	gpl_cp << "plot 0 notitle lc rgb 'white'" << endl;
	gpl_cp.close();
#ifdef CONNECTION_PLOT_CREATION
	system("gnuplot connections_arrows.gpl"); // resulting PDF file needs several MB of disk space
#endif
#endif //CONNECTION_PLOT

	txt_cn.close();
	createNetworkColorPlot(gpl_cn, Nl, -1.0, 3, "inc_connection", "", true, "# of incoming connections");
	gplscript << "gnuplot inc_connection_map.gpl" << endl;

	cout << "Number of E->E connections: " << total_c_count_exc << " (expected: " << pc * (pow(Nl*Nl, 2)-(Nl*Nl)) << ")" << endl;
	cout << "Number of I->E connections: " << total_c_count_inh << " (expected: " << pc * pow(Nl, 2) * pow(Nl_inh, 2) << ")" << endl;

// ==============================================================================================================================
	// (Light/current) stimulation loop

	for(int i = 0; i < m; i++)
	{

		// Reset network
		net.reset();
#ifdef STIMULUS_COMPARISON
		net2.reset();
#endif

		// Clear stimuli
		lst.clear();
#ifdef STIMULUS_COMPARISON
		cst.clear();
#endif
#ifdef SPIKES_PLOT
		int total_spike_num = 0; // total number of spikes in currently considered interval
		bool* spike = new bool[rasterplot.size()]; // spiking activity of specific neurons in currently considered interval
		for (int srn = 0; srn < rasterplot.size(); srn++)
			spike[srn] = false;
#endif

		// Time loop
		for (int j = 0; j <= n+offset; j++)
		{
			// Update percentage
			double p_new = double(round(double(i*(n+offset) + j) / double(m*(n+offset)) * 1000.0)) / 10.0; // round to first decimal place
			if (p_new > p)
			{
				p = p_new;
				printf("\rProgress: %.1f %% completed.", p);
				fflush(stdout);
			}

			// Stimuli start at offset time
			if (j == offset)
			{
				// Set light stimulus
				lst.addRectPulse(i*dE, 0, pulse_length);
				net.setGaussianLightStimulus(lst, center, center); // set amplitude of light stimulus in the center of the network

				// Set current stimulus
#ifdef STIMULUS_COMPARISON
				cst.addRectPulse(net2.getSteadyCurrent(i*dE, 0), 0, pulse_length); // use steady current (same stimulus as for light is not possible!)
				net2.setGaussianCurrentStimulus(cst, center, center); // set amplitude of light stimulus in the center of the network
#endif
			}

			// Calculate next step for Network (and implicitly for Neurons and Channelrhodopsin)
#ifdef SPIKES_PLOT
			total_spike_num +=
#endif
				net.processTimeStep(j); // necessarily negative time steps during offset - those spikes will not be counted
#ifdef SPIKES_PLOT
			// Detect spikes in this interval for raster plot
			for (int srn = 0; srn < rasterplot.size(); srn++)
			{
				if (net.getActivity(rasterplot[srn])) // if once in the interval a spike occurs
					spike[srn] = true;
			}

			// Write spike data to file
			if ((i == 0) && (j-offset) % int(raster_bin/dt) == 0)
			{
				txt_spikes << (j-offset)*dt << "\t\t\t" << total_spike_num; // time and number of spikes within interval raster_bin/dt
				total_spike_num = 0; // reset total spike count

				// Write raster plot data
				for (int srn = 0; srn < rasterplot.size(); srn++)
				{
					txt_spikes << "\t\t\t" << spike[srn]; // activity of specific neurons within interval ("has been active or not")
					spike[srn] = false;
				}
				txt_spikes << endl;
			}
#endif


#ifdef STIMULUS_COMPARISON
			net2.processTimeStep(j);
#endif

		} // end of for(j)
#ifdef SPIKES_PLOT
		delete[] spike;
#endif
// ==============================================================================================================================
		// Compute mean firing rate, standard deviation for the firing rates and fit the Gaussian distribution
			//  FITTING DOES NOT YET WORK PROPERLY FOR INHIB. POPULATION

		double mfr_light = 0.; // mean firing rate of the whole excitatory population in the light-stimulated network
		double mfr_light_inh = 0.; // mean firing rate of the whole inhibitory population in the light-stimulated network
		double sdfr_light = 0.; // standard deviation of the firing rate of the whole excitatory population in the light-stimulated network
		double sdfr_light_inh = 0.; // standard deviation of the firing rate of the whole inhibitory population in the light-stimulated network
		int sdfr_light_part1 = 0; // part 1 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_light_part2 = 0; // part 2 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_light_part1_inh = 0; // part 1 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_light_part2_inh = 0; // part 2 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int denom_light = 4; // denominator for sum of slice results

		double sigma_gauss_light = 0.; // standard deviation of Gaussian firing rate distribution evoked by light stimulus (not to be confused with standard deviation 												 // of firing rates!)
		double sigma_gauss_err_light = 0.; // error of the standard deviation of Gaussian firing rate distribution evoked by light stimulus
		double A_gauss_light = 0.; // amplitude of Gaussian firing rate distribution evoked by light stimulus									 // of firing rates!)
		double A_gauss_err_light = 0.; // error of the amplitude of Gaussian firing rate distribution evoked by light stimulus
		double b_gauss_light = 0.; // baseline of Gaussian firing rate distribution evoked by light stimulus										 // of firing rates!)
		double b_gauss_err_light = 0.; // error of the baseline of Gaussian firing rate distribution evoked by light stimulus
		double sigma_gauss_light_inh = 0.; // standard deviation of Gaussian firing rate distribution in inhibitory population evoked by light stimulus (not to be
													  // confused with standard deviation of firing rates!)
		double sigma_gauss_err_light_inh = 0.; // error of the standard deviation of Gaussian firing rate distribution in inhibitory population evoked by light stimulus
#ifdef STIMULUS_COMPARISON
		double mfr_current = 0.; // mean firing rate of the whole excitatory population in the current-stimulated network
		double mfr_current_inh = 0.; // mmean firing rate of the whole inhibitory population in the current-stimulated network
		double sdfr_current = 0.; // standard deviation of the firing rate of the whole excitatory population in the light-stimulated network
		double sdfr_current_inh = 0.; // standard deviation of the firing rate of the whole inhibitory population in the light-stimulated network
		int sdfr_current_part1 = 0; // part 1 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_current_part2 = 0; // part 2 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_current_part1_inh = 0; // part 1 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_current_part2_inh = 0; // part 2 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int denom_current = 4; // denominator for sum of slice results

		double sigma_gauss_current = 0.; // standard deviation of Gaussian firing rate distribution evoked by current stimulus (not to be confused with standard 													// deviation of firing rates!)
		double sigma_gauss_err_current = 0.; // error of the standard deviation of Gaussian firing rate distribution evoked by current stimulus
		double sigma_gauss_current_inh = 0.; // standard deviation of Gaussian firing rate distribution in inhibitory population evoked by current stimulus (not to be
															// confused with standard deviation of firing rates!)
		double sigma_gauss_err_current_inh = 0.; // error of the standard deviation of Gaussian firing rate distribution in inhibitory population evoked by current stimulus

#endif
		double fit_radius[Nl*Nl]; // x-values for fitting (radius)
		double fit_radius_inh[Nl_inh*Nl_inh]; // x-values for fitting (radius) in inh. population
		double fit_fr[Nl*Nl]; // y-values for fitting (firing rate of light-stimulated network)
		double fit_fr_inh[Nl_inh*Nl_inh]; // y-values for fitting (firing rate of light-stimulated network) in inh. population
#ifdef STIMULUS_COMPARISON
		double fit_fr2[Nl*Nl]; // y-values for fitting (firing rate of current-stimulated network)
#endif
		double fit_fr_err[Nl*Nl]; fill_n(fit_fr_err, Nl*Nl, 0.2); // data error estimates for firing rates -- estimate by 0.2
		double fit_fr_err_inh[Nl*Nl]; fill_n(fit_fr_err_inh, Nl_inh*Nl_inh, 0.2); // data error estimates for firing rates -- estimate by 0.2

		double b = 0.; // baseline shift of the Gaussian - should be zero, just as in the Gaussian stimulus
		double var_init = pow(net.getGaussSigma(),2); // use variance of Gaussian stimulus distribution as initial value
		double par_light[3] = {10., var_init, b}; // set initial values for A & var and fixed value for b (for light stimulation)
		double par_light_inh[3] = {0., var_init, b}; // set initial values for A & var and fixed value for b (for light stimulation)
		double err_light[3]; // parameter error array (for light stimulation)
#ifdef STIMULUS_COMPARISON
		double par_curr[3] = {0., var_init, b}; // set initial values for A & var and fixed value for b (for current stimulation)
		double err_curr[3]; // parameter error array (for current stimulation)
#endif
		double chisq = 0.; // chi squared of fit function

		for (int m=0; m<pow(Nl,2)+pow(Nl_inh,2); m++)
		{
			int nu = net.getSpikeCount(m); // get number of spikes of neuron m in light-stimulated network
#ifdef STIMULUS_COMPARISON
			int nu2 = net2.getSpikeCount(m); // get number of spikes of neuron m in current-stimulated network
#endif

			if (m < pow(Nl,2)) // neuron m is in excitatory population
			{
				mfr_light += double(nu);
				sdfr_light_part1 += pow(nu,2);
				sdfr_light_part2 += nu;

				fit_radius[m] = sqrt(pow(rowNN(m)-center,2) + pow(colNN(m)-center,2)); // store radius in fit data array
				fit_fr[m] = nu/(t_max/1000.);	// store firing rate in fit data array

#ifdef STIMULUS_COMPARISON
				mfr_current += double(nu2);
				sdfr_current_part1 += pow(nu2,2);
				sdfr_current_part2 += nu2;

				fit_fr2[m] = nu2/(t_max/1000.);	// store firing rate in fit data array
#endif
			}
			else // neuron m is in inhibitory population
			{
				int mprime = m - pow(Nl,2);

				mfr_light_inh += double(nu);
				sdfr_light_part1_inh += pow(nu,2);
				sdfr_light_part2_inh += nu;

				fit_radius_inh[mprime] = sqrt(pow(rowInh(mprime)-center_inh,2) + pow(colInh(mprime)-center_inh,2)); // store radius in fit data array
				fit_fr_inh[mprime] = nu/(t_max/1000.);	// store firing rate in fit data array

#ifdef STIMULUS_COMPARISON
				mfr_current_inh += double(nu2);
				sdfr_current_part1_inh += pow(nu2,2);
				sdfr_current_part2_inh += nu2;
#endif
			}
		}

		if (i > 0) // computation of Gaussian sigma only makes sense if there is a stimulus
		{
			// findFit for excitatory population
			gsl_multifit_function_fdf func; // fit function plus its Jacobian (fdf)
			fitdata_pairs fd = {(size_t)Nl*Nl, fit_radius, fit_fr, fit_fr_err}; // structure of arrays for the data points (x, y, yerror)

			// set the properties of the fit function
			func.f = &gauss_f; // the fit function
			func.df = &gauss_df_A_var; // the Jacobian of the fit function
			func.fdf = &gauss_fdf_A_var; // function calling _f and _df functions
			func.n = Nl*Nl; // number of data points
			func.p = 3; // number of parameters (both fit parameters and constant parameters)
			func.params = &fd; // fitdata structure (contains arrays for (x, y, yerror))

			if (findFit(func, par_light, err_light, &chisq)) // do the fit procedure
			{
				sigma_gauss_light = sqrt(abs(par_light[1])); // read out sigma_gauss
				sigma_gauss_err_light = sqrt(abs(err_light[1])); // read out the sigma_gauss error
				A_gauss_light = par_light[0]; // read out A_gauss
			 	A_gauss_err_light = err_light[0]; // read out the A_gauss error
				b_gauss_light = par_light[2]; // read out b_gauss
				b_gauss_err_light = err_light[2]; // read out the b_gauss error

				logf << "Irradiance " << double(i)*dE << ": fit succeeded"
				     << ", A = " << A_gauss_light << " +- " << A_gauss_err_light
				     << ", sigma = " << sigma_gauss_light << " +- " << sigma_gauss_err_light
				     << ", b = " << b_gauss_light << " +- " << b_gauss_err_light
				     << ", chi_squared = " << chisq << endl;
			}
			else
			{
				logf << "Irradiance " << double(i)*dE << ": fit error" << endl;
			}

			if (Nl_inh > 0)
			{
				// findFit for inhibitory population (does not yet work properly!)
				fitdata_pairs fd_inh = {(size_t)Nl_inh*Nl_inh, fit_radius_inh, fit_fr_inh, fit_fr_err_inh}; // structure of arrays for the data points (x, y, yerror)

				// set the properties of the fit function
				func.f = &gauss_f; // the fit function
				func.df = &gauss_df_A_var; // the Jacobian of the fit function
				func.fdf = &gauss_fdf_A_var; // function calling _f and _df functions
				func.n = Nl_inh*Nl_inh; // number of data points
				func.p = 3; // number of parameters (both fit parameters and constant parameters)
				func.params = &fd_inh; // fitdata structure (contains arrays for (x, y, yerror))

				if (findFit(func, par_light_inh, err_light, NULL)) // do the fit procedure
				{
					sigma_gauss_light_inh = sqrt(abs(par_light_inh[1])); // read out sigma_gauss
					sigma_gauss_err_light_inh = sqrt(abs(err_light[1])); // read out the sigma_gauss error

					logf << "Irradiance " << double(i)*dE << " (inh. population): fit succeeded"
					     << ", A = " << A_gauss_light << " +- " << A_gauss_err_light
					     << ", sigma = " << sigma_gauss_light << " +- " << sigma_gauss_err_light
			 		     << ", b = " << b_gauss_light << " +- " << b_gauss_err_light
					     << ", chi_squared = " << chisq << endl;
				}
				else
				{
					logf << "Irradiance " << double(i)*dE << ": fit error (inh. population)" << endl;
				}
			}
		}

		mfr_light /= (double(Nl*Nl)*t_max / 1000.0); // compute the mean firing rate from the number of spikes
		if (Nl_inh > 0)
			mfr_light_inh /= (double(Nl_inh*Nl_inh)*t_max / 1000.0); // compute the mean firing rate from the number of spikes
		sdfr_light = sqrt( ( double(sdfr_light_part1) - double(pow(sdfr_light_part2, 2)) / pow(Nl, 2) ) / double(Nl*Nl-1) ) / (t_max / 1000.0); // compute standard 																														// deviation of the firing rates according to Steiner's translation theorem
		if (Nl_inh > 0)
			sdfr_light_inh = sqrt( ( double(sdfr_light_part1_inh) - double(pow(sdfr_light_part2_inh, 2)) / pow(Nl_inh, 2) ) /
										double(Nl_inh*Nl_inh-1) ) / (t_max / 1000.0); // compute standard deviation of the firing rates according to Steiner's translation 																										 // theorem
#ifdef STIMULUS_COMPARISON
		mfr_current /= (double(Nl*Nl)*t_max / 1000.0);
		sdfr_current = sqrt( ( double(sdfr_current_part1) - double(pow(sdfr_current_part2, 2)) / pow(Nl, 2) ) / double(Nl*Nl-1) ) / (t_max / 1000.0);
		if (Nl_inh > 0)
			sdfr_current_inh = sqrt( ( double(sdfr_current_part1_inh) - double(pow(sdfr_current_part2_inh, 2)) / pow(Nl_inh, 2) )
										/ double(Nl_inh*Nl_inh-1) ) / (t_max / 1000.0);
#endif

// ==============================================================================================================================
		// Output in data files

		// Write mean firing rate, standard deviation of the firing rate and Gauss sigma of the firing rate (over irradiance/current) in data file
		txt_msf << fixed << double(i)*dE << "\t\t\t" //1 Light intensity
				    << dtos(A_gauss_light, 6, true) << "\t\t\t" //2 Gaussian amplitude
				    << dtos(A_gauss_err_light, 6, true) << "\t\t\t" //3 Fit error for Gaussian amplitude of the firing rate
				    << dtos(sigma_gauss_light, 6, true) << "\t\t\t" //4 Gaussian sigma of the firing rate
				    << dtos(sigma_gauss_err_light, 6, true) << "\t\t\t" //5 Fit error for Gaussian sigma of the firing rate
				    << dtos(b_gauss_light, 6, true) << "\t\t\t" //6 Gaussian baseline of the firing rate
				    << dtos(b_gauss_err_light, 6, true) << "\t\t\t" //7 Fit error for Gaussian baseline of the firing rate
				    << dtos(mfr_light, 6, true) << "\t\t\t" //8 Mean firing rate
				    << dtos(sdfr_light, 6, true) << "\t\t\t" //9 Standard deviation of the mean firing rate
				    << dtos(mfr_light_inh, 6, true) << "\t\t\t" //10 Mean firing rate (inh)
				    << dtos(sdfr_light_inh, 6, true); //11 Standard deviation of the mean firing rate (inh.)


#ifdef STIMULUS_COMPARISON
		txt_msf << "\t\t\t" << net2.getSteadyCurrent(i*dE, 0) << "\t\t\t" //12 Current amplitude
		        << dtos(mfr_current, 6, true) << "\t\t\t" //13 Mean firing rate
				    << dtos(sdfr_current, 6, true) << "\t\t\t" //14 Standard deviation of the mean firing rate
				    << dtos(sigma_gauss_current, 6, true) << "\t\t\t" //15 Gaussian sigma of the firing rate
				    << dtos(sigma_gauss_err_current, 6, true) << "\t\t\t" //16 Fit error for Gaussian sigma of the firing rate
				    << dtos(mfr_current_inh, 6, true) << "\t\t\t" //17 Mean firing rate (inh. population)
				    << dtos(sdfr_current_inh, 6, true) << "\t\t\t" //18 Standard deviation of the mean firing rate (inh. population)
				    << dtos(sigma_gauss_current_inh, 6, true) << "\t\t\t" //19 Gaussian sigma of the firing rate (inh. population)
				    << dtos(sigma_gauss_err_current_inh, 6, true); //20 Fit error for Gaussian sigma of the firing rate (inh. population)
#endif

		txt_msf << endl;


		// Write contour plot file
		if (contour != NULL)
		{
			*contour << fixed << frequency << "\t\t\t" << 100.*pc << "\t\t\t" << tau_syn << "\t\t\t" << J << "\t\t\t" << double(i)*dE << "\t\t\t"
			         << dtos(A_gauss_light, 6, true) << "\t\t\t"
 			         << dtos(sigma_gauss_light, 6, true) << "\t\t\t"
			         << dtos(b_gauss_light, 6, true) << "\t\t\t"
			         << dtos(mfr_light, 6, true);

#ifdef STIMULUS_COMPARISON
			*contour << "\t\t\t" << dtos(mfr_current, 6, true)
			         << "\t\t\t" << dtos(sigma_gauss_current, 6, true)
			         << "\t\t\t" << dtos(mfr_current_inh, 6, true)
			         << "\t\t\t" << dtos(sigma_gauss_current_inh, 6, true);
#endif
			*contour << endl;
		}

#ifdef SEEK_I_CONST
		*seekic = mfr_light; // shall be done only once because in "SEEK_I_CONST" mode only irradiance E=0 is used
#endif


// ==============================================================================================================================
		// Preparations for output scripts for firing rate map plots, Gauss sigma fit plots and irradiance map plot (for specific stimulus amplitudes)
		if (considerAmplitude(i))
		{
			for (int m=0; m < pow(Nl,2)+pow(Nl_inh,2); m++)
			{
				double fr = 1000. * double(net.getSpikeCount(m)) / t_max;
#ifdef STIMULUS_COMPARISON
				double fr2 = 1000. * double(net2.getSpikeCount(m)) / t_max;
#endif

				if (m < pow(Nl,2))	// neuron m is in excitatory population
				{
					if (fr > max_firing_rate)
						max_firing_rate = fr;
					if (fr < min_firing_rate)
						min_firing_rate = fr;
#ifdef STIMULUS_COMPARISON
					if (fr2 > max_firing_rate)
						max_firing_rate = fr2;
					if (fr2 < min_firing_rate)
						min_firing_rate = fr2;
#endif

					txt_fr[current_E] << fixed << rowNN(m) << "\t\t\t" << colNN(m) << "\t\t\t" << fr
#ifdef STIMULUS_COMPARISON
											<< "\t\t\t" << fr2
#endif
											<< endl;

					if (current_E == stimplot.size()-1) // data output for light stimulus plot (for last specific amplitude)
					{
						txt_stim << fixed << rowNN(m) << "\t\t\t" << colNN(m) << "\t\t\t" << net.getLightStimulusAt(offset, m) << endl;
						if ((m+1) % Nl == 0)
							txt_stim << endl; // another free line necessary for contour plot
					}
					if ((m+1) % Nl == 0)
						txt_fr[current_E] << endl; // another free line necessary for contour plot
				}
				else // neuron m is in inhibitory population
				{
					if (fr > max_firing_rate_inh)
						max_firing_rate_inh = fr;
					if (fr < min_firing_rate_inh)
						min_firing_rate_inh = fr;
#ifdef STIMULUS_COMPARISON
					if (fr2 > max_firing_rate_inh)
						max_firing_rate_inh = fr2;
					if (fr2 < min_firing_rate_inh)
						min_firing_rate_inh = fr2;
#endif
					int m_eff = m - pow(Nl, 2);
					txt_fr_inh[current_E] << fixed << rowInh(m_eff) << "\t\t\t" << colInh(m_eff) << "\t\t\t" << fr
#ifdef STIMULUS_COMPARISON
												 << "\t\t\t" << fr2
#endif
												 << endl;

					if ((m_eff+1) % Nl_inh == 0)
						txt_fr_inh[current_E] << endl; // another free line necessary for contour plot
				}
			}

			// Close firing rate map data files
			txt_fr[current_E].close();
			txt_fr_inh[current_E].close();

			// Create fit plot (Gaussian FR distribution) of excitatory population for light stimulus
			createGaussFitPlot(gpl_fit_light[current_E], A_gauss_light, sigma_gauss_light, b_gauss_light, i*dE, center, "_light");
			gplscript << "gnuplot gauss_sigma_fit_" << dtos(i*dE,2) << "_light.gpl" << endl;

#ifdef STIMULUS_COMPARISON
			// Create fit plot (Gaussian FR distribution) of excitatory population for current stimulus
			createGaussFitPlot(gpl_fit_light[current_E], par_curr[0], sigma_gauss_current, par_curr[2], i*dE, center, "_current");
			gplscript << "gnuplot gauss_sigma_fit_" << dtos(i*dE,2) << "_current.gpl" << endl;
#endif

			// Create irradiance plot for light stimulus
			if (current_E == stimplot.size()-1)
			{
				txt_stim.close(); // close data file
				createNetworkColorPlot(gpl_light_stim, Nl, stimplot.back(), 3, "irradiance", "", true, "E / mW/mm²", 0., stimplot.back());
				gplscript << "gnuplot irradiance_" << dtos(stimplot.back(),2) << "_map.gpl" << endl;
#ifdef STIMULUS_COMPARISON
				createNetworkColorPlot(gpl_curr_stim, Nl, stimplot.back(), 4, "currentstim", "", true, "E / mW/mm²", 0., stimplot.back());
				gplscript << "gnuplot currentstim_" << dtos(stimplot.back(),2) << "_map.gpl" << endl;
#endif
			}
			current_E++;
		}
	} // end of for(i)
	cout << endl;

// ==============================================================================================================================
	// Create output scripts for firing rate map plots (do this here to have the final min_firing_rate and max_firing_rate)
	for (int e=0; e<stimplot.size(); e++)
	{
		// Create firing rate plot of excitatory population for light stimulus
		createNetworkColorPlot(gpl_lfr[e], Nl, stimplot[e], 3, "fr_exc", "_light", true, "{/Symbol n} / Hz", min_firing_rate, max_firing_rate);
		gplscript << "gnuplot fr_exc_" << dtos(stimplot[e],2) << "_map_light.gpl" << endl;

		// Create firing rate plot of inhibitory population for light stimulus
		createNetworkColorPlot(gpl_lfr_inh[e], Nl_inh, stimplot[e], 3, "fr_inh", "_light", true, "{/Symbol n} / Hz", min_firing_rate_inh, max_firing_rate_inh);
		gplscript << "gnuplot fr_inh_" << dtos(stimplot[e],2) << "_map_light.gpl" << endl;

#ifdef STIMULUS_COMPARISON
		// Create firing rate plot of excitatory population for external current stimulus
		createNetworkColorPlot(gpl_cfr[e], Nl, stimplot[e], 4, "fr_exc", "_current", true, "{/Symbol n} / Hz", min_firing_rate, max_firing_rate);
		gplscript << "gnuplot fr_exc_" << dtos(stimplot[e],2) << "_map_current.gpl" << endl;

		// Create firing rate plot of inhibitory population for external current stimulus
		createNetworkColorPlot(gpl_cfr_inh[e], Nl_inh, stimplot[e], 4, "fr_inh", "_current", true, "{/Symbol n} / Hz", min_firing_rate_inh, max_firing_rate_inh);
		gplscript << "gnuplot fr_inh_" << dtos(stimplot[e],2) << "_map_current.gpl" << endl;
#endif

	}

// ==============================================================================================================================
	// Create plot for spikes per millisecond
#ifdef SPIKES_PLOT
	gpl_spikes << "set term pdf enhanced font \"Sans, 20\" color solid lw 2.5" << endl
				  << "set output '" << dateStr("_spike_num.pdf'") << endl
				  << "set xlabel \"time / ms\"" << endl
				  << "set ylabel \"spikes/ms\"" << endl
				  << "set tmargin at screen 0.95" << endl
				  << "set bmargin at screen 0.23" << endl << endl
				  << "plot [x=0:" << t_max << "] '" << dateStr("_spikes.txt") << "' \\" << endl
				  << 		"\tusing 1:2 notitle with lines";
	txt_spikes.close();
	gpl_spikes.close();
	system("gnuplot spike_num.gpl");
	gplscript << "gnuplot spike_num.gpl" << endl;
	gplscript << "gnuplot spike_raster_plot.gpl" << endl;

	// Create spike raster plot
	createSRP(gpl_srp, rasterplot, t_max);
#endif

// ==============================================================================================================================
	// Create plots for mean firing rate and spatial Gauss sigma over irradiance/current stimulus
	createFROverStimPlot(gpl_mf, dE, E_max, PLOT_TYPE_MEAN);
	gplscript << "gnuplot mean_fr.gpl" << endl;
	createFROverStimPlot(gpl_mf, dE, E_max, PLOT_TYPE_GAUSS_AMP);
	gplscript << "gnuplot gauss_amp_fr.gpl" << endl;
	createFROverStimPlot(gpl_sf, dE, E_max, PLOT_TYPE_GAUSS_SIGMA);
	gplscript << "gnuplot gauss_sigma_fr.gpl" << endl;

	// Close files that have remained open
	txt_msf.close();
	gplscript.close(); // script is not executed, it is only created in case it is needed later on

	// Call GNUplot to process the files
	system("gnuplot mean_fr.gpl");
	system("gnuplot gauss_sigma_fr.gpl");

// ==============================================================================================================================
	// Display final information and save parameters
	int tsec = timeMeasure(false);
	char el_time [32];
	sprintf(el_time, "Elapsed time: %d min %d sec", int(floor(tsec / 60)), int(tsec % 60));
	cout << el_time << endl;
	saveParams(el_time); // additional informaton: elapsed time, threshold frequency for constant illumination
	if (contour != NULL)
		*contour << endl; // attach empty line to contour plot file
	if (chdir("../..") == -1)
		showChDirErrMessage();
	cout << separator << endl;

// ==============================================================================================================================
	// Free allocated memory
	delete[] txt_fr;
	delete[] txt_fr_inh;
	delete[] gpl_lfr;
	delete[] gpl_lfr_inh;
	delete[] gpl_fit_light;
#ifdef STIMULUS_COMPARISON
	delete[] gpl_cfr;
	delete[] gpl_cfr_inh;
	delete[] gpl_fit_current;
#endif
	logf.close();
	return 0;
}

#ifdef SEEK_I_CONST
/*** setSeekICVar ***
 * Sets the variable to communicate with NetworkBatch class for seeking I_const */
void setSeekICVar(double *_seekic)
{
	seekic = _seekic;
}
#endif

/*** setParams ***
 * Sets the simulation parameters on given values and resets network(s) */
void setParams(double _t_pulse, double _frequency, double _const_current, double _tau_syn, double _J,
               double _J_ee, double _J_ei, double _J_ie, double _J_ii, double _pc, double _spatial_sigma)
{
	t_pulse = _t_pulse;
	frequency = _frequency;
	pc = _pc;
	tau_syn = _tau_syn;
	J = _J;
	net.setConstCurrent(_const_current);
	net.setSynTimeConstant(_tau_syn);
	net.setCouplingStrengths(_J_ee*J, _J_ei*J, _J_ie*J, _J_ii*J);
	net.setConnections(_pc);
	net.setGaussSigma(_spatial_sigma);
	net.reset();
#ifdef STIMULUS_COMPARISON
	net2.setConstCurrent(_const_current);
	net2.setSynTimeConstant(_tau_syn);
	net2.setCouplingStrengths(_J_ee*J, _J_ei*J, _J_ie*J, _J_ii*J);
	net2.setConnections(_pc);
	net2.setGaussSigma(_spatial_sigma);
	net2.reset();
#endif
}

/*** Constructor ***
 * Sets all parameters on given values and calls constructors for Neuron instances *
 * - _contour: ofstream of contour plot data file or NULL */
NetworkSimulation(const int _Nl, const int _Nl_inh, const double _dt, const double _dE, const double _raster_bin, double _t_offset, double _t_max,
						double _E_max, const vector<double> _stimplot, const vector<int> _rasterplot, ofstream *_contour)
	: offset(int(ceil(_t_offset / dt))), Nl(_Nl), Nl_inh(_Nl_inh), net(_dt, offset, _Nl, _Nl_inh), dt(_dt), dE(_dE), raster_bin(_raster_bin),
          center(double(_Nl)/2.+0.5), center_inh(double(_Nl_inh)/2.+0.5), stimplot(_stimplot), rasterplot(_rasterplot)
#ifdef STIMULUS_COMPARISON
	  , net2(_dt, offset, _Nl, _Nl_inh)
#endif
{
	contour = _contour;
	t_max = _t_max;
	E_max = _E_max;
}

};
