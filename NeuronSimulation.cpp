/**************************************************************************
 * Simulation of a single neuron with Channelrhodopsin-2 channels
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

#include "Tools.hpp"
#include "Plots.hpp"
#include "Network.hpp"

//#define STIMULUS_COMPARISON // performs a comparison between pulsed light stimulus and pulsed current stimulus
#define PHOTOCURRENT_SATURATION // creates a plot showing the saturation of the photocurrent for increasing irradiance (STIMULUS_COMPARISON needs to be set)
#define PLOT_FR_OPROB // creates a plot of instantaneous firing rate and open probability over time for the amplitudes specified in stimplot
//#define MARK_MAXIMA // marks the maxima in instantaneous firing rate plots
//#define MARK_HALFLIFE  // marks the half-life values in instantaneous firing rate plots
#define TAU_ADAPT_FWHM	// seeks tau_adapt, maximum, baseline and FWHM values
#define PLOT_AV_PEAKS  // plots the averaged peaks and visualized FWHM of the firing rate and the open probability for all amplitudes (needs def. of TAU_ADAPT_FWHM)

#define SLIDING_WINDOW // use sliding window to compute instantaneous firing rates (mind the correlations between time bins)
//#define FIXED_WINDOW // use fixed windows to compute instantaneous firing rates (mind that the firing rates are dependent on the window location)

#define fr_window_size 1000 // size of the numerical window to compute the instantaneous firing rate, in units of time bins (must be even)
#define half_value_deviation 0.04 // maximum allowed relative deviation from half the maximum value for finding half-life values
#define sat_tolerance 0.0001  // tolerance for finding saturation value (the greater, the more maxima are considered to be within saturation)
#define num_trials 900 // ATTENTION! -> number of trials to pass through for averaging over the firing rate (srqt(num_trials) HAS TO BE a natural number)
#define sat_fit_limit 150000 // range in time steps in which maxima are considered beyond saturation for tau_adapt fitting
#define freq_max_skip 15.0 // frequency in Hz below which the first maximum will NOT be skipped for tau_adapt fit
#define light_max_skip 7.5 // irradiance in mW/mm² above which the first maximum will NOT be skipped for tau_adapt fit

/*** NeuronSimulation class ***
 * simulates one neuron, has an instances of Neuron class inside a Network object */
class NeuronSimulation {

private:

/*** Simulation parameters ***/
const double dt; // ms, one time step for numerical integration
const double dE; // mW/mm^2, irradiance accuracy (in case of STIMULUS_COMPARISON: also dimension nA)
const int dE_average; // number of irradiance bins to be averaged over (1: no averaging) -> save only every [dE_average] value
const int offset; // timesteps of the temporal offset up to the stimulus onset - is not included in t_max
double t_max;  // ms, total duration of simulation
double E_max; // mW/mm^2, irradiance to go up to (in case of STIMULUS_COMPARISON: also dimension nA)
double t_pulse; // ms, duration of one stimulus pulse
double frequency; // Hz, frequency for stimulus pulses
Network neur; // instance of one network of unconnected neurons
#ifdef STIMULUS_COMPARISON
Network neur2; // instance of another network of unconnected neurons
#endif
string cplotfile; // filename of color plot data file (is needed to have time stamp from first call of simulate())

/*** Output parameters ***/
const vector<double> stimplot; // stimulus amplitudes for which data and plots will be created (only multiples of dE with precision up to .2 are processed)
const double plot_start; // time from which the plots for a specific stimulus amplitude shall be rendered
const double plot_end; // time up to which the plots for a specific stimulus amplitude shall be rendered - mind that firing rates can only be plotted
							  // up to t_max+t_offset - fr_window_size/2

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
	ofstream f ("params.txt");
	f << dateStr("") << endl; // time stamp
	f << endl;

	// Parameters
	f << "Simulation parameters:" << endl;
	f << "dt = " << dt << " ms" << endl;
	f << "dE = " << dE << " mW/mm^2" << endl;
	f << "dE_average = " << dE_average << endl;
	f << "t_max = " << t_max << " ms (" << int(ceil(t_max / dt)) << " steps)" << endl;
	f << "t_offset = " << offset*dt << " ms (" << offset << " steps)" << endl;
	f << "E_max = " << E_max << " mW/mm^2 (" << int(ceil(E_max / dE)) << " steps)" << endl;
	f << "t_pulse = " << t_pulse << " ms" << endl;
	f << "frequency = " << frequency << " Hz" << endl;
	f << "fr_window_size = " << fr_window_size << " (" << (fr_window_size*dt) << " ms, "
#ifdef SLIDING_WINDOW
	  << "sliding)" << endl;
#elif defined FIXED_WINDOW
	  << "fixed)" << endl;
#endif
	f << "half_value_deviation = " << half_value_deviation << endl;
	f << "sat_tolerance = " << sat_tolerance << endl;
	f << "sat_fit_limit = " << sat_fit_limit <<  " (" << (sat_fit_limit*dt) << " ms)" << endl;
	f << "freq_max_skip = " << freq_max_skip << " Hz" << endl;
	f << "light_max_skip = " << light_max_skip << " mW/mm^2" << endl;
	f << "num_trials = " << num_trials << endl;
	neur.saveNeuronParams(&f); // only save the neuron parameters, since no real network is used here
	f << endl;

	// Additional information
	f << "Constant stimulation beyond " << 1000.0/t_pulse << " Hz" << endl;
	f << str << endl;
	f.close();
}

#ifdef TAU_ADAPT_FWHM
/*** getPeakFWHM ***
 * Retrieves the FWHM of an averaged peak that is saved in a data file (in fact, either of the averaged firing rate *
 * peak or the averaged open prob. peak). The FWHM is obtained looking for the two values that are *
 * closest to half of the height minus the baseline, the positions of these values are returned *
 * in hm_left and hm_right - in case of ambiguities, always chooses the most left value for hm_left *
 * and the most right value for hm_right *
 * - peakfile: file containing the peak data (text file, first column must be time)
 * - col: column of the file in which the y values to be used are located
 * - amp: the peak's absolute amplitude and its time step position
 * - baseline: y value baseline (and the time position where it occurred), everything above is investigated
 * - hm_left: the variable to return the left "half maximum" in units of time steps
 * - hm_right: the variable to return the right "half maximum" in units of time steps
 * - return: true if determining the values was successful */
bool getPeakFWHM(string peakfile, int col, idpair amp, idpair baseline, int &hm_left, int &hm_right)
{
	ifstream pf(peakfile.c_str(), ios::in);
	string buf; // buffer to read one line
	const char *sep_chars = {"\t "}; // possible characters to separate columns
	double max = (amp.y-baseline.y);  // the maximum of the peak relative to the baseline
	double last_y = 0., last_last_y = 0.; // the two last y values
	int line = 0; // line in file (equivalent to time step)

	hm_left = -1;
	hm_right = -1;

	if (!pf.is_open()) // check if file was opened successfully
		return false;

	while (getline(pf, buf)) // while end of file has not yet been reached
	{
		int pos_start, pos_end; // start and end position of a column in the buffer
		double y; // y value

		pos_end = buf.find_first_of(sep_chars);
		//double t = stod(buf.substr(0, pos_end)); // read time value (first column)

		for (int i=1; i<=col; i++) // find the wanted column
		{
			if ( (pos_start = buf.find_first_not_of(sep_chars, pos_end+1)) == -1 ) // find start of column
			{
				hm_left = 0; hm_right = 0;
				return false; // column is missing
			}

			if ( (pos_end = buf.find_first_of(sep_chars, pos_start+1)) == -1 ) // find end of column
				pos_end = buf.length(); // this was the last column in this line
		}

		y = stod(buf.substr(pos_start, pos_end-pos_start)); // read the y value from the wanted column - one place of precision gets lost!
		y -= baseline.y; // subtract the baseline

		if (( line >= 2 ? isCloser(last_y, max/2., last_last_y) : true ) // is closer to half than previous value?
				 		 && ( isCloser(last_y, max/2., y) )) // is closer to half than next value?
		{
			// only accept values with a deviation of less than half_value_deviation to half the maximum value
			if ( abs(last_y - max/2.) / (max/2.) < half_value_deviation ) {

				if (baseline.x > amp.x) // possible order of occurence: LMRB, MRBL (L = hm_left, M = maximum, R = hm_right, B = baseline)
				{
					if (hm_left < 0 && line <= amp.x) // search as far left from maximum as possible (LMRB)
						hm_left = line-1;
					else if (hm_left < 0 && line > baseline.x) // search right from baseline again for hm_left (MRBL)
						hm_left = line-1;
					else if (line > amp.x && line <= baseline.x) // hm_right lies in between maximum and baseline, go as far right as possible
						hm_right = line-1;
				}
				else // possible order of occurence: BLMR, RBLM
				{
					if (hm_left < 0 && line > baseline.x && line <= amp.x) // hm_left lies in between baseline and maximum, stay as far left as possible
						hm_left = line-1;
					else if (line > amp.x) // search as far right from maximum as possible (BLMR)
						hm_right = line-1;
					else if (line <= baseline.x) //  search left from baseline again for hm_right (RBLM)
						hm_right = line-1;
				}

			}
		}
		last_last_y = last_y;
		last_y = y;
		line++;
	}
	pf.close();

	if (hm_left >= 0 && hm_right >= 0) // if both parts were found
		return true;

	hm_left = 0; hm_right = 0;
	return false;
}

/*** getTauEstimate ***
 * Computes the adaptation time in an averaging way (is used as initial value for GSL fit or *
 * even as final value in case the fit does not yield a proper result) *
 * - maxima: vector of maxima idpairs of the considered quantity that shall be used for computation of tau_adapt
 * - fit_max: the idpair of the fit maximum of the considered quantity
 * - Q_steady: the steady-state y-value of the considered quantity Q (in fact either firing rate or OSP)
 * - return: the tau_adapt value, -1. in case of failure */
double getTauEstimate(vector<idpair> &maxima, idpair &fit_max, double &Q_steady)
{
	double Q_zero = fit_max.y - Q_steady; // for averaging over tau_adapt; equals Q(t=0) - Q_steady
	int num_tau = 0; // count of tau's used for computation
	double mean_tau = 0.; // mean tau value (averaged over the computed tau's of all fit points)

	if (maxima.size() > 0)
	{
		for (int j=0; j<maxima.size(); j++)
		{
			// compute tau_adapt:
			if (  (maxima.at(j).y > Q_steady) 	 // do neither consider those maxima for which (O(t) - Q_steady) becomes <= 0
				&& (maxima.at(j).y < fit_max.y) ) // nor values greater or equal to the fit maximum (otherwise, log becomes 0)
															 // (all other unacceptable values must have been excluded earlier from the maxima vector)
			{
				// Formula: tau_adapt = (t - t_0) / ln(Q_zero / (Q(t) - Q_steady))
				double tau = double(maxima.at(j).x - fit_max.x) * dt / log(Q_zero / (maxima.at(j).y - Q_steady));
				if (!std::isnan(tau)) // nan occurs usually for Q_steady > maxima[j] (first if-clause)
				{
					mean_tau += tau;
					num_tau++;
				}
			}
		}
	}
	else
		return -1.;

	if (num_tau > 0)
		mean_tau /= num_tau; // average over number of tau values
	else
		return -1.;

	return mean_tau;
}
#endif

public:

/*** simulate ***
 * Runs the neuron simulation *
 * - working_dir: working directory *
 * - first_sim: tells if this is the first simulation run by the batch code */
int simulate(string working_dir, bool first_sim)
{
	// Start time measurement
	timeMeasure(true);

	// Constants
	const int m = int(ceil(E_max / dE)); // number of irradiance steps
	const int n = int(ceil(t_max / dt)); // number of time steps
	const int pulse_length = int(ceil(t_pulse / dt)); // number of time steps of one pulse
	const int period_steps = int(floor((1000/frequency)/dt)); // number of time steps of one period
	const int period_num = int(floor(t_max/(1000/frequency))); // number of full periods in simulation runtime

	const string path = working_dir + "/f=" + dtos(frequency, 2); // path to working directory
	const string separator = getSeparator(); cout << separator << endl; // string of characters for a separator in command line

	// Try to change directory
	mkDir(path);
	if (chdir(path.c_str()) == -1) {
		showChDirErrMessage();
		cout << separator << endl;
		return -1;
	}

	// Output with general information
	cout << "\x1b[33mSimulating neuron for " << t_max << " ms (" << dateStr("", !first_sim) << ")\x1b[0m" << endl; // new time stamp here, possibly
	cout << "Stimulus: rectangle, t_pulse = " << t_pulse << " ms, \x1b[34mfrequency = " << frequency << " Hz\x1b[0m" << endl;

	// Declarations and initializations
	double *spikes_lst = new double[m+1]; // average spikes per trial neuron evoked by light stimulus
#ifdef STIMULUS_COMPARISON
	double *spikes_cst = new double[m+1]; // average spikes per trial neuron evoked by current stimulus
#endif
	idpair global_peak [m+1]; // array of the global total current peaks for each stimulus amplitude (for light stimulation), are initially set to zero by constructor
	idpair global_min [m+1]; // array of the global total current minima for each stimulus amplitude (for light stimulation), are initially set to zero by constructor
#ifdef STIMULUS_COMPARISON
	idpair global_peak2 [m+1]; // array of the global total current peaks for each stimulus amplitude (for current stimulation)
	idpair global_min2 [m+1]; // array of the global total current minima for each stimulus amplitude (for current stimulation)
#endif
	idpair global_peak_IChR2 [m+1]; // array of the global ChR2 current peaks for each stimulus amplitude (for light stimulation)
	idpair global_fr_max[m+1]; // contains the highest value of the firing rate temporal course (is initially set to zero in loop)
	idpair global_O_max[m+1]; // contains the highest value of the open prob. temporal course (is initially set to zero in loop)

	ofstream of(working_dir + "/" + cplotfile, ofstream::app); // color plot data file for appending (was previously solely written by gnuplot)
	ofstream data(dateStr("_data.txt")); // data file for creating plots over stimulus amplitude
	ofstream gpl_if("firingrate.gpl"); // gnuplot file for creating a plot of firing rate over irradiance/current amplitude in a PDF file
	ofstream* txt_stim = new ofstream[stimplot.size()]; // data for creating plots for specific stimulus amplitudes
	ofstream* gpl_stim = new ofstream[2*stimplot.size()]; // gnuplot files for creating current (1st half) / voltage (2nd half) plots for specific light stimulus 																				// amplitudes
#ifdef STIMULUS_COMPARISON
	ofstream* gpl_stim2 = new ofstream[2*stimplot.size()]; // gnuplot files for creating current (1st half) / voltage (2nd half) plots for specific
																			 // light stimulus amplitudes
#ifdef PHOTOCURRENT_SATURATION
	ofstream gpl_ichr2("light_current.gpl"); // gnuplot file for creating a plot of light current over light stimulus amplitude in a PDF file
#endif
#endif
	Stimulus lst = Stimulus(period_steps); // light stimulus
#ifdef STIMULUS_COMPARISON
	Stimulus cst = Stimulus(period_steps); // current stimulus for neur2
#endif
	double p = 0.; // percentage of process completion
	int current_E = 0; // specifies current position (index) in stimplot[] array
	double fr; // contains the current averaged instantaneous firing rates
#ifdef PLOT_FR_OPROB
	ofstream* txt_fr_oprob = new ofstream[stimplot.size()]; // data for creating firing rate/open probability comparison
	ofstream* gpl_fr_oprob = new ofstream[stimplot.size()]; // gnuplot files for creating firing rate/open probability comparison
#endif
	const double f_threshold = (1000.0/t_pulse) - 1.0; // threshold frequency - constant illumination from this frequency upward
	vector<double> oprob_old(fr_window_size/2, 0.); // vector containing the last <fr_window_size/2> open probabilities (.front() at t, .back() at t')
	vector<double> cprob_old(fr_window_size/2, 0.); // vector containing the last <fr_window_size/2> closed probabilities (.front() at t, .back() at t')
	vector<double> peak_fr(period_steps, 0.); // vector containing all summed up firing rates of one stimulus period ("the average firing rate peak")
	vector<double> peak_O(period_steps, 0.); // vector containing all summed up open probabilities of one stimulus period ("the average O(t) peak")
#ifdef PLOT_AV_PEAKS
	ofstream* gpl_peak_fr = new ofstream[m]; // plot files for average firing rate peaks
	ofstream* gpl_peak_O = new ofstream[m]; // plot files for average open probability peaks
#endif
	ofstream* txt_peak = new ofstream[m]; // data files for average firing rate/open probability peaks
	double fr_steady[m+1]; // firing rate steady-state approximation for each irradiance
	double O_steady[m+1]; // open prob. steady-state approximation for each irradiance
	int sat_begin[m+1]; // marks first time step of saturation for each irradiance (is determined considering the open probabilities, not the firing rates!)
	ofstream gplscript("gpl"); // shell script for calling all the gnuplot scripts (to re-generate the plots, e.g. for another gnuplot version)
	vector<idpair> maxima_fr; // vector that contains all temporal firing rate maxima (as a pair of values)
	vector<idpair> maxima_O; // vector that contains all temporal open prob. maxima (as a pair of values)
	double tau_adapt_fr[m+1], tau_adapt_O[m+1]; // tau_adapt values for each stimulus amplitude
	idpair fit_fr_max[m+1]; // array of fit function maxima for tau_adapt_fr fit for each stimulus frequency
	idpair fit_O_max[m+1]; // array of fit function maxima for tau_adapt_O fit for each stimulus frequency
	ofstream errlog("err.log"); // output text file for error messages

	// Prepare files for specific stimulus plots
	for (int i=0; i<2*stimplot.size(); i++)
	{
		char filename[100];

		if (i < stimplot.size())
		{
			sprintf(filename, "_stim_%.2f.txt", stimplot[i]); // only up to stimplot.size()
			txt_stim[i].open(dateStr(filename)); // only up to stimplot.size()
			sprintf(filename, "lightstimulus_I_%.2f_PDFplot.gpl", stimplot[i]); // current plot
		}
		else
			sprintf(filename, "lightstimulus_V_%.2f_PDFplot.gpl", stimplot[i-stimplot.size()]); // voltage plot
		gpl_stim[i].open(filename);
#ifdef STIMULUS_COMPARISON
		if (i < stimplot.size())
			sprintf(filename, "currentstimulus_I_%.2f_PDFplot.gpl", stimplot[i]); // current plot
		else
			sprintf(filename, "currentstimulus_V_%.2f_PDFplot.gpl", stimplot[i-stimplot.size()]); // voltage plot
		gpl_stim2[i].open(filename);
#endif
		if ( (i < stimplot.size() && !txt_stim[i].is_open()) || !gpl_stim[i].is_open()
#ifdef STIMULUS_COMPARISON
				|| !gpl_stim2[i].is_open()
#endif
			)
		{
			cout << "Unable to open file!" << endl << separator << endl;
			return -1;
		}

#ifdef PLOT_FR_OPROB
		if (i < stimplot.size())
		{
			sprintf(filename, "_fr_oprob_%.2f.txt", stimplot[i]);
			txt_fr_oprob[i].open(dateStr(filename));
			sprintf(filename, "fr_oprob_%.2f_PDFplot.gpl", stimplot[i]);
			gpl_fr_oprob[i].open(filename);
			if (!txt_fr_oprob[i].is_open() || !gpl_fr_oprob[i].is_open())
			{
				cout << "Unable to open file!" << endl << separator << endl;
				return -1;
			}
		}
#endif
	}

	// Check if files have been opened properly
  	if ( !gplscript.is_open() || !of.is_open() || !data.is_open() || !gpl_if.is_open() || !errlog.is_open() )
	{
		cout << "Unable to open file!" << endl << separator << endl;
		return -1;
	}

	// (Light/current) stimulation loop
	for(int i = 0; i <= m; i++)
	{
		idpair peak_height_fr(0,std::numeric_limits<double>::lowest()),
                       peak_height_O(0,std::numeric_limits<double>::lowest()); // maxima of the averaged peaks
		idpair peak_base_fr(0,std::numeric_limits<double>::max()),
		       peak_base_O(0,std::numeric_limits<double>::max()); // baselines of the averaged peaks
		idpair local_max_fr(0,std::numeric_limits<double>::lowest()), local_max_oprob(0,std::numeric_limits<double>::lowest()); // temporary variable for local maxima

		// Reset neuron(s), global current minima and spike count number(s)
		neur.reset();
		spikes_lst[i] = 0.;
		global_min[i] = idpair(0, neur.getConstCurrent(0));
#ifdef STIMULUS_COMPARISON
		neur2.reset();
		spikes_cst[i] = 0.;
		global_min2[i] = idpair(0, neur2.getConstCurrent(0));
#endif

		// Clear stimuli
		lst.clear();
#ifdef STIMULUS_COMPARISON
		cst.clear();
#endif

		// Clear vectors of maxima
		vector<idpair>().swap(maxima_fr);
		vector<idpair>().swap(maxima_O);

		// Clear open and closed probability vectors
		oprob_old.assign(fr_window_size/2, 0.);
		cprob_old.assign(fr_window_size/2, 0.);

		// Clear average firing rate/open probability peak
		peak_fr.assign(period_steps, 0.);
		peak_O.assign(period_steps, 0.);

		if (i > 0)
		{
		// No saturation in the beginning
			sat_begin[i] = 0;
			O_steady[i] = 0.;

#ifdef TAU_ADAPT_FWHM
		// Open average peak files
#ifdef PLOT_AV_PEAKS
			gpl_peak_fr[i-1].open(string("fr_peak_") + dtos(i*dE,2) + string(".gpl"));
			gpl_peak_O[i-1].open(string("O_peak_") + dtos(i*dE,2) + string(".gpl"));
#endif
			txt_peak[i-1].open(dateStr(string("_peak_") + dtos(i*dE,2) + string(".txt")));

			if (
#ifdef PLOT_AV_PEAKS
				 !gpl_peak_fr[i-1].is_open() || !gpl_peak_O[i-1].is_open() ||
#endif
				 !txt_peak[i-1].is_open() )
			{
				cout << "Unable to open file!" << endl << separator << endl;
				return -1;
			}
#endif
		}

		// Set global firing rate/open prob. maximum to zero
		global_fr_max[i] = idpair(0, 0.);
		global_O_max[i] = idpair(0, 0.);

		// Update percentage
		double p_new = double(round(double(i) / double(m) * 1000.0)) / 10.0; // round to first decimal place
		if (p_new > p)
		{
			p = p_new;
			printf("\rProgress: %.1f %% completed.", p);
			fflush(stdout);
		}

		// Time loop
		for (int j = 0; j <= n+offset+fr_window_size/2; j++)
		{
			// Stimuli start at offset time
			if (j == offset)
			{
				// Set light stimulus
				lst.addRectPulse(i*dE, 0, pulse_length);
				//lst.printShape("a", dt);
				neur.setUniformLightStimulus(lst);

				// Set current stimulus
#ifdef STIMULUS_COMPARISON
				cst.addRectPulse(neur2.getSteadyCurrent(i*dE, 0), 0, pulse_length); // use I_s as stimulation current
				neur2.setUniformCurrentStimulus(cst);
#endif
			}

			// Calculate next step for all Neuron trials (including ChR2)
			int num_spikes_lst = neur.processTimeStep(j);
#ifdef STIMULUS_COMPARISON
			int num_spikes_cst = neur2.processTimeStep(j);
#endif

			// Count spikes per trial neuron (only those after initial offset)
			if (j >= offset)
			{
				spikes_lst[i] += double(num_spikes_lst) / num_trials;
#ifdef STIMULUS_COMPARISON
				spikes_cst[i] += double(num_spikes_cst) / num_trials;
#endif
			}

			// Set averaged firing rate to zero (because firing rates are added to this variable)
			fr = 0.;

			// Shift elements of open and closed probability vectors for this time step
			oprob_old.insert(oprob_old.begin(), 0.); // insert 0. at first position in vector
			oprob_old.erase(oprob_old.end()-1); // delete last element of the vector
			cprob_old.insert(cprob_old.begin(), 0.); // insert 0. at first position in vector
			cprob_old.erase(cprob_old.end()-1); // delete last element of the vector

			// Loop over trials (to determine instantaneous firing rates)
			for (int t = 0; t < num_trials; t++)
			{
				// Update open and closed probability vectors
				oprob_old[0] += neur.getOpenProb(t); // add open prob. of trial t to open probability vector
				cprob_old[0] += neur.getClosedProb(t); // add closed prob. of trial t to closed probability vector

#ifdef SLIDING_WINDOW
				// Compute instantaneous firing rate using a sliding window
				// 	(what is actually computed is the firing rate fr_window_size/2 in the past)
				if (j >= fr_window_size/2) 			// start when half the window time has passed
				{

#elif defined FIXED_WINDOW
				// Compute instantaneous firing rate using fixed windows
				// 	(what is actually computed is the spike-count rate for the last fr_window_size)
				if ( (j >= fr_window_size)	// start if at least one window time has passed
				  && (j % fr_window_size == 0) ) // only if a full window is completed
				{

#endif
					int s = neur.getSpikeCount(t) + 1; // the number of the currently considered spike, starting from the one that occurred latest
					int fws_eff; // the effective window size (only of importance for the beginning)

					if (j < fr_window_size)
						fws_eff = j; // reduced window size (in the very first time steps)
					else
						fws_eff = fr_window_size; // normal window size (most of the time)

					int time_s, time_to_s = 0;

					do // determine number s of first spike within the window
					{
						time_s = neur.getSpikeTime(--s, t); // get time of the s-th spike of trial t
						time_to_s = abs(j - time_s); // get time that has passed since the s-th spike of trial t
					} while(time_s >= 0 && time_to_s < fws_eff); // loop as long as spikes are found and these spikes are within the window

					fr += (neur.getSpikeCount(t) - s) / (fr_window_size*dt/1000); // compute the instantaneous firing rate - important: division by fr_window_size!
				}

			} // End of for(t)

			// Divide by number of trials to obtain averaged open and closed probability
			oprob_old[0] /= num_trials;
			cprob_old[0] /= num_trials;

			// Divide by number of trials to obtain averaged instantaneous firing rate
			fr /= num_trials;


			if (i > 0) // for i=0, there is no illumination and thus no saturation/peaks
			{
			// Find global maxima, collect local maxima, detect saturation and compute average steady-state peak
				if ( ((j-offset) >= fr_window_size/2)  // start when offset is over and half the window time has passed
				  && ((j-fr_window_size/2) < (period_num*period_steps + offset)) ) // do this as long as t' is still within a FULL period
				{

					int index = (j-offset-fr_window_size/2) % period_steps; // time bin in stimulus period related to t' time bin

					if ( sat_begin[i] ) // in saturation: compute average steady-state peak
					{
						peak_fr[index] += fr; // sum up stimulus periods: add firing rate to related stimulus period time bin
						peak_O[index] += oprob_old.back(); // sum up stimulus periods: add open probability to related stimulus period time bin
					}

					if (fr > global_fr_max[i].y) // find global firing rate maximum
						global_fr_max[i] = idpair(j-fr_window_size/2, fr);

					if (oprob_old.back() > global_O_max[i].y) // find global open probability maximum
						global_O_max[i] = idpair(j-fr_window_size/2, oprob_old.back());


					if (index == 0) // Beginning of current period
					{
						local_max_oprob.y = -1.;
						local_max_fr.y = -1.;
					}
					if (oprob_old.back() > local_max_oprob.y) // find local open probability maximum (ignores step j = n; within tolerance "negligible")
						local_max_oprob = idpair(j-fr_window_size/2, oprob_old.back());
					if (fr > local_max_fr.y)	// find local firing rate maximum
						local_max_fr = idpair(j-fr_window_size/2, fr);

					// End of current period, save maxima, detect saturation
					// --> constant illumination is covered as well, for there is one maximum chosen in every period!
					if (index == period_steps-1)
					{
						maxima_fr.push_back(local_max_fr);
						maxima_O.push_back(local_max_oprob);

						// Detect saturation (only done for O(t), firing rate is assumed to be in steady state as well) - rough estimate of O_steady is used,
						// preciser values for O_steady and fr_steady are obtained from averaged peaks
						if ( !sat_begin[i] )
						{
							if ( abs(O_steady[i] - local_max_oprob.y) < sat_tolerance ) // is last maximum value close enough to this one?
								sat_begin[i] = j - fr_window_size/2; // saturation "starts" at beginning of next stimulus period (in terms of t' !!!)
							else
								O_steady[i] = local_max_oprob.y; // use this maximum value for further checks
						}
					}

					// Last passage of a FULL stimulus period: print firing rate peak in file
					if ( (period_num*period_steps + offset - period_steps) <= (j-fr_window_size/2) ) // check that t' is beyond the beginning of the last FULL
									  																							// period
					{
						int div = (period_num*period_steps + offset - sat_begin[i]) / period_steps; // division always yields an integer number
						peak_fr[index] /= div; // divide by number of full periods within saturation
						peak_O[index] /= div; // divide by number of full periods within saturation

						if (peak_fr[index] > peak_height_fr.y) // determine firing rate peak maximum
							peak_height_fr = idpair(index, peak_fr[index]);

						if (peak_O[index] > peak_height_O.y) // determine open prob. peak maximum
							peak_height_O = idpair(index, peak_O[index]);

						if (peak_fr[index] < peak_base_fr.y) // determine firing rate peak minimum
							peak_base_fr = idpair(index, peak_fr[index]);

						if (peak_O[index] < peak_base_O.y) // determine open prob. peak minimum
							peak_base_O = idpair(index, peak_O[index]);

#ifdef FIXED_WINDOW
						if ((j-offset) % fr_window_size == 0) // to make sure an output is written only every fr_window_size step
#endif
							txt_peak[i-1] << fixed << index*dt << "\t\t\t" << peak_fr[index] << "\t\t\t" << peak_O[index] << endl; // print into file
					}
				}
			}

			// Plots for specific stimulus amplitudes
			if (considerAmplitude(i)) {
				double lSt = 0.0;
#ifdef STIMULUS_COMPARISON
				double cSt = 0.0;
#endif
				if (j >= offset) // before offset, there is no stimulus
				{
					lSt = lst.getStimulusAt(j-offset);
#ifdef STIMULUS_COMPARISON
					cSt = cst.getStimulusAt(j-offset);
#endif
				}
#ifdef STIMULUS_COMPARISON
				txt_stim[current_E] << fixed << j*dt << "\t\t\t" // writing to a file causes slow-down, but few memory is used
										  << neur.getCurrent(0) << "\t\t\t" << neur.getVoltage(0) << "\t\t\t" << lSt << "\t\t\t" << neur.getChR2Current(0) << "\t\t\t"
										  << neur2.getCurrent(0) << "\t\t\t" << neur2.getVoltage(0) << "\t\t\t" << cSt << endl;
#else
				txt_stim[current_E] << fixed << j*dt << "\t\t\t" // writing to a file causes slow-down, but few memory is used
										  << neur.getCurrent(0) << "\t\t\t" << neur.getVoltage(0) << "\t\t\t" << lSt << endl;
#endif
#ifdef PLOT_FR_OPROB
				if (j >= fr_window_size/2) // write the instantaneous firing rate and the open and closed state probabilities (half the window size ago, at t') in the file
					txt_fr_oprob[current_E] << fixed << (j-fr_window_size/2)*dt << "\t\t\t" << fr << "\t\t\t" << oprob_old.back() << "\t\t\t" << cprob_old.back() << endl;
					// writing to a file causes slow-down, but few memory is used
#endif
				if (j==n+offset+fr_window_size/2)
					current_E++;
			}

			// Seek global total current peak for light stimulation
			if (global_peak[i].y < neur.getCurrent(0))
			{
				global_peak[i] = idpair(j, neur.getCurrent(0));
			}
			// Seek global total current minimum for light stimulation
			else if (global_min[i].y > neur.getCurrent(0))
			{
				global_min[i] = idpair(j, neur.getCurrent(0));
			}

#ifdef STIMULUS_COMPARISON
			// Seek global total current peak for current stimulation
			if (global_peak2[i].y < neur2.getCurrent(0))
			{
				global_peak2[i] = idpair(j, neur2.getCurrent(0));
			}
			// Seek global total current minimum for current stimulation
			else if (global_min2[i].y > neur2.getCurrent(0))
			{
				global_min2[i] = idpair(j, neur2.getCurrent(0));
			}
#endif
			// Seek global ChR2 current peak for light stimulation
			if (global_peak_IChR2[i].y < neur.getChR2Current(0))
			{
				global_peak_IChR2[i] = idpair(j, neur.getChR2Current(0));
			}


		} // End of for(j)

		// Close average peak file and find FWHM, mean peak height, baseline
		int hm_left_fr, hm_right_fr, hm_left_O, hm_right_O;
#ifdef TAU_ADAPT_FWHM
		if (i > 0)
		{
			txt_peak[i-1].close();

			if (frequency < f_threshold) // pulsed illumination
			{
				if (!getPeakFWHM(dateStr(string("_peak_") + dtos(i*dE,2) + string(".txt")), 1, peak_height_fr, // get FWHM of averaged firing rate peak
								     peak_base_fr, hm_left_fr, hm_right_fr))
					errlog << "Error while trying to find FWHM of firing rate peak at irradiance " << i*dE << endl;
				fr_steady[i] = peak_height_fr.y; // use averaged height as new steady state value
			}
			else // constant illumination
			{
				fr_steady[i] = (peak_height_fr.y+peak_base_fr.y) / 2; // use mean between averaged height and baseline (should be very close)
																					   // as new steady state value
				hm_left_fr = 0; hm_right_fr = 0; // there are no peaks, thus no FWHM
			}

#ifdef PLOT_AV_PEAKS
			createPeakPlot(gpl_peak_fr[i-1], period_steps*dt, i*dE, PLOT_TYPE_FR, peak_height_fr.y, // create plot of the firing rate peak with visualized FWHM
								peak_height_fr.x*dt, peak_base_fr.y, hm_left_fr*dt, hm_right_fr*dt);
			gplscript << "gnuplot fr_peak_" << dtos(i*dE,2) << ".gpl" << endl;
#endif

			if (frequency < f_threshold) // pulsed illumination
			{
				if (!getPeakFWHM(dateStr(string("_peak_") + dtos(i*dE,2) + string(".txt")), 2, peak_height_O, // get FWHM of averaged open prob. peak
								     peak_base_O, hm_left_O, hm_right_O))
					errlog << "Error while trying to find FWHM of open prob. peak at irradiance " << i*dE << endl;
				O_steady[i] = peak_height_O.y; // use averaged height as new steady-state value
			}
			else // constant illumination
			{
				O_steady[i] = (peak_height_O.y+peak_base_O.y) / 2; // use mean between averaged height and baseline (should be the same, though)
																				   // as new steady-state value
				hm_left_O = 0; hm_right_O = 0; // there are no peaks, thus no FWHM
			}

#ifdef PLOT_AV_PEAKS
			createPeakPlot(gpl_peak_O[i-1], period_steps*dt, i*dE, PLOT_TYPE_O, peak_height_O.y, // create plot of the open prob. peak with visualized FWHM
								peak_height_O.x*dt, peak_base_O.y, hm_left_O*dt, hm_right_O*dt);
			gplscript << "gnuplot O_peak_" << dtos(i*dE,2) << ".gpl" << endl;
#endif
		}
#endif

		double err_fr[4]; // FR fit parameter error array
		double err_O[4]; // O(t) fit parameter error array

#ifdef TAU_ADAPT_FWHM
		// Fit tau_adapt
		if (i > 0)
		{
			if (maxima_fr.size() != period_num)
				errlog << "Number of fr. maxima (" << maxima_fr.size() << ") and number of periods (" << period_num << ") do not match!" << endl;
			if (maxima_O.size() != period_num)
				errlog << "Number of O(t) maxima (" << maxima_O.size() << ") and number of periods (" << period_num << ") do not match!" << endl;

			// Determine fit maximum
			fit_fr_max[i] = global_fr_max[i]; // by default, fit function maximum is the global maximum
			fit_O_max[i] = global_O_max[i]; // by default, fit function maximum is the global maximum
			if ( (frequency >= freq_max_skip) && (i*dE <= light_max_skip) ) // frequency (or rather number of maxima) beyond certain threshold and irradiance below
																								 // another threshold - skip maxima
			{
				fit_fr_max[i] = firstIdpair(global_fr_max[i].x, maxima_fr, 1);
				fit_O_max[i] = firstIdpair(global_O_max[i].x, maxima_O, 1);
			}

			vector<double> yerr_for_fit; // data points error vector
			int iter; // for return value of findFit

			// erase maxima that shall not be used for fitting
			for (int mint = 0; mint < maxima_fr.size(); mint++)
			{
				if (maxima_fr.at(mint).x < fit_fr_max[i].x) // do not use maxima before fit maximum
					maxima_fr.erase(maxima_fr.begin()+mint--);
				else if (maxima_fr.at(mint).x > sat_begin[i]+sat_fit_limit) // do not use maxima that lie more than sat_fit_limit beyond saturation
					maxima_fr.erase(maxima_fr.begin()+mint--);
			}
			for (int mint = 0; mint < maxima_O.size(); mint++)
			{
				if (maxima_O.at(mint).x < fit_O_max[i].x) // do not use maxima before fit maximum
					maxima_O.erase(maxima_O.begin()+mint--);
				else if (maxima_O.at(mint).x > sat_begin[i]+sat_fit_limit) // do not use maxima that lie more than sat_fit_limit beyond saturation
					maxima_O.erase(maxima_O.begin()+mint--);
			}

			// set the properties of the fit function
			gsl_multifit_function_fdf func; // fit function plus its Jacobian (fdf)
			func.f = &exp_f; // the fit function
			func.p = 4; // number of parameters (both fit parameters and constant parameters)

			// fit tau_adapt_fr from firing rate maxima
			tau_adapt_fr[i] = getTauEstimate(maxima_fr, fit_fr_max[i], fr_steady[i]); // get estimate for tau_adapt using old method (cf. ChannelSimulation.cpp)
			if (maxima_fr.size() > func.p)
			{
				double par_fr[4] = {fit_fr_max[i].y-fr_steady[i], // A (var)
										  tau_adapt_fr[i]/dt, // tau (var)
										  fr_steady[i], // b (fixed)
										  double(fit_fr_max[i].x)}; // t_0 (fixed) // parameter array for fit, containing initial values
				vector<double>(maxima_fr.size(), 0.1).swap(yerr_for_fit); // fill yerr_for_fit with as many 0.1's as the size of maxima_fr
				fitdata_idpairs fd_fr = {maxima_fr.size(), &maxima_fr[0], &yerr_for_fit[0]}; // fit data for fitting firing rate maxima
				func.df = &exp_df_A_tau; // the Jacobian of the fit function (fit via A and tau)
				func.fdf = &exp_fdf_A_tau; // function calling _f and _df functions (fit via A and tau)
				func.n = maxima_fr.size(); // number of data points
				func.params = &fd_fr; // fitdata structure (contains arrays for (x, y, yerror))
				iter = findFit(func, par_fr, err_fr, NULL); // do the fit procedure
				if (iter > 0)
				{
					fit_fr_max[i].y = par_fr[0] + par_fr[2];

					errlog << "FR fit converged after " << iter << " iterations (tau = " << par_fr[1]*dt << "+-" << err_fr[1]*dt
								 << ", initial tau = " << tau_adapt_fr[i] << ", A = " << par_fr[0] << "+-" << err_fr[0] << ", b = " << par_fr[2] << "+-"
								 << err_fr[2] << ")!" << std::endl;
					tau_adapt_fr[i] = par_fr[1] * dt;
				}
				else
				{
					errlog << "FR fit did not converge (initial tau = " << tau_adapt_fr[i] << ")!" << endl;
					tau_adapt_fr[i] = 0.;
				}
			}
			else
			{
				errlog << "No FR fit possible (initial tau = " << tau_adapt_fr[i] << ")!" << endl;
				tau_adapt_fr[i] = 0.;
			}

			// fit tau_adapt_O from open prob. maxima
			tau_adapt_O[i] = getTauEstimate(maxima_O, fit_O_max[i], O_steady[i]); // get estimate for tau_adapt using old method (cf. ChannelSimulation.cpp)
			if (maxima_O.size() > func.p)
			{
				double par_O[4] = {fit_O_max[i].y-O_steady[i], // A (var)
										 tau_adapt_O[i]/dt, // tau (var)
										 O_steady[i], // b (fixed)
										 double(fit_O_max[i].x)}; // t_0 (fixed) // parameter array for fit, containing initial values
				vector<double>(maxima_O.size(), 0.02).swap(yerr_for_fit); // fill yerr_for_fit with as many 0.02's as the size of maxima_O
				fitdata_idpairs fd_O = {maxima_O.size(), &maxima_O[0], &yerr_for_fit[0]}; // fit data for fitting open prob. maxima
				func.df = &exp_df_A_tau; // the Jacobian of the fit function (fit via A and tau)
				func.fdf = &exp_fdf_A_tau; // function calling _f and _df functions (fit via A and tau)
				func.n = maxima_O.size(); // number of data points
				func.params = &fd_O; // fitdata structure (contains arrays for (x, y, yerror))
				iter = findFit(func, par_O, err_O, NULL); // do the fit procedure
				if (iter > 0)
				{
					fit_O_max[i].y = par_O[0] + O_steady[i];
					errlog << "O(t) fit converged after " << iter << " iterations (tau = " << par_O[1]*dt << "+-" << err_O[1]*dt
								 << ", initial tau = " << tau_adapt_O[i] << ", A = " << par_O[0] << "+-" << err_O[0] << ")!" << std::endl;
					tau_adapt_O[i] = par_O[1] * dt;
				}
				else
				{
					errlog << "O(t) fit did not converge (initial tau = " << tau_adapt_O[i] << ")!" << endl;
					tau_adapt_O[i] = 0.;
				}
			}
			else
			{
				errlog << "No O(t) fit possible (initial tau = " << tau_adapt_O[i] << ")!" << endl;
				tau_adapt_O[i] = 0.;
			}
		} // end of tau_adapt fits
#endif

		// Average over certain number of irradiance bins (no averaging for dE_average = 1)
		if ((i % dE_average) == 0)
		{
			double fr_irradiance = 0.0;
			int fri_div = dE_average;
#ifdef STIMULUS_COMPARISON
			double fr_current = 0.0;
			int frc_div = dE_average;
#endif

			// Sum up
			for (int r=0; r < dE_average; r++)
			{
				if (!std::isnan(fr_irradiance) && spikes_lst[i-r] > EPSILON)
					fr_irradiance += spikes_lst[i-r];
				else
					fri_div--;
#ifdef STIMULUS_COMPARISON
				if (!std::isnan(fr_current) && spikes_cst[i-r] > EPSILON)
 					fr_current += spikes_cst[i-r];
				else
					frc_div--;
#endif
			}

			fr_irradiance = ((fri_div > 0) ? fr_irradiance / double(fri_div) : 0.0) / (t_max / 1000);
#ifdef STIMULUS_COMPARISON
			fr_current = ((frc_div > 0) ? fr_current / double(frc_div) : 0.0) / (t_max / 1000);
#endif

			// Write data file (dtos(,,true) ensures that no zeros are printed)
			data << fixed << double(i)*dE
				  << "\t\t\t" << dtos(fr_irradiance, 6, true);
#ifdef STIMULUS_COMPARISON
			data << "\t\t\t" << neur2.getSteadyCurrent(i*dE, 0)
				  << "\t\t\t" << dtos(fr_current, 6, true);
#endif
			data << "\t\t\t" << dtos(tau_adapt_fr[i], 6, true)
				  << "\t\t\t" << err_fr[1]*dt
				  << "\t\t\t" << dtos(tau_adapt_O[i], 6, true)
				  << "\t\t\t" << err_O[1]*dt << endl;
		}

		if (i > 0)
		{
			// Write color plot data file
			of << dtos(frequency,1) << "\t\t"
			   << dtos(i*dE,2) << "\t\t"
				<< peak_height_fr.y << "\t\t"
				<< ((hm_right_fr >= hm_left_fr) ? (hm_right_fr-hm_left_fr)*dt : (period_steps+hm_right_fr-hm_left_fr)*dt) << "\t\t"
				<< global_fr_max[i].y << "\t\t"
				<< peak_base_fr.y << "\t\t"
				<< dtos(tau_adapt_fr[i], 6, true) << "\t\t"
				<< peak_height_O.y << "\t\t"
				<< ((hm_right_O >= hm_left_O) ? (hm_right_O-hm_left_O)*dt : (period_steps+hm_right_O-hm_left_O)*dt) << "\t\t"
				<< global_O_max[i].y << "\t\t"
				<< peak_base_O.y << "\t\t"
				<< dtos(tau_adapt_O[i], 6, true) << endl;
		}
	} // End of for(i)
	cout << endl;
	of << endl; // add empty line to color plot

	// Output of plots for one specific stimulus amplitude
	cout.precision(3); // set number of significant digits in output
	for (int i=0; i<2*stimplot.size(); i++)
	{
		int E;
		char gplcall[100];
		char amp_str[16];

		if (i < stimplot.size())
		{

			// === Plot for firing rate & open probability vs. time ===
#ifdef PLOT_FR_OPROB
			E = int(round(stimplot[i]/dE));
			sprintf(amp_str, "%.2f", stimplot[i]);
			sprintf(gplcall, "gnuplot fr_oprob_%s_PDFplot.gpl", amp_str);

			gpl_fr_oprob[i] << "set term pdf enhanced font \"Sans, 20\" color lw 2.5" << endl;
			gpl_fr_oprob[i] << "set output '" << dateStr("_fr_oprob_") << amp_str << ".pdf'" << endl;
			gpl_fr_oprob[i] << "set xlabel \"t / ms\"" << endl;
			gpl_fr_oprob[i] << "set ylabel \"{/Symbol n}_{inst} / Hz\"" << endl;
			gpl_fr_oprob[i] << "set y2label \"Open probability\" offset -1.5" << endl;
			gpl_fr_oprob[i] << "set key samplen 2" << endl << "set key top right" << endl;
			gpl_fr_oprob[i] << "set ytics nomirror" << endl;
			gpl_fr_oprob[i] << "set y2tics 0.1" << endl << endl;

			gpl_fr_oprob[i] << "# Set margins" << endl;
			gpl_fr_oprob[i] << "set tmargin at screen 0.95" << endl;
			gpl_fr_oprob[i] << "set bmargin at screen 0.23" << endl;
			gpl_fr_oprob[i] << "set rmargin at screen 0.84" << endl;
			gpl_fr_oprob[i] << "set lmargin at screen 0.16" << endl << endl;


			// Maximum points
			if (E == m) // only for highest considered stimulus amplitude
			{
				gpl_fr_oprob[i] << endl << "# Maxima" << endl;

				for (int j=0; j<maxima_O.size(); j++)
				{
					char color[8];
					if ( (maxima_O.at(j).x < sat_begin[E]) // for maxima/half-lives not in saturation choose yellow
						|| (maxima_O.at(j).x < fit_O_max[E].x) ) // for all maxima "left" of global maximum as well
						strcpy(color, "yellow");
					else if ((j % 2) == 0) // for even j's choose blue...
						strcpy(color, "blue");
					else 					// for odd ones dark green
						strcpy(color, "#2F8E20");
					// full circles for maxima
	#define MARK_MAXIMA
	#ifndef MARK_MAXIMA
					gpl_fr_oprob[i] << "# "; // comment out the creation of maximum points
	#endif
					gpl_fr_oprob[i] << "set object circle at second " << maxima_O.at(j).x*dt << "," << maxima_O.at(j).y
						  << " radius char 0.15 fillcolor rgb '" << color << "' fillstyle solid noborder" << endl;
				}
			}

			// Fit maximum point in the open prob. plots - also without MARK_POINTS set
			gpl_fr_oprob[i] << endl << "# Fit maximum" << endl;
			gpl_fr_oprob[i] << "set object circle at second " << fit_O_max[E].x*dt << "," << fit_O_max[E].y
					  << " radius char 0.15 fillcolor rgb 'purple' fillstyle solid noborder" << endl; // mind that this has not the same height as the first maximum circle
																																 // because the fit changes the y-value of fit_O_max!

			// Saturation marking by dashed line
			if (E > 0)
			{
				gpl_fr_oprob[i] << endl << "# Saturation value" << endl;
				gpl_fr_oprob[i] << "set arrow from " << sat_begin[E]*dt << ",graph 0 to " << sat_begin[E]*dt << ",graph 1 nohead lc rgb 'blue' dt 3" << endl;
			}

			gpl_fr_oprob[i] << endl << "set yrange [0:" << global_fr_max[E].y*1.1 << "]" << endl
			                << "set y2range [0:1]" << endl; // alternative: global_O_max[E].y*1.1
			gpl_fr_oprob[i] << "plot [x=" << plot_start << ":" << plot_end
								 << "]" << " '" << dateStr("_fr_oprob_") << amp_str << ".txt'"
							 	 << " using 1:2 title '{/Symbol n}_{inst}(t)' dt 1 lc rgb 'red' with lines, "
							 	 << "'' using 1:3 title 'O(t)' axes x1y2 dt 1 lc rgb '#00e673' with lines, "
							 	 << "'' using 1:4 title 'C(t)' axes x1y2 dt 4 lc rgb '#264d73' with lines";
			gpl_fr_oprob[i] << ",  " << fit_fr_max[E].y - fr_steady[E] << " * exp(-(x-" << double(fit_fr_max[E].x)*dt << ")/" << tau_adapt_fr[E] << ") + "
								 << fr_steady[E] << " notitle"
								 << ",  " << fit_O_max[E].y - O_steady[E] << " * exp(-(x-" << double(fit_O_max[E].x)*dt << ")/" << tau_adapt_O[E] << ") + "
								 << O_steady[E] << " notitle axes x1y2" << endl;
			gpl_fr_oprob[i].close();
			system(gplcall);
			gplscript << gplcall << endl;
#endif // (def PLOT_FR_OPROB)


		// === Plot for light stimulus ===

			// current plot
			sprintf(gplcall, "gnuplot lightstimulus_I_%s_PDFplot.gpl", amp_str);
		}
		else // voltage plot
		{
			E = int(round(stimplot[i-stimplot.size()]/dE));
			sprintf(amp_str, "%.2f", stimplot[i-stimplot.size()]);
			sprintf(gplcall, "gnuplot lightstimulus_V_%s_PDFplot.gpl", amp_str);
			cout << "Created plot for a light stimulus amplitude of " << stimplot[i-stimplot.size()] << " mW/mm^2." << endl;
		}

		// Set gnuplot terminal
		gpl_stim[i] << "set term pdf enhanced font \"Sans, 20\" color lw 2.5" << endl;

		// Current plot
		if (i < stimplot.size())
		{
			gpl_stim[i] << "set output '" << dateStr("_lightstimulus_I_") << amp_str << ".pdf'" << endl;
			gpl_stim[i] << "set ylabel \"I_{total} / nA\"" << endl;
			gpl_stim[i] << "set yrange [" << (global_min[E].y < 0 ? 1.1 : 0.9)*global_min[E].y << ":" << 1.1*global_peak[E].y << "]" << endl;
			if (E > EPSILON)
			{
				gpl_stim[i] << "set y2range [0.0:" << 1.1*stimplot[i] << "]" << endl;
				gpl_stim[i] << "set y2tics " << stimplot[i] << endl;
			}
		}
		// Voltage plot
		else
		{
			gpl_stim[i] << "set output '" << dateStr("_lightstimulus_V_") << amp_str << ".pdf'" << endl;
			gpl_stim[i] << "set ylabel \"V / mV\"" << endl;
			gpl_stim[i] << "set yrange [" << -73 << ":" << 38 << "]" << endl;
			if (E > EPSILON)
			{
				gpl_stim[i] << "set y2range [0.0:" << 1.1*stimplot[i-stimplot.size()] << "]" << endl;
				gpl_stim[i] << "set y2tics " << stimplot[i-stimplot.size()] << endl;
			}
		}

		gpl_stim[i] << "set xlabel \"t / ms\"" << endl;
		if (E > EPSILON)
			gpl_stim[i] << "set y2label \"E / mW/mm^2\"" << endl;
		gpl_stim[i] << "set ytics nomirror" << endl;
		gpl_stim[i] << "#set object circle at first " << global_peak[E].x*dt << "," << global_peak[E].y
				      << " radius char 0.15 fillcolor rgb 'purple' fillstyle solid noborder" << endl << endl;

		gpl_stim[i] << "set tmargin at screen 0.95" << endl;
		gpl_stim[i] << "set bmargin at screen 0.23" << endl << endl;

		// Current plot
		if (i < stimplot.size())
		{
			gpl_stim[i] << "plot [x=" << plot_start << ":" << plot_end << "] '";
			if (E > EPSILON)
				gpl_stim[i] << dateStr("_stim_") << amp_str << ".txt' using 1:4 notitle axes x1y2 dt 1 lc rgb \"gray\" with filledcu x1, '"; // light stimulus
#ifdef STIMULUS_COMPARISON
			// additional curve for I_ChR2 in arbitrary units
			//gpl_stim[i] << dateStr("_stim_") << amp_str << ".txt' using 1:" << stimplot[i]/global_peak_IChR2[E].y
			//														  << "*$5 notitle axes x1y2 dt 3 lc rgb \"blue\" with lines, '"; // ChR2 current
#endif
			gpl_stim[i] << dateStr("_stim_") << amp_str << ".txt' using 1:2 notitle dt 1 lc rgb \"red\" with lines" << endl; // total current
		}
		// Voltage plot
		else
		{
			gpl_stim[i] << "plot [x=" << plot_start << ":" << plot_end<< "] '";
			if (E > EPSILON)
				gpl_stim[i] << dateStr("_stim_") << amp_str << ".txt' using 1:4 notitle axes x1y2 dt 1 lc rgb \"gray\" with filledcu x1, '"; // light stimulus

			gpl_stim[i] << dateStr("_stim_") << amp_str << ".txt' using 1:3 notitle dt 1 lc rgb \"red\" with lines, " // voltage
							<< neur.getSteadyVoltage(stimplot[i-stimplot.size()], 0) << " notitle lc rgb \"green\"" << endl; // steady-state voltage
		}
		gpl_stim[i].close();
		system(gplcall);
		gplscript << gplcall << endl;

		// === Plot for external current stimulus ===
#ifdef STIMULUS_COMPARISON
		char gplcall2[100];
		// Current plot
		if (i < stimplot.size())
			sprintf(gplcall2, "gnuplot currentstimulus_I_%s_PDFplot.gpl", amp_str);
		// Voltage plot
		else
			sprintf(gplcall2, "gnuplot currentstimulus_V_%s_PDFplot.gpl", amp_str);
		gpl_stim2[i] << "set term pdf enhanced font \"Sans, 20\" color lw 2.5" << endl;
		// Current plot
		if (i < stimplot.size())
		{
			gpl_stim2[i] << "set output '" << dateStr("_currentstimulus_I_") << amp_str << ".pdf'" << endl;
			gpl_stim2[i] << "set ylabel \"I_{total} / nA\"" << endl;
			gpl_stim2[i] << "set yrange [" << (global_min2[E].y < 0 ? 1.1 : 0.9)*global_min2[E].y << ":" << 1.1*global_peak2[E].y << "]" << endl;

			// Check if global I_ChR2 peak is within plot range
			if ( global_peak_IChR2[E].x*dt >= plot_start && global_peak_IChR2[E].x*dt <= plot_end )
			{
				gpl_stim2[i] << "set y2range [0.0:" << 1.1* ( global_peak_IChR2[E].y > neur2.getSteadyCurrent(stimplot[i], 0) ? global_peak_IChR2[E].y
																	: neur2.getSteadyCurrent(stimplot[i], 0) ) << "]" << endl;
			}
			else
			{
				gpl_stim2[i] << "set y2range [0.0:" << 1.5*neur2.getSteadyCurrent(stimplot[i], 0) << "]" << endl;
			}
			gpl_stim2[i] << "set y2tics " << round(neur2.getSteadyCurrent(stimplot[i], 0)*1000.0)/1000.0 << endl;
		}
		// Voltage plot
		else
		{
			gpl_stim2[i] << "set output '" << dateStr("_currentstimulus_V_") << amp_str << ".pdf'" << endl;
			gpl_stim2[i] << "set ylabel \"V / mV\"" << endl;
			gpl_stim2[i] << "set yrange [" << -73 << ":" << 38 << "]" << endl;
			gpl_stim2[i] << "set y2range [0.0:" << 1.1*neur2.getSteadyCurrent(stimplot[i-stimplot.size()], 0) << "]" << endl;
			gpl_stim2[i] << "set y2tics " << round(neur2.getSteadyCurrent(stimplot[i-stimplot.size()], 0)*1000.0)/1000.0 << endl;
		}
		gpl_stim2[i] << "set xlabel \"t / ms\"" << endl;
		gpl_stim2[i] << "set y2label \"I_{cst}, I_{ChR2} / nA\"" << endl;
		gpl_stim2[i] << "set ytics nomirror" << endl << endl;

		gpl_stim2[i] << "set tmargin at screen 0.95" << endl;
		gpl_stim2[i] << "set bmargin at screen 0.23" << endl << endl;

		// Current plot
		if (i < stimplot.size())
		{
			gpl_stim2[i] << "set key samplen 2" << endl << "set key top right" << endl; // key
			gpl_stim2[i] << "plot [x=" << plot_start << ":" << plot_end;
			gpl_stim2[i] << "] '" << dateStr("_stim_") << amp_str << ".txt' using 1:8 title 'Stim.' axes x1y2 dt 1 lc rgb \"gray\" with filledcu x1, '" // current stim.
										 << dateStr("_stim_") << amp_str << ".txt' using 1:5 title 'I_{ChR2}' axes x1y2 dt 1 lc rgb \"blue\" with lines, '" // ChR2 current
										 << dateStr("_stim_") << amp_str << ".txt' using 1:6 title 'I_{total}' dt 1 lc rgb \"red\" with lines" << endl;  // total current
		}
		// Voltage plot
		else
		{
			gpl_stim2[i] << "plot [x=" << plot_start << ":" << plot_end;
			gpl_stim2[i] << "] '" << dateStr("_stim_") << amp_str << ".txt' using 1:8 notitle axes x1y2 dt 1 lc rgb \"gray\" with filledcu x1, '" // current stimulus
										 << dateStr("_stim_") << amp_str << ".txt' using 1:7 notitle dt 1 lc rgb \"red\" with lines" << endl; // voltage
		}
		gpl_stim2[i].close();
		system(gplcall2);
		gplscript << gplcall2 << endl;
#endif
	}

	// Write gnuplot file for firing rate over irradiance/current stimulus
	gpl_if << "set term pdf enhanced font \"Sans, 20\" color lw 2.5" << endl;
	gpl_if << "set output '" << dateStr("_firingrate.pdf'") << endl;
	gpl_if << "set log x" << endl;
	gpl_if << "set xlabel \"Ê / mW/mm^2\"" << endl;
#ifdef STIMULUS_COMPARISON
	gpl_if << "set x2label \"Î_{cst} / nA\"" << endl;
	gpl_if << "set key top left" << endl;
	gpl_if << "set key samplen 2" << endl;
	gpl_if << "set log x2" << endl;
	gpl_if << "set x2range [" << neur2.getSteadyCurrent(0.9*dE*dE_average, 0) << ":" << neur2.getSteadyCurrent(E_max, 0) << "]" << endl;
	gpl_if << "set xtics nomirror" << endl;
	gpl_if << "set x2tics" << endl;
#endif
	gpl_if << "set ylabel \"{/Symbol n} / Hz\"" << endl << endl;
	gpl_if << "set tmargin at screen 0.95" << endl;
	gpl_if << "set bmargin at screen 0.23" << endl << endl;
	gpl_if << "plot [x=" << 0.9*dE*dE_average << ":" << E_max << "] '" << dateStr("_data.txt'")
#ifdef STIMULUS_COMPARISON
			 << " using 1:2 title \"Ê\" with lines, '" << dateStr("_data.txt'") << " using 3:4 title \"Î_{cst}\" axes x2y1 with lines" << endl;
#else
			 << " using 1:2 notitle with lines" << endl;
#endif
	gpl_if.close();

#ifdef STIMULUS_COMPARISON
#ifdef PHOTOCURRENT_SATURATION
	// Write gnuplot file for photocurrent over light stimulus
	gpl_ichr2 << "set term pdf enhanced font \"Sans, 20\" color lw 2.5" << endl;
	gpl_ichr2 << "set output '" << dateStr("_light_current.pdf'") << endl;
	gpl_ichr2 << "set xlabel \"E / mW/mm^2\"" << endl;
	gpl_ichr2 << "set ylabel \"I_{ChR2,sat} / nA\"" << endl << endl;
	gpl_ichr2 << "set tmargin at screen 0.95" << endl;
	gpl_ichr2 << "set bmargin at screen 0.23" << endl << endl;
	gpl_ichr2 << "plot [x=" << dE*dE_average << ":" << E_max << "] '" << dateStr("_data.txt'") << " using 1:3 notitle with lines, "
				 << neur.getSaturationCurrent(0) << " notitle dt 2" << endl;
	gpl_ichr2.close();
#endif
#endif

	// Close files that have remained open
	data.close();
	of.close();

	// Call GNUplot to process the files
	system("gnuplot firingrate.gpl");
	gplscript << "gnuplot firingrate.gpl" << endl;
#ifdef STIMULUS_COMPARISON
	system("gnuplot light_current.gpl");
	gplscript << "gnuplot light_current.gpl" << endl;
#endif
	gplscript.close(); // script is not executed, it is only created in case it is needed later on
	errlog.close();

	// Display final information and save parameters
	int tsec = timeMeasure(false);
	char el_time [32];
	sprintf(el_time, "Elapsed time: %d min %d sec", int(floor(tsec / 60)), int(tsec % 60));
	cout << el_time << endl;
	saveParams(el_time); // additional informaton: elapsed time, threshold frequency for constant illumination
	if (chdir("../..") == -1)
		showChDirErrMessage();
	cout << separator << endl;

	// Delete allocated memory
	delete[] spikes_lst;
	delete[] txt_stim;
	delete[] gpl_stim;
	delete[] txt_fr_oprob;
	delete[] gpl_fr_oprob;
#ifdef STIMULUS_COMPARISON
	delete[] spikes_cst;
	delete[] gpl_stim2;
#endif
	return 0;
}

/*** setParams ***
 * Sets the simulation parameters on given values and resets neuron(s) */
void setParams(double _t_max, double _E_max, double _t_pulse, double _frequency, double _const_current)
{
	t_max = _t_max;
	E_max = _E_max;
	t_pulse = _t_pulse;
	frequency = _frequency;
	neur.setConstCurrent(_const_current);
	neur.setCouplingStrengths(0., 0., 0., 0.);
	neur.reset();
	neur.setConnections(0.);
#ifdef STIMULUS_COMPARISON
	neur2.setConstCurrent(_const_current);
	neur2.setCouplingStrengths(0., 0., 0., 0.);
	neur2.reset();
	neur2.setConnections(0.);
#endif
}

/*** Constructor ***
 * Sets all parameters on given values and calls constructors for Neuron instances */
NeuronSimulation(const double _dt, const double _dE, const int _dE_average, double _t_offset, double _t_max, double _E_max, const vector<double> _stimplot,
					  string _cplotfile, double _plot_start, double _plot_end)
	: offset(int(ceil(_t_offset / dt))), neur(_dt, offset, int(sqrt(num_trials)), 0), dt(_dt), dE(_dE), dE_average(_dE_average), stimplot(_stimplot), plot_start(_plot_start), plot_end(_plot_end)
#ifdef STIMULUS_COMPARISON
		, neur2(_dt, offset, int(sqrt(num_trials)), 0)
#endif
{
	cplotfile = _cplotfile;
	t_max = _t_max;
	E_max = _E_max;
}

};
