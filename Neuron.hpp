/**************************************************************************
 * Leaky Integrate-and-Fire model of a single neuron with ChR2 channels
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

#include <random>

using namespace std;

#include "ChR2_3state.hpp"

#define TYPE_EXC 1
#define TYPE_INH 2
//#define INWARD_RECT // specifies if voltage-dependent inward rectification shall be considered

/*** Neuron class ***
 * Represents one neuron with a stochastically modelled population of ChR2 */
class Neuron {

private:

/*** Computational parameters ***/
double dt; // ms, one time step for numerical simulation
int tb_offset; // the number of offset timesteps before stimulation

/*** State variables ***/
double V; // mV, the current membrane potential
double I_ChR2; // nA, the current evoked by the ChR2 channels
double I_cst; // nA, the externally applied stimulus current
double I_ext; // nA, the current evoked by external synaptic inputs (computed using an OU process with mean 0)
double I_int; // nA, the synaptic input from other network neurons affecting this neuron
double refractory; // ms, time span until refractory period is over
vector<int> spike_history; // vector of all spike times (in units of time steps) in the process since the last reset
int inh_incoming; // number of incoming inhibitory connections in a network
int exc_incoming; // number of incoming excitatory connections in a network
int outgoing; // number of outgoing connections in a network
Stimulus cst; // current stimulus for this neuron
Stimulus lst; // light stimulus for this neuron
bool cst_set; // true if a current stimulus has been set since the last reset
bool lst_set; // true if a light stimulus has been set since the last reset
ChR2 ch; // the stochastic ChR2 channel instance
minstd_rand0 rg; // default uniform generator for random numbers (seed is chosen in constructor)
normal_distribution<double> norm_dist; // normal distribution to obtain Gaussian white noise, constructed in Neuron class constructor

protected:

/*** Physical parameters ***/
double tau_mem; // ms, the membrane time constant
double tau_ref; // ms, refractory period - has to be at least one time step!
double G_mem; // µS, conductance of the cell membrane
double V_rev; // mV, the reversal potential of the neuron
double V_reset; // mV, the reset potential of the neuron
double V_th; // mV, the threshold potential of the neuron
double V_spike; // mV, the height of an action potential
int N_ChR2; // expression level/number of ChR2 channels in this neuron
#ifdef INWARD_RECT
double V_inw1; // mV, constant for voltage-dependent inward rectification
double V_inw2; // mV, constant for voltage-dependent inward rectification
#endif
double tau_OU; // ms, correlation time of the Ornstein-Uhlenbeck process
double sigma_WN; // nA s^1/2, standard deviation of the Gaussian white noise driving the OU process
double I_const; // nA, the mean of the external synaptic inputs (i.e., of the OU process)
int type; // the type of this neuron (inhibitory/excitatory - TYPE_INH/TYPE_EXC)


/*** normalRandomNumber ***
 * Returns a random number drawn from a normal distribution with standard deviation 1 and mean 0 *
 * - return: the random number of type double (technically in units of sqrt(s)) */
double normalRandomNumber()
{
	double nd = norm_dist(rg);

	return nd;
}

public:

/*** saveNeuronParams ***
 * Saves all the neuron parameters (including the channel parameters) to a given file */
void saveNeuronParams(ofstream *f) const
{
	*f << endl;
	*f << "Neuron parameters:" << endl;
	*f << "tau_mem = " << tau_mem << " ms" << endl;
	*f << "tau_ref = " << tau_ref << " ms" << endl;
	*f << "G_mem = " << G_mem << " µS" << endl;
	*f << "V_rev = " << V_rev << " mV" << endl;
	*f << "V_reset = " << V_reset << " mV" << endl;
	*f << "V_th = " << V_th << " mV" << endl;
	*f << "V_spike = " << V_spike << " mV" << endl;
	*f << "N_ChR2 = " << N_ChR2 << endl;
#ifdef INWARD_RECT
	*f << "V_inw1 = " << V_inw1 << " mV" << endl;
	*f << "V_inw2 = " << V_inw2 << " mV" << endl;
#endif
	*f << "tau_OU = " << tau_OU << " ms" << endl;
	*f << "sigma_WN = " << sigma_WN << " nA s^1/2" << endl;
	*f << "I_const = " << I_const << " nA" << endl;

	ch.saveChannelParams(f);
}

/*** getNumberIncoming ***
 * Returns the number of either inhibitory or excitatory incoming connections to this neuron *
 * from other neurons in a network *
 * - int type: the type of incoming connections (inh./exc.)
 * - return: the number of incoming connections */
int getNumberIncoming(int type) const
{
	if (type == TYPE_INH)
		return inh_incoming;
	else if (type == TYPE_EXC)
		return exc_incoming;
}

/*** getNumberOutgoing ***
 * Returns the number of connections outgoing from this neuron to other *
 * neurons in a network *
 * - return: the number of outgoing connections */
int getNumberOutgoing() const
{
	return outgoing;
}

/*** incNumberIncomingInh ***
 * Increases the number of incoming inhibitory connections to this neuron (only to be used while *
 * a network is being built) */
void incNumberIncomingInh()
{
	inh_incoming++;
}

/*** incNumberIncomingExc ***
 * Increases the number of incoming excitatory connections to this neuron (only to be used while *
 * a network is being built) */
void incNumberIncomingExc()
{
	exc_incoming++;
}

/*** incNumberOutgoing ***
 * Increases the number of outgoing connections from this neuron (only to be used while *
 * a network is being built) */
void incNumberOutgoing()
{
	outgoing++;
}

/*** getVoltage ***
 * Returns the membrane potential of the neuron *
 * - return: the membrane potential in mV */
double getVoltage() const
{
	return V;
}

/*** getCurrent ***
 * Returns total current affecting the neuron *
 * - return: the instantaneous current in nA */
double getCurrent() const
{
	return I_ChR2+I_cst+I_const+I_ext;
}

/*** getChR2Current ***
 * Returns current evoked by the ChR2 molecules *
 * - return: the instantaneous ChR2 channel current in nA */
double getChR2Current() const
{
	return I_ChR2;
}

/*** getStimulusCurrent ***
 * Returns current evoked by external stimulation *
 * - return: the instantaneous current stimulus in nA */
double getStimulusCurrent() const
{
	return I_cst;
}

/*** getFluctCurrent ***
 * Returns fluctuating current evoked by external synaptic inputs (OU process around zero)*
 * - return: the instantaneous fluctuating external synaptic current in nA */
double getFluctCurrent() const
{
	return I_ext;
}

/*** getConstCurrent ***
 * Returns the mean current elicited by external inputs *
 * - return: the constant current in nA */
double getConstCurrent() const
{
	return I_const;
}

/*** getSynapticCurrent ***
 * Returns the internal synaptic current that arrived in the previous time step *
 * - return: the synaptic current in nA */
double getSynapticCurrent() const
{
	return I_int;
}

/*** getActivity ***
 * Returns true if the neuron is spiking in this instant of duration dt *
 * - return: whether neuron is firing or not */
bool getActivity() const
{
	if (V == V_spike)
		return true;
	else
		return false;
}

/*** getSpikeTime ***
 * Returns the spike time for a given spike number *
 * - int n: the number of the spike (in temporal order, starting with 1)
 * - return: the spike time in units of time bins for the n-th spike (or -1 if it does not exist) */
int getSpikeTime(int n) const
{
	if (n <= spike_history.size() && n >= 1)
		return spike_history.at(n-1);
	else
		return -1;
}

/*** getSpikeCount ***
 * Returns the number of spikes that have occurred since the last reset *
 * - return: the number of spikes */
int getSpikeCount() const
{
	return spike_history.size();
}

/*** getSteadyVoltage ***
 * Returns the steady-state voltage (V_ss) for a constant light stimulus of given irradiance *
 * - double irradiance: the irradiance value in mW/mm^2 *
 * - return: the steady-state voltage at this irradiance in mV */
double getSteadyVoltage(double irradiance) const
{
#ifdef THREE_STATE_CHR2
	double ph_inf, gamma_r;
	ph_inf = ch.getPhotonInflux(irradiance);
	gamma_r = ch.getGammaR();

	return (1.0/(14.0*ch.gamma_d0*(ph_inf + gamma_r)*G_mem))*(1250.0*ph_inf*gamma_r*(G_mem + ch.G_ChR2*1e-9*N_ChR2) + ch.gamma_d0*(ph_inf + gamma_r)*
			 (7.0*I_const + G_mem*(760.0 + 7.0*V_rev)) - sqrt(50000.0*ph_inf*gamma_r*(76.0*ph_inf*ch.gamma_d0 + 125.0*ph_inf*gamma_r + 76.0*ch.gamma_d0*gamma_r)*
			 ch.G_ChR2*1e-9*G_mem*N_ChR2 + pow(ch.gamma_d0*gamma_r*(7.0*I_const + G_mem*(-760+7*V_rev)) + ph_inf*(-1250*gamma_r*(G_mem -
			 ch.G_ChR2*1e-9*N_ChR2) + ch.gamma_d0*(7.0*I_const + G_mem*(-760+7.0*V_rev))), 2))); // mind the 1e-9!
#endif
}

/*** getSteadyCurrent ***
 * Returns the steady-state current (I_ss) for a constant light stimulus of given irradiance *
 * - double irradiance: the irradiance value in mW/mm^2 *
 * - return: the steady-state current at this irradiance in nA */
double getSteadyCurrent(double irradiance) const
{
#ifdef THREE_STATE_CHR2
	double ph_inf, gamma_r;
	ph_inf = ch.getPhotonInflux(irradiance);
	gamma_r = ch.getGammaR();

	return -((1.0/(14.0*ch.gamma_d0*(ph_inf + gamma_r)))*(-1250.0*ph_inf*gamma_r*(G_mem + ch.G_ChR2*1e-9*N_ChR2) + ch.gamma_d0*(ph_inf + gamma_r)*
			   (7.0*I_const + G_mem*(-760.0 + 7.0*V_rev)) + sqrt(50000.0*ph_inf*gamma_r*(76.0*ph_inf*ch.gamma_d0 + 125.0*ph_inf*gamma_r + 76.0*ch.gamma_d0*gamma_r)*
			   ch.G_ChR2*1e-9*G_mem*N_ChR2 + pow(ch.gamma_d0*gamma_r*(7.0*I_const + G_mem*(-760+7*V_rev)) + ph_inf*(-1250*gamma_r*(G_mem -
			   ch.G_ChR2*1e-9*N_ChR2) + ch.gamma_d0*(7.0*I_const + G_mem*(-760+7.0*V_rev))), 2)))); // mind the 1e-9!
#else
	return 0;
#endif
}

/*** getSaturationCurrent ***
 * Returns the saturation current for a constant light stimulus in the infinite irradiance limit *
 * - return: the saturation current for the set parameters in nA */
double getSaturationCurrent() const
{
#ifdef THREE_STATE_CHR2
	double gamma_r = ch.getGammaR();
	return -((1.0/(14.0*ch.gamma_d0))*(-760.0*ch.gamma_d0*G_mem - 1250.0*gamma_r*G_mem + 7.0*ch.gamma_d0*I_const - 1250.0*gamma_r*ch.G_ChR2*1e-9*N_ChR2 +
   		   7.0*ch.gamma_d0*G_mem*V_rev + sqrt(3800000.0*ch.gamma_d0*gamma_r*ch.G_ChR2*1e-9*G_mem*N_ChR2 +
     		   6250000.0*pow(gamma_r,2)*ch.G_ChR2*1e-9*G_mem*N_ChR2 + pow(-1250.0*gamma_r*(G_mem - ch.G_ChR2*1e-9*N_ChR2) +
       	   ch.gamma_d0*(7.0*I_const + G_mem*(-760.0 + 7.0*V_rev)),2)))); // mind the 1e-9!
#else
	return 0;
#endif
}

/*** getOpenProb ***
 * Returns the probability that a ChR2 channel is open *
 * - return: the open probability */
double getOpenProb() const
{
	return ch.getO();
}

/*** getClosedProb ***
 * Returns the probability that a ChR2 channel is closed *
 * - return: the closed probability */
double getClosedProb() const
{
	return ch.getC();
}

/*** processTimeStep ***
 * Processes one time step (of duration dt) for the neuron *
 * - int tb_step: time step at which to evaluate stimulus (< 0 before stimulus onset) *
 * - double _I_int [optional]: synaptic current evoked by other neurons, only of importance in a network */
void processTimeStep(int tb_step, double _I_int = 0.0)
{
	double sigma_OU = sigma_WN / sqrt(2.*tau_OU/1000.); // leads to the sigma of the OU process; times 1000 because required unit is [s]
#ifdef INWARD_RECT
	double G = (1-exp(-V/V_inw1))/(V/V_inw2); // from Grossman et al., 2011
#else
	double G = 1.;
#endif

	if (lst_set)
		ch.setIrradiance(lst.getStimulusAt(tb_step-tb_offset)); // set the current light stimulus value as irradiance for the channel
	ch.processTimeStep(tb_step, V); // process time step for channels
	if (cst_set)
		I_cst = cst.getStimulusAt(tb_step-tb_offset); // get stimulus current in nA
	if (V < 0.) // no current during action potential
		I_ChR2 = - N_ChR2 * 1e-9 * ch.G_ChR2 * G * V * ch.getCalcO(); // compute channel current in nA (therefore, multiply G_ChR2 by 10^-9)
	else
		I_ChR2 = - N_ChR2 * 1e-9 * ch.G_ChR2 * G * V_reset * ch.getCalcO(); // use V_reset as voltage for time bin in which a spike occurred

#ifndef LIF_FR_TEST
#ifdef DELTA_SYNAPSES
	I_ext = normalRandomNumber() * sqrt(1000./dt) * sigma_WN; // sqrt(1/dt) has been added in revised version in 2017
#else
	I_ext = I_ext * exp(-dt/tau_OU) + normalRandomNumber() * sqrt(1. - exp(-2.*dt/tau_OU)) * sigma_OU; // compute external synaptic input in nA
#endif
#endif
	I_int = _I_int;

	V = V * exp(-dt/tau_mem) + (V_rev + (I_const + I_ext + I_cst + I_ChR2 + I_int) / G_mem) * (1. - exp(-dt/tau_mem)); // compute mem. pot. in mV (analytical solution)

	if (refractory > 0.0) // if during refractory period
	{
		V = V_reset;
		refractory -= dt;
	}
	else if (V >= V_th) // threshold crossing
	{
		V = V_spike;
		refractory = tau_ref;
		spike_history.push_back(tb_step);
	}
}

/*** setCurrentStimulus ***
 * Sets a current stimulus for the neuron *
 * - Stimulus& _cst: shape of one stimulus period */
void setCurrentStimulus(const Stimulus& _cst)
{
	cst = Stimulus(_cst); // use copy constructor
	cst_set = true;
}

/*** setLightStimulus ***
 * Sets a light stimulus for the ChR2 molecules in the neuron *
 * - Stimulus& _lst: shape of one stimulus period */
void setLightStimulus(const Stimulus& _lst)
{
	lst = _lst; // use assignment operator
	lst_set = true;
}

/*** getCurrentStimulusAt ***
 * Returns the current stimulus magnitude at a given time *
 * - tb_step: time step at which to evaluate stimulus
 * - return: stimulus at given time (if stimulus is not set, 0.0) */
double getCurrentStimulusAt(int tb_step) const
{
	if (cst_set)
		return cst.getStimulusAt(tb_step-tb_offset);
	else
		return 0.0;
}

/*** getLightStimulusAt ***
 * Returns the light stimulus magnitude at a given time *
 * - tb_step: time step at which to evaluate stimulus
 * - return: stimulus at given time (if stimulus is not set, 0.0) */
double getLightStimulusAt(int tb_step) const
{
	if (lst_set)
		return lst.getStimulusAt(tb_step-tb_offset);
	else
		return 0.0;
}

/*** multiplyLightStimulus ***
 * Multiplies the set light stimulus by a real number *
 * - double r: number to multiply */
void multiplyLightStimulus(double r)
{
	lst.multiplyBy(r);
}

/*** multiplyCurrentStimulus ***
 * Multiplies the set current stimulus by a real number *
 * - double r: number to multiply */
void multiplyCurrentStimulus(double r)
{
	cst.multiplyBy(r);
}

/*** setConstCurrent ***
 * Sets the constant current (mean of the OU process) to a newly defined value *
 * - double _I_const: constant current in nA */
void setConstCurrent(double _I_const)
{
	I_const = _I_const;
}

/*** setTauOU ***
 * Sets the time constant of the Ornstein-Uhlenbeck process *
 * (synaptic time constant of assumed input synapses) *
 * - double _tau_OU: the synaptic time constant */
void setTauOU(double _tau_OU)
{
	tau_OU = _tau_OU;
}

/*** setType ***
 * Sets the type of this neuron (inhbitory/excitatory) *
 * - int _type: the neuron type */
void setType(int _type)
{
	type = _type;
}

/*** getType ***
 * Returns the type of this neuron (inhbitory/excitatory) *
 * - return: the neuron type */
int getType()
{
	return type;
}

/*** resetConnections ***
 * Resets the number of connections */
void resetConnections()
{
	inh_incoming = 0;
	exc_incoming = 0;
	outgoing = 0;
}

/*** reset ***
 * Resets neuron to initial state */
void reset()
{
	vector<int>().swap(spike_history); // additionally to clearing the vector, reallocates it
	V = V_rev;
	I_ChR2 = 0.0; // assume closed channel
	I_cst = 0.0;
	I_ext = 0.0;
	I_int = 0.0;
	refractory = 0.0; // neuron ready to fire
	cst_set = false;
	lst_set = false;

	rg.seed(getClockSeed()); // set new seed by clock's epoch
	norm_dist.reset(); // reset the normal distribution for random numbers
	ch.reset();
}

/*** Constructor ***
 * Sets all parameters on experimentally determined values *
 * - double _dt: the duration *
 * - int _tb_offset: the number of offset timesteps before stimulation */
Neuron(const double _dt, const int _tb_offset) :
	dt(_dt), tb_offset(_tb_offset), ch(_dt), norm_dist(0.0,1.0), rg(getClockSeed())
{
	tau_mem = 10.0; // from Dayan & Abbott, fig. 5.
	tau_ref = 3.0; // estimated to model afterhypolarization in combination with V_reset
	G_mem = 0.1; // from Dayan & Abbott, fig. 5.
	V_rev = -65.0; // from Dayan & Abbott, fig. 5.
	V_reset = -70.0; // estimated to model afterhypolarization in combination with tau_ref (see also Dayan & Abbott, p. 4)
	V_th = -55.0; // from Dayan & Abbott, p. 162
	V_spike = 35.0; // estimated (only relevant for visualization)
	N_ChR2 = 60000; // low estimation
	//N_ChR2 = 300000; // high estimation
#ifdef INWARD_RECT
	V_inw1 = 40.; // from Grossman et al., 2011
	V_inw2 = 15.; // from Grossman et al., 2011
#endif
	sigma_WN = 0.010; // estimated
	I_const = 0.116; // estimated (just before steady-state is above threshold), is overwritten by setConstCurrent()
	setTauOU(5.0); // estimated

	reset();
	resetConnections();
}

/*** Destructor ***
 * Frees the allocated memory */
~Neuron()
{

}


};
