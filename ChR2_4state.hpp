/**************************************************************************
 * Four-state model of a single channelrhodopsin-2 channel
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

#define FOUR_STATE_CHR2
#define FIT_WITH_MAX // for initial conditions C1=1, C2=0 (otherwise, FIT_WITHOUT_MAX)
#define h 6.62606957e-5 // 10^{-29} J*s, Planck constant, 10^{-5} instead of 10^{-34} for canceling out cross section and lambda magnitude in computation
#define c 299792458.0 // m/s, speed of light

#include "Stimulus.hpp"

/*** Channelrhodopsin-2 class ***
 * Stochastically models one ChR2 channel molecule with four states O1, O2, *
 * C1 and C2 (of which two are open and two are closed) */
class ChR2
{
	friend class Neuron; // to allow Neuron class access to parameters

private:
	// state variables:
	double C1; // probability for first closed state
	double O1; // probability for first open state
	double C2; // probability for second closed state
	double O2; // probability for second open state
	double p; // activation rate probability function for the ion channel
	//int t_light; // time step at which light occurred last
	int t_light_begin; // time step at which light began to shine

protected:
	// parameters:
	double tau_ChR; // ms, activation time constant for the ion channel
	double ret_cross_section; // 10^{-20} m^2, retinal cross-section
	double G_ChR2; // fS, single channel conductance
	double gamma; // ratio of the conductances of O1 and O2
	double lambda_max; // nm, the wavelength of maximum absorption
	double wloss; // loss factor for channel environment
	double gamma_d1; // 1/s, closing rate O1 -> C1
	double gamma_d2; // 1/s, closing rate O2 -> C2
	double gamma_th; // 1/s, thermal transition rate C2 -> C1
	double e_12d; // 1/s, basic transition rate O1 -> O2
	double e_21d; // 1/s, basic transition rate O2 -> O1
	double e_12; // 1/s, transition rate O1 -> O2, dependent on irradiance
	double e_21; // 1/s, transition rate O2 -> O1, dependent on irradiance
	double efficiency1; // probability that a photon absorbed by the retinal activates the channel in state O1
	double efficiency2; // probability that a photon absorbed by the retinal activates the channel in state O2

	// others:
	const double dt; // ms, duration of one time step
	double photon_influx1; // 1/s, effective photons per second the channel is exposed to in state C1
	double photon_influx2; // 1/s, effective photons per second the channel is exposed to in state C2
	double S_of_theta; // sigmoidal function of the optical stimulation protocol

public:
	/*** saveChannelParams ***
	 * Saves all the channel parameters to a given file */
	void saveChannelParams(ofstream *f) const
	{
		*f << endl;
		*f << "ChR2 parameters:" << endl;
		*f << "tau_ChR = " << tau_ChR << " ms" << endl;
		*f << "ret_cross_section = " << ret_cross_section << "e-20 m^2" << endl;
		*f << "G_ChR2 = " << G_ChR2 << " fS" << endl;
		*f << "gamma = " << gamma << endl;
		*f << "lambda_max = " << lambda_max << " nm" << endl;
		*f << "wloss = " << wloss << endl;
		*f << "gamma_d1 = " << gamma_d1 << " 1/s" << endl;
		*f << "gamma_d2 = " << gamma_d2 << " 1/s" << endl;
		*f << "gamma_th = " << gamma_th << " 1/s" << endl;
		*f << "e_12d = " << e_12d << " 1/s" << endl;
		*f << "e_21d = " << e_21d << " 1/s" << endl;
		*f << "e_12 = " << e_12 << " 1/s" << endl;
		*f << "e_21 = " << e_21 << " 1/s" << endl;
		*f << "efficiency1 = " << efficiency1 << endl;
		*f << "efficiency2 = " << efficiency2 << endl;
		*f << "photon_influx1(E_max) = " << photon_influx1 << " 1/s" << endl;
		*f << "photon_influx2(E_max) = " << photon_influx2 << " 1/s" << endl;
		*f << "S_of_theta = " << S_of_theta << " 1/s" << endl;
	}

	/*** getO1 ***
	 * - return: the current value of the first open state probability */
	double getO1() const
	{
		return O1;
	}

	/*** getO2 ***
	 * - return: the current value of the second open state probability */
	double getO2() const
	{
		return O2;
	}

	/*** getC1 ***
	 * - return: the current value of the first closed state probability */
	double getC1() const
	{
		return (1-O1-O2-C2);
	}

	/*** getC2 ***
	 * - return: the current value of the second closed state probability */
	double getC2() const
	{
		return C2;
	}

	/*** getO ***
	 * Returns the total probability to be in an open state *
	 * - return: the current value of the total open state probability */
	double getO() const
	{
		return O1 + O2;
	}

	/*** getCalcO ***
	 * Returns the probability to be in an open state, needed for calculation *
	 * - return: the current value of the weighted sum of the first and second open state probability */
	double getCalcO() const
	{
		return O1 + gamma*O2;
	}

	/*** getC ***
	 * Returns the total probability to be in a closed state *
	 * - return: the current value of the total closed state probability */
	double getC() const
	{
		return getC1() + C2;
	}


	/*** setIrradiance ***
	 * Calculates the photon influx for a given irradiance and therefrom calculates the irradiance-dependent
	 * quantities (for not each step has to be calculated)
	 * - double light_intensity: irradiance in mW/mm^2 */
	void setIrradiance(double irradiance)
	{
		photon_influx1 = efficiency1 * (irradiance * 1000.0 * ret_cross_section * lambda_max) / (h * c * wloss);
		photon_influx2 = efficiency2 * (irradiance * 1000.0 * ret_cross_section * lambda_max) / (h * c * wloss);
		e_12 = e_12d + 5.0 * log(1 + irradiance/0.024); // from Williams et al., 2013
		e_21 = e_21d + 4.0 * log(1 + irradiance/0.024); // from Williams et al., 2013
		S_of_theta = 0.5 * (1 + tanh(120 * (100 * irradiance - 0.1))); // equals 1 for all values greater than 0.002 mW/mm^2
  	}

	/*** setPhotonInflux ***
	 * Sets the photon influxes directly (for frequency response check only)
	 * - double _photon_influx0: brutto photon influx in 1/s */
	void setPhotonInflux(double _photon_influx0)
	{
		double irradiance = _photon_influx0 * (h * c * wloss) / (1000.0 * ret_cross_section * lambda_max); // needed to compute following quantities
		e_12 = e_12d + 5.0 * log(1 + irradiance/0.024); // from Williams et al., 2013
		e_21 = e_21d + 4.0 * log(1 + irradiance/0.024); // from Williams et al., 2013
		S_of_theta = 0.5 * (1 + tanh(120 * (100 * irradiance - 0.1))); // equals 1 for all values greater than 0.002 mW/mm^2

		photon_influx1 = efficiency1 * _photon_influx0;
		photon_influx2 = efficiency2 * _photon_influx0;
  	}

   	/*** processTimeStep ***
	 * Processes one time step (of duration dt) for the channel's state probabilities
	 * - int t_step: time step at which to evaluate stimulus (only of importance when LightStimulus is set) *
	 * - double mem_pot [optional]: the current membrane potential (not used in this model) */
	void processTimeStep(int t_step, double mem_pot = -70.0)
	{
		double gamma_a1, gamma_a2; // 1/s, activation rates C1 -> O1 and C2 -> O2, respectively
		double p; // activation probability

		// Activation function computation, cf. Nikolic et al. 2009 Appendix 1
		if (photon_influx1 >= 1e-30) // light occurred
		{
			if (t_light_begin == 0)
				t_light_begin = t_step;
			p = (1 - exp(- (t_step-t_light_begin)*dt / tau_ChR)); // calculation in milliseconds (dt, tau_ChR in ms)
		}
		else // no light
		{
			t_light_begin = 0;
			p = 0.0;
		}

		gamma_a1 = photon_influx1 * p;
		gamma_a2 = photon_influx2 * p;

		double dO1 = (gamma_a1 - (gamma_a1 + gamma_d1 + e_12) * O1 + (e_21 - gamma_a1) * O2 - gamma_a1 * C2) * dt/1000.0;
		double dO2 = (gamma_a2 * C2 - (gamma_d2 + e_21) * O2 + e_12 * O1) * dt/1000.0;
		double dC2 = (gamma_d2 * O2 - (gamma_a2 + gamma_th) * C2) * dt/1000.0;

		C2 = C2 + dC2;
  		O1 = O1 + dO1;
		O2 = O2 + dO2;
  	}

	/*** reset ***
	 * Resets probabilities to initial state ('closed' state has probability 1) */
	void reset()
	{
		C1 = 1;
		C2 = 0;
		O1 = 0;
		O2 = 0;
		p = 0;
		photon_influx1 = 0;
		photon_influx2 = 0;
		t_light_begin = 0;
	}

	/*** Complete constructor ***
 	 * Sets all parameters on given values */
  	ChR2(double _tau_ChR, double _gamma_d1, double _gamma_d2, double _gamma_th, double _e_12d, double _e_21d, double _efficiency1, double _efficiency2,
		  const double _dt) : dt(_dt)
	{
		tau_ChR = _tau_ChR;
		gamma_d1 = _gamma_d1;
		gamma_d2 = _gamma_d2;
		gamma_th = _gamma_th;
		e_12d = _e_12d;
		e_21d = _e_21d;
		efficiency1 = _efficiency1;
		efficiency2 = _efficiency2;
		reset();
	}

	/*** Default constructor ***
	 * Sets all parameters on experimentally determined values */
  	ChR2(const double _dt) : dt(_dt)
	{
		/* Parameters from Williams et al., 2013, tab. 1 */
		tau_ChR = 1.3;
		ret_cross_section = 12.0;
		lambda_max = 470.0;
		wloss = 1.3;
		gamma_d1 = 75.0 + 43.0 * tanh(-(-70.0 + 20.0) / 20.0); // = 117.4 at membrane potential -70 mV
		gamma_d2 = 50.0;
		gamma_th = 4345.87e-5 * exp(-0.0211539274 * (-70.0)); // = 0.191 at membrane potential -70 mV
		e_12d = 11.0;
		e_21d = 8.0;
		efficiency1 = 0.8535;
		efficiency2 = 0.14;

		/* Parameters from FR multifit *
		gamma_d1 = 35.5232;
		gamma_d2 = 111.564;
		gamma_th = 20.1497;
		e_12d = 18.2932;
		e_21d = 4.49165e-05;
		efficiency1 = 0.351144;
		efficiency2 = 0.148877;*/

		/* Parameters from Nikolic et al., 2009 *
		tau_ChR = 0.2; // Nikolic et al., 2009, fig. 9
		ret_cross_section = 1.20; // Nikolic et al., 2009, p. 404
		lambda_max = 470.0; // Nikolic et al.
		wloss = 1.1; // Nikolic et al., 2009, fig. 4
		gamma_d1 = 110.0; // Nikolic et al., 2009, fig. 4
		gamma_d2 = 25.0; // Nikolic et al., 2009, fig. 4
		gamma_th = 4.0; // Williams et al. (2013), tab. 1 ("used by others")
		e_12d = 10.0; // Nikolic et al., 2009, fig. 4
		e_21d = 15.0; // Nikolic et al., 2009, fig. 4
		efficiency1 = 0.5; // Nikolic et al., 2009, fig. 4
		efficiency2 = 0.15; // Nikolic et al., 2009, fig. 4 */

		G_ChR2 = 100.0; // estimated (cf. Lin 2011)
		gamma = 0.1; // from Williams et al., 2013

		reset();
	}

};
