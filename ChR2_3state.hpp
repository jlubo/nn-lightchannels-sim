/**************************************************************************
 * Three-state model of a single channelrhodopsin-2 channel
 **************************************************************************
 *
 * Copyright (C) Jannik Luboeinski 2020
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

#define THREE_STATE_CHR2
#define FIT_WITH_MAX
#define h 6.62606957e-5 // 10^{-29} J*s, Planck constant, 10^{-5} instead of 10^{-34} for canceling out cross section and lambda magnitude in computation
#define c 299792458.0 // m/s, speed of light

#include "Stimulus.hpp"

/*** Channelrhodopsin-2 class ***
 * Stochastically models one ChR2 channel molecule with three states 'closed', 'open' and
 * 'desensitized', which is initially in 'closed' state */
class ChR2
{
	friend class Neuron; // to allow Neuron class access to parameters

private:
	// state variables:
	double C; // probability for closed state
	double O; // probability for open state
	double D; // probability for desensitized state
	int t_light_begin; // time step at which light began to shine

protected:
	// parameters:
	double tau_ChR; // ms, activation time constant for the ion channel
	double ret_cross_section; // 10^{-20} m^2 - see remark on "h", retinal cross-section
	double G_ChR2; // fS, single channel conductance
	double lambda_max; // nm - see remark on "h", the wavelength of maximum absorption
	double wloss; // loss factor for channel environment
	double gamma_d0; // 1/s, desensitization rate (O -> D) for -70 mV
	double gamma_r0; // 1/s, recovery rate (D -> C)
	double efficiency; // probability that a photon absorbed by the retinal activates the channel

	// others:
	double dt; // ms, duration of one time step
	double photon_influx; // 1/s, effective photons per second the channel is exposed to


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
		*f << "lambda_max = " << lambda_max << " nm" << endl;
		*f << "wloss = " << wloss << endl;
		*f << "gamma_d0 = " << gamma_d0 << " 1/s" << endl;
		*f << "gamma_r0 = " << gamma_r0 << " 1/s" << endl;
		*f << "efficiency = " << efficiency << endl;
		*f << "photon_influx(E_max) = " << photon_influx << " 1/s" << endl;
	}

	/*** getO ***
	 * - return: the current value of the open state probability */
	double getO() const
	{
		return O;
	}

	/*** getCalcO ***
	 * Returns the probability to be in an open state, needed for calculation *
	 * - return: the current value of the open state probability */
	double getCalcO() const
	{
		return O;
	}

	/*** getC ***
	 * - return: the current value of the closed state probability */
	double getC() const
	{
		return (1-O-D);
	}

	/*** getD ***
	 * - return: the current value of the desensitized state probability */
	double getD() const
	{
		return D;
	}

	/*** getGammaD ***
	 * Computes and returns the voltage-dependent gamma_d parameter *
	 * - double mem_pot [optional]: the current membrane potential, if not set, -70 mV is used *
	 * - return: gamma_d in 1/s */
	double getGammaD(double mem_pot = -70.0) const
	{
		return gamma_d0 * (1.0 - 0.0056*(mem_pot + 70)); // from Tchumatchenko et al., 2013
	}

	/*** getGammaR ***
	 * Computes and returns the light-dependent gamma_r parameter as suggested by Nikolic et al. *
	 * - double mem_pot: the current membrane potential *
	 * - return: gamma_r in 1/s */
	double getGammaR() const
	{
		return gamma_r0; // light-assisted recovery might be implemented here (see Nikolic et al., 2009)
	}

	/*** getPhotonInflux ***
	 * Computes and returns the photon influx for a given irradiance
	 * - double irradiance: irradiance in mW/mm^2 *
	 * - return: netto photon influx in 1/s */
	double getPhotonInflux(double irradiance) const
	{
		return efficiency * (irradiance * 1000.0 * ret_cross_section * lambda_max) / (h * c * wloss);
  	}

	/*** setIrradiance ***
	 * Calculates the photon influx for a given irradiance (for not each step has to be calculated)
	 * - double irradiance: irradiance in mW/mm^2 */
	void setIrradiance(double irradiance)
	{
		photon_influx = efficiency * (irradiance * 1000.0 * ret_cross_section * lambda_max) / (h * c * wloss);
  	}

	/*** setPhotonInflux ***
	 * Sets the photon influx (for frequency response check only)
	 * - double _photon_influx: netto photon influx in 1/s */
	void setPhotonInflux(double _photon_influx)
	{
		photon_influx = _photon_influx;
  	}

   	/*** processTimeStep ***
	 * Processes one time step (of duration dt) for the channel's state probabilities
	 * - int tb_step: time step at which to evaluate stimulus (only of importance when LightStimulus is set) *
	 * - double mem_pot [optional]: the current membrane potential (needed to compute gamma_d) */
	void processTimeStep(int tb_step, double mem_pot = -70.0)
	{
		double gamma_d = getGammaD(mem_pot); // deactivation rate
		double gamma_r = getGammaR(); // recovery rate
		double p; // activation probability

		// Activation function computation, cf. Nikolic et al. 2009 Appendix 1
		if (photon_influx >= 1e-30) // light occurred
		{
			if (t_light_begin == 0)
				t_light_begin = tb_step;
			p = (1 - exp(- (tb_step-t_light_begin)*dt / tau_ChR)); // calculation in milliseconds (dt, tau_ChR in ms)
		}
		else // no light
		{
			t_light_begin = 0;
			p = 0.0;
		}

		double dO = (p*photon_influx - (p*photon_influx+gamma_d)*O - p*photon_influx*D) * dt/1000.0;
		double dD = (gamma_d * O - gamma_r * D) * dt/1000.0;

  		O = O + dO;
		D = D + dD;
  	}

	/*** reset ***
	 * Resets probabilities to initial state ('closed' state has probability 1) */
	void reset()
	{
		C = 1;
		O = 0;
		D = 0;
		photon_influx = 0.0;
		t_light_begin = 0;
	}

	/*** Default constructor ***
	 * Sets all parameters on experimentally determined values */
  	ChR2(const double _dt) : dt(_dt)
	{
		tau_ChR = 1.3; // Williams et al., 2013
		ret_cross_section = 12.0; // Williams et al., 2013
		G_ChR2 = 100.0; // estimated (cf. Lin 2011)
		lambda_max = 470.0; // Williams et al., 2013
		wloss = 1.3; // Williams et al., 2013
		gamma_d0 = 126.74; // Tchumatchenko et al., 2013, H134R mutant
		gamma_r0 = 8.38; // Tchumatchenko et al., 2013, H134R mutant
		/* wild-type ChR2:
		gamma_d0 = 236.35; // Tchumatchenko et al., 2013
		gamma_r0 = 3.60; // Tchumatchenko et al., 2013 */
		efficiency = 0.5; // Nikolic et al., 2009
		reset();
	}

};
