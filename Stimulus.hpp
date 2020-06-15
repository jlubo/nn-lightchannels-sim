/**************************************************************************
 * Template for defining a stimulus that is delivered to neurons
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

/*** Stimulus class ***
 * Models one period of a stimulus (of arbitrary dimension) - different pulse shapes can be added, *
 * when at one and the same time several pulses overlap, their magnitudes *
 * are simply summed up */
class Stimulus
{
	double* shape; // stimulus magnitude values for all time steps
	int n; // number of time steps of a period

	/*** freeShape ***
	 * Frees the memory reserved for the shape array, in case there is reserved memory */
	void freeShape()
	{
		if (shape != NULL)
		{
			delete[] shape;
			shape = NULL;
		}
	}

public:
	/*** printShape ***
	 * Writes stimulus shape data into a file with the "real" times and creates *
	 * a plot using a GNUplot script *
	 * - f: the file name without suffix
	 * - dt: duration of one time step in ms (for conversion to "real" time)
    * - return: whether suceeded or not */
	bool printShape(const char* name, double dt) const
	{
		char txt[64], gpl[64], pdf[64];
		ofstream f;

		// prepare file names
		strcpy(txt, name);
		strcat(txt, ".txt");
		strcpy(gpl, name);
		strcat(gpl, ".gpl");
		strcpy(pdf, name);
		strcat(pdf, ".pdf");

		// open data file
		f.open(txt);
		if (!f.is_open())
			return false;

		// write data file
		f.precision(3);
		for(int i=0; i<n; i++)
			f << fixed << double(i)*dt << "\t\t\t" << shape[i] << endl;
		f.close();

		// open script file
		f.open(gpl);
		if (!f.is_open())
			return false;

		// write script file
		f << "set term pdf enhanced font \"Sans, 9\" color solid lw 3" << endl;
		f << "set output '" << pdf << "'" << endl;
		f << "set xlabel \"t / ms\"" << endl;
		f << "set ylabel \"stim. dim.\"" << endl;
		f << "plot [x=" << -5*dt << ":" << (n+5)*dt << "][y=" << 0.1*max(shape, n) << ":" << 1.1*max(shape, n) << "] '"
		  << txt << "' using 1:2 notitle with lines" << endl;
		f.close();

		strcpy(txt, "gnuplot ");
		strcat(txt, gpl);
		system(txt);
		return true;
	}

	/*** getStimulusAt ***
	 * Returns the stimulus magnitude at a given time (assumes periodicity *
	 * of stimulus) *
	 * - t_step: time step at which to evaluate stimulus
    	 * - return: stimulus at given time */
	double getStimulusAt(int t_step) const
	{
		return shape[t_step % n];
	}

	/*** getLength ***
	 * Returns the stimulus length as number of time bins *
    	 * - return: stimulus length */
	int getLength() const
	{
		return n;
	}

	/*** addRectPulse ***
	 * Adds a rectangle pulse of a given length and at a given time *
	 * to the stimulus *
	 * - magnitude: height of the rectangle pulse
	 * - offset: time in steps (within a period) at which rectangle pulse shall start
	 * - len: time length of the rectangle pulse in steps
    	 * - return: true if successful */
	bool addRectPulse(double magnitude, int offset, int len)
	{
		if (offset >= n)
			return false;
		for(int i=0; i<(len > n-offset ? n-offset : len); i++)
		{
			shape[offset+i] += magnitude;
		}
		return true;
	}

	/*** addSinePulse ***
	 * Adds a sine pulse (a whole period) of a given period and at a given time *
	 * to the stimulus *
	 * - magnitude: amplitude of sine function
	 * - offset: time in steps (within a period) at which sinus pulse shall start
	 * - period: time length of one sinus period in steps
    	 * - return: true if successful */
	bool addSinePulse(double magnitude, int offset, int period)
	{
		if (offset >= n)
			return false;
		for(int i=0; i<(period > n-offset ? n-offset : period); i++)
		{
			shape[offset+i] += magnitude * sin( (2.0*M_PI/period) * i );
		}
		return true;
	}

	/*** multiplyBy ***
	 * Multiplies all shape magnitudes by a real number *
	 * - r: number to multiply with */
	void multiplyBy(double r)
	{
		for(int i=0; i<n; i++)
		{
			shape[i] *= r;
		}
	}

	/*** clear ***
	 * Clears all stimulus content */
	void clear()
	{
		for (int i = 0; i < n; i++)
		{
			shape[i] = 0.0;
		}
	}

	/*** Empty constructor ***
	 * Just for syntactic reasons */
  	Stimulus() : n(0)
	{
		shape = NULL; // to declare that shape has not been assigned memory
	}

	/*** Principal constructor ***
	 * Sets main characteristics of stimulus *
	 * - double n: total time step count of one period */
  	Stimulus(double _n) : n(_n)
	{
		if (n > 0)
			shape = new double[n];
		else
			shape = NULL;
		clear();
	}

	/*** Copy constructor ***
    	 * Reads shape and properties from a Stimulus variable *
   	 * and copies them into a new Stimulus (basically does *
	 * the same as assignment operator) *
	 * - Stimulus& st: reference to a variable of type Stimulus */
  	Stimulus(const Stimulus& _st) : n( _st.getLength() )
	{
		if (n > 0)
			shape = new double[n];
		else
			shape = NULL;
		for (int i = 0; i < n; i++)
		{
			shape[i] = _st.getStimulusAt(i);
		}
	}

	/*** Assignment operator ***
	 * Reads shape and properties from a Stimulus variable *
    	 * and copies them into a new Stimulus (basically does *
	 * the same as copy constructor) *
	 * - Stimulus& st: reference to a variable of type Stimulus */
	Stimulus& operator=(const Stimulus& _st)
	{
		if (this == &_st) return *this;

		n = _st.getLength();
		freeShape();

		if (n > 0)
			shape = new double[n];
		for (int i = 0; i < n; i++)
		{
			shape[i] = _st.getStimulusAt(i);
		}

		return *this;
	}


	/*** Destructor ***
	 * Frees reserved memory */
  	~Stimulus()
	{
		freeShape();
	}

};
