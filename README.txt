SIMULATION OF NEURAL NETWORKS AND SINGLE NEURONS STIMULATED VIA LIGHT-SENSITIVE CHANNELS
========================================================================================

Copyright (C) Jannik Luboeinski 2015-2020
GSL v2.4: Copyright (C) Mark Galassi, Jim Davies, James Theiler, Brian Gough, Reid Priedhorsky, Gerard Jungman, Michael Booth, Fabrice Rossi, Simone Piccardi, Carlo Perassi, Ho-Jin Dan, Szymon Jaroszewicz, Nicolas Darnis, Tuomo Keskitalo, Ivo Alxneit, Jason H. Stover, Patrick Alken, Rhys Ulerich, Pavel Holoborodko, Pedro Gonnet 2017 

E-mail: jannik.lubo[at]gmx.de

================================================================================================================
LICENSE:
========

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along with this program (in "LICENSE.txt").
If not, see <https://www.gnu.org/licenses/>.


================================================================================================================
MAIN REMARKS:
=============

The package provided here serves to run customizable simulations of neural networks and single neurons receiving
stimulation through light-activated channels, and to measure the spatial firing rate response in networks
as well as the firing rate and channel dynamics of single neurons under pulsed light stimulation. However, the code 
offers even more features, some of which will be mentioned in the following. Please feel free to contact the authors
on questions about this code and on further investigations that can be done with it.

In addition to the code, this package also contains pre-compiled, ready-to-run binaries for Linux and Windows
(Ubuntu 18.04 LTS and Windows 10), linked with GNU Scientific Library 2.4, which can be found in "bin/*".

This software has originally been developed in the scope of the following study:

 Luboeinski and Tchumatchenko, Network Neuroscience, 2020

Please consider to cite this work if you use the code or binaries for scientific purposes. 


================================================================================================================
CODE REMARKS:
=============

Features of the model can be adapted in a modular way as described in "CUSTOM_MODELS.txt".

The code has been developed for Linux, and also for Windows operating systems. Although it is mostly platform-independent, 
file system operations may have to be adapted if you use another O/S. Also, if you want to use GNU Scientific Library
(which the code employs by default), you need to have the library installed.

Fitting of a Gaussian curve to the spatial firing rate distribution of the excitatory population is implemented
using the GNU Scientific Library. Fitting results should always be checked for convergence (see "gslfit.log"
and "*_gauss_sigma_fit_*.pdf" files). If gnuplot is installed, an additional gnuplot fit will automatically be
created in "*_gauss_sigma_fit_*.pdf" files for those light intensities defined in the "stimplot" variable in
"NetworkBatch.cpp".

Fitting of a Gaussian curve to the spatial firing rate distribution of the inhibitory population has only
partially been implemented (only important for local connectivities).

Adaptation time constants of the firing rate of single neurons and the open state probability of channels can
automatically be determined from the traces. The fitting method is implemented using GNU Scientific Library.
The results depend strongly on the initial conditions and should therefore be considered with caution.

The code contains a stimulus comparison mode to directly perform a comparison between the impact of electric
current stimulation and the impact of light stimulation (to be switched on using '#define STIMULUS_COMPARISON').
The implementation of this mode, however, is still incomplete and has not been tested properly.

Finally, a short remark on the purpose of the main code files: "NetworkSimulation.cpp"/"NeuronSimulation.cpp" contain
most of the code necessary for simulations with varying light intensities, while "NetworkBatch.cpp"/"NeuronBatch.cpp"
are used to loop over varying stimulus frequencies, network connectivities or coupling strengths and to process
command-line arguments. Furthemore, "NetworkBatch.cpp" automatically adjusts the mean external input current for
certain parameter settings to achieve a mean firing rate of 5 Hz without stimulation.

Files:

NetworkBatch.cpp: main functions for a batch of network simulations with different stimulus parameters
NeuronBatch.cpp: main functions for a batch of single-neuron simulations with different stimulus parameters
NetworkSimulation.cpp: class for performing network simulations
NetworkSimulation.cpp: class for performing single-neuron simulations
Network.cpp: class describing a network
Neuron.cpp: class describing a neuron
ChR2_3state.hpp: class describing a channel probabilistically using a three-state model
ChR2_4state.hpp: class describing a channel probabilistically using a four-state model (optionally)
Plots.hpp: collection of plotting functions
Stimulus.hpp: class for defining a stimulus protocol
Tool.hpp: collection of utility functions


================================================================================================================
COMPILE AND LINK IN LINUX:
==========================

Requirements:

 g++ needs to be installed (tested with g++ version 4.6.3)

 GNU Scientific Library (tested with GSL version 2.4) needs to be installed; alternatively, GSL functions in
 "Tools.hpp" can be commented and gnuplot or some other program be used for fitting


For network simulations, run:

 g++ -std=c++0x -D LINUX "NetworkBatch.cpp" -o "net.out" -lgsl -lgslcblas


For single-neuron simulations, run:

 g++ -std=c++0x -D LINUX "NeuronBatch.cpp" -o "sn.out" -lgsl -lgslcblas


To include the GSL static libraries into the binary, add:

 -static

See also the bash script "compile.sh". Pre-compiled binaries for Linux (tested with Ubuntu 18.04 LTS) can be found in the directory "bin/".


================================================================================================================
COMPILE AND LINK IN WINDOWS:
============================

Requirements:

 g++ (MinGW) needs to be installed (tested with g++ version 4.6.3) and "[MinGW]\bin" bound to the PATH environment variable

 GNU Scientific Library (tested with GSL version 2.4) needs to be installed in the MinGW directory, for example using MSYS; alternatively, 
 GSL functions in "Tools.hpp" can be commented and gnuplot or some other program be used for fitting


For network simulations, run:

 g++ -std=c++0x -D WINDOWS "NetworkBatch.cpp" -o "net.exe" -lgsl -lgslcblas


For single-neuron simulations, run:

 g++ -std=c++0x -D WINDOWS "NeuronBatch.cpp" -o "sn.exe" -lgsl -lgslcblas


To include the GSL static libraries into the binary, add:

 -static

See also the bash script "compile.bat". Pre-compiled binaries for Linux (tested with Ubuntu 18.04 LTS) can be found in the directory "bin/".


================================================================================================================
RUN SIMULATIONS IN LINUX:
=========================

If you want plots to be created automatically, which is recommended, please make sure that gnuplot (tested with version 5.0.1) 
is installed and bound to the system command "gnuplot".

For network simulations, enter the shell and run:

 ./net.out


For single-neuron simulations, enter the shell and run:

 ./sn.out


Command-line arguments can be passed to set parameter values, for example:

 ./net.out -pc=0.01 -J_start=0.2 -J_end=2.0

All possible command-line arguments can be found in "NetworkBatch.cpp" and "NeuronBatch.cpp".


================================================================================================================
RUN SIMULATIONS IN WINDOWS:
===========================

If you want plots to be created automatically, which is recommended, please make sure that gnuplot (tested with version 5.0.1)
is installed and the gnuplot path is bound to the PATH environment variable.

For network simulations, simply run:

 net.exe


For single-neuron simulations, simply run:

 sn.exe


Command-line arguments can be passed to set parameter values, for example:

 net.exe -pc=0.01 -J_start=0.2 -J_end=2.0

All possible command-line arguments can be found in "NetworkBatch.cpp" and "NeuronBatch.cpp".
