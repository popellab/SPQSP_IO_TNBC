# QSP FILES - SBML CONVERSION

This section first examines the Systems Biology Markup Language (SBML) files and converts them to C++ files which are compatible with SUNDIAL/CVODE package. 
All required files for SBML conversion are available at https://github.com/popellab/SPQSP_IO/tree/main/sbml_cvode and https://github.com/popellab/SPQSP_IO_TNBC (Wang_QSP_Model.xml and vct_simulation.xml). 
Nevertheless, in order to use the spatial QSP model, this step is not necessary since the SBML conversion has already been performed by the developer.

1.	Start program to at package directory. <br />
Replace <pkg_dir> with package directory.<br />
`$ cd <pkg_dir>`<br />
`$ python converter_gui.py`

2.	Import SBML model file (Wang_QSP_Model.xml).<br />

3.	Click “Load setting” button and load the default configuration file (vct_simulation.xml).<br />
a.	The configuration help to specify simulation time, relative, and absolute tolerance (ignore the hash mismatch error). <br />
b.	Recommend: Click “Validate units” and select “convert units”. <br />

4.	Deactivate “Use extra adjustable variables” and click “Analyze model” to apply the configuration to C++ files. Select the output directory and type in the Namespace (Cancer_VCT).<br />

5.	Export files (type ODE_system in the Export space) that are necessary for spQSP model, including:<br />
ODE_system.h/ODE_system.cpp: QSP model class header and implementation files (derived from class CVODEBase).<br />
Param.h/Param.cpp: model parameter classes (derived from class ParamBase).<br />
ODE_system_params.xml: model parameter file.  <br />

The generated outputs must be slightly modified to be able to run the spatial QSP model. Their final form, ready to use, can be found at https://github.com/popellab/SPQSP_IO_TNBC/tree/main/TNBC/SP_QSP_TNBC/ode and https://github.com/popellab/SPQSP_IO_TNBC/tree/main/TNBC/TNBC_sim/resource (param_all.xml).

-----------

# Spatial QSP SIMULATION

(The following instructions are similar to the guidelines written by Shuming Zhang for a different spatial QSP model: https://github.com/popellab/spQSP-omics-2021)

## Ubuntu Operating System Configuration (Only required for Windows User)

This step helps to setup the Ubuntu operating system in for Windows user via virtual machine.
1. The virtual machine host: VirtualBox is available at https://www.virtualbox.org/
2 . The Ubuntu Desktop image (Latest version 20.04.2) is available at   http://www.releases.ubuntu.com/20.04/
3. Enter the “Oracle VM VirtualBox Manager”, press “New” to create the virtual machine with all default settings. (Recommend allocate 20 GB for storage and 2 GB of RAM)

**Notice**: All following operations should be done in the Linux operating systems (the virtual machine), NOT Windows.

## Required Library Installation
Libraries: **SUNDIALS**: version:4.0.1; **Boost**: version 1.70.0

-SUNDIALS <br />
1. Download is available at: https://computing.llnl.gov/projects/sundials/sundials-software <br />
The following files are downloaded: 
 `sundials-4.0.1.tar.gz`

2. Decompress Archieve: <br />
`$ tar xzf sundials-4.0.1.tar.gz`

3. Install cmake if not already available: <br />
`$ sudo apt install cmake-curses-gui`

4.	Create install and build directories: <br />
`$ mkdir -p ~/lib/sundials-4.0.1`
`$ mkdir -p ~/Downloads/sundials-build`
`$ cd ~/Downloads/sundials-build`

5. Configuration <br />
`$ ccmake ~/Downloads/sundials-4.0.1` <br />
Press c key to enter configuration interface
Set install directory: CMAKE_INSTALL_PREFIX set to `~/lib/sundials-4.0.1`
Set example install directory: EXAMPLE_INSTALL_PATH set to `~/lib/sundials-4.0.1/examples`
Press c repeatedly to process configuration; press g to generate Makefile and exit configuration.
6. Build <br />
From `~/Downloads/sundials-build/` <br />
`$ make` <br />
`$ make install`

-Boost Version 1.70.0 <br />

1. Source code available at: https://www.boost.org/users/history/version_1_70_0.html <br />
The following files are downloaded: <br />
`boost_1_70_0.tar.gz`

2.	Decompress the archive: <br />
`$ tar xzf boost_1_70_0.tar.gz` <br />

Official instructions is available at:
https://www.boost.org/doc/libs/1_70_0/more/getting_started/unix-variants.html

3. Building separately-compiled boost libraries <br />
`$ cd ~/Downloads/boost_1_70_0` <br />
`$ ./bootstrap.sh --prefix=$HOME/lib/boost_1_70_0` <br />
`$ ./b2 install` <br />

## Model Simulation 
The Makefile of this model is available at: `~/TNBC_omics/TNBC_single/linux/` <br />
`$ make TNBC_s_sim` <br />
`$ ./TNBC_s_sim -h` <br />
For all options to configure the simulation
