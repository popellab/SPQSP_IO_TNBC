#include "Param.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <math.h>
#include "../ode/Param.h"
#include "../ode/ODE_system.h"

// parameters from QSP module. 
extern CancerVCT::Param qsp_params;

#define QP(x) CancerVCT::ODE_system::get_class_param(x)
#define AVOGADROS 6.022140857E23 

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

namespace pt = boost::property_tree;

static int SEC_PER_DAY = 86400;
static int HOUR_PER_DAY = 24;

// must match the order in the ParamInt and ParamFloat enums
#define PARAM_DESCRIPTION_FIELD_COUNT 3
const char* _description[][PARAM_DESCRIPTION_FIELD_COUNT] =
{
	//{"fullpath", "desc", "constraint"}

	//------------------------ float --------------------------//
	/* QSP */
	{ "Param.QSP.simulation.weight_qsp", "", "prob" },
	{ "Param.QSP.simulation.t_steadystate", "days", "pos" },
	{ "Param.QSP.simulation.t_resection", "days", "pos" },
	/* ABM */
	//environmental
	{ "Param.ABM.Environment.SecPerSlice", "", "pos" },
	{ "Param.ABM.Environment.ScalingFactor", "scaling factor", "pos" },
	{ "Param.ABM.Environment.Sources", "sources", "pos" },
	{ "Param.ABM.Environment.recSiteFactor", "number of adhesion site per port voxel", "pos" },
	{ "Param.ABM.Environment.adhSiteDens", "total adhesion site density on tumor vasculature, mol/m^3", "pos" },
	//pharmacokinetics	
	{ "Param.ABM.Pharmacokinetics.nivoDoseIntervalTime", "days", "pos" },
	{ "Param.ABM.Pharmacokinetics.nivoDose", "mole/m^3", "pos" },
	{ "Param.ABM.Pharmacokinetics.durvDoseIntervalTime", "days", "pos" },
	{ "Param.ABM.Pharmacokinetics.durvDose", "mole/m^3", "pos" },	
	{ "Param.ABM.Pharmacokinetics.entDoseIntervalTime", "days", "pos" },
	{ "Param.ABM.Pharmacokinetics.entDose", "mole/m^3", "pos" },	
	{ "Param.ABM.Pharmacokinetics.ipiDoseIntervalTime", "days", "pos" },
	{ "Param.ABM.Pharmacokinetics.ipiDose", "mole/m^3", "pos" },	
	//T cell
	{ "Param.ABM.TCell.lifespanMean", "days", "pos" },
	{ "Param.ABM.TCell.lifespanSD", "days", "pos" },
	{ "Param.ABM.TCell.moveProb", "", "pr" },
	{ "Param.ABM.TCell.IL2_release_time", "amount of time to release IL2 after stimulation, sec", "pos" },
	{ "Param.ABM.TCell.IL2_prolif_th", "accumulative IL2 exposure to proliferate, sec*ng/mL", "pos" },
	{ "Param.ABM.TCell.IFNg_release_time", "amount of time to release IFNg after stimulation, sec", "pos" },
	// Treg
	{ "Param.ABM.Treg.moveProb", "", "pr" },
	// MDSC
	{ "Param.ABM.MDSC.moveProb", "", "pr" },	
	//cancer cell
	{"Param.ABM.CancerCell.asymmetricDivProb", "", "pr"},
	{"Param.ABM.CancerCell.senescentDeathRate", "per day", "pos"},
	{"Param.ABM.CancerCell.moveProb_csc", "", "pr"},
	{"Param.ABM.CancerCell.moveProb", "", "pr"},
	{"Param.ABM.CancerCell.mincc", "", "pos"},
	{"Param.ABM.CancerCell.IFNgUptake", "per sec", "pos"},
	//agent chemokine interaction
	{"Param.ABM.cell.PDL1_th", "percent of max PDL1_syn to be detectable", "prob"},
	{"Param.ABM.cell.IFNg_PDL1_half", "c of IFNg to induce PDL1 to half maximal level, ng/mL", "pos"},
	{"Param.ABM.cell.IFNg_PDL1_n", "hill coef for PDL1 expression", "pos"},
	{"Param.ABM.cell.PDL1_halflife", "halflife of PDL1 expression, days", "pos"},
	/* molecular level */
	// diffusion grid
	{"Param.Molecular.biofvm.IFNg.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.IFNg.release","ng/sec, one cell", "pos"},
	{"Param.Molecular.biofvm.IFNg.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.IL_2.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.IL_2.release","ng/sec, one cell", "pos"},
	{"Param.Molecular.biofvm.IL_2.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.CCL2.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.CCL2.release","ng/sec, one cell", "pos"},
	{"Param.Molecular.biofvm.CCL2.decayRate","1/sec", "pos"},	
	{"Param.Molecular.biofvm.CCL2.molecularWeight","kDa", "pos"},		
	{"Param.Molecular.biofvm.ArgI.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.ArgI.release","ng/sec, one cell", "pos"},
	{"Param.Molecular.biofvm.ArgI.decayRate","1/sec", "pos"},	
	{"Param.Molecular.biofvm.ArgI.molecularWeight","kDa", "pos"},
	{"Param.Molecular.biofvm.NO.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.NO.release","ng/sec, one cell", "pos"},
	{"Param.Molecular.biofvm.NO.decayRate","1/sec", "pos"},	
	{"Param.Molecular.biofvm.NO.molecularWeight","kDa", "pos"},

	//------------------------ int ----------------------------//
	{"Param.ABM.Environment.Tumor.XSize", "", "pos"},
	{"Param.ABM.Environment.Tumor.YSize", "", "pos"},
	{"Param.ABM.Environment.Tumor.ZSize", "", "pos"},
	{"Param.ABM.Environment.Tumor.VoxelSize", "voxel resolution, microns", "pos"},
	{"Param.ABM.Environment.Tumor.nr_T_voxel", "", "pos"},
	{"Param.ABM.Environment.Tumor.nr_T_voxel_C", "", "pos"},
	{"Param.ABM.Environment.ShuffleInterval", "", "pos"},
	{"Param.ABM.Environment.gridshiftInterval", "", "pos"},
	{"Param.ABM.TCell.div_interval", "", "pos"},
	{"Param.ABM.TCell.div_limit", "", "pos"},	
	{"Param.ABM.Treg.div_limit_exp", "", "pos"},		
	{"Param.ABM.CancerCell.progenitorDivMax", "", "pos"},
	{"Param.Molecular.stepPerSlice","", "pos"},

	// ---------------------- bool -----------------------------//
	{"Param.QSP.simulation.use_resection", "", "" },
	{"Param.Molecular.allMolecularOff", "", ""},
	{"Param.Molecular.diffusionOff", "", ""},
	{"Param.ABM.Pharmacokinetics.nivoOn", "", "" },
	{"Param.ABM.Pharmacokinetics.durvOn", "", "" },
	{"Param.ABM.Pharmacokinetics.entOn", "", "" },
	{"Param.ABM.Pharmacokinetics.ipiOn", "", "" },

};

Param::Param()
	:ParamBase()
{
	setupParam();
}

/*! Setup parameter storage
	instantiation of pure virtual member of the base class.
	setup description vector;
	initialize parameter value vectors with 0/false, 
	with size determined by enums, 
	so that other base class members can access vector sizes
*/
void Param::setupParam(){

	size_t nrExternalParam = PARAM_FLOAT_COUNT + PARAM_INT_COUNT + PARAM_BOOL_COUNT;
	for (size_t i = 0; i < nrExternalParam; i++)
	{
		_paramDesc.push_back(std::vector<std::string>(_description[i], 
			_description[i]+ PARAM_DESCRIPTION_FIELD_COUNT));
	}
	_paramFloat = std::vector<double>(PARAM_FLOAT_COUNT, 0);
	_paramInt= std::vector<int>(PARAM_INT_COUNT, 0);
	_paramBool= std::vector<bool>(PARAM_BOOL_COUNT, false);
	_paramFloatInternal = std::vector<double>(PARAM_FLOAT_INTERNAL_COUNT, 0);
	_paramIntInternal = std::vector<int>(PARAM_INT_INTERNAL_COUNT, 0);
	_paramBoolInternal = std::vector<bool>(PARAM_BOOL_INTERNAL_COUNT, false);
}

/*! Calculate internal parameters
*/
void Param::processInternalParams(){
	
	_paramFloatInternal[PARAM_AVOGADROS] = AVOGADROS;

	//micrometer to cm
	_paramFloatInternal[PARAM_VOXEL_SIZE_CM] = _paramInt[PARAM_VOXEL_SIZE] / 1e4;

	_paramFloatInternal[PARAM_T_CELL_LIFE_MEAN_SLICE] = _paramFloat[PARAM_T_CELL_LIFE_MEAN]
		/ _paramFloat[PARAM_SEC_PER_TIME_SLICE] * SEC_PER_DAY;

	_paramFloatInternal[PARAM_T_CELL_LIFE_SD_SLICE] = _paramFloat[PARAM_T_CELL_LIFE_SD]
		/ _paramFloat[PARAM_SEC_PER_TIME_SLICE] * SEC_PER_DAY;

	_paramFloatInternal[PARAM_PDL1_DECAY_SLICE] = std::exp(-_paramFloat[PARAM_PDL1_DECAY_DAY]
		* _paramFloat[PARAM_SEC_PER_TIME_SLICE] / SEC_PER_DAY);		

	_paramBoolInternal[PARAM_MOLECULAR_MODULES_ON]
		= !getVal(PARAM_ALL_MOLECULAR_OFF);

	_paramBoolInternal[PARAM_DIFFUSION_ON] 
		= getVal(PARAM_MOLECULAR_MODULES_ON) && !getVal(PARAM_DIFFUSION_OFF);
	
	_paramBoolInternal[PARAM_IFN_SINK_ON] 
		= getVal(PARAM_DIFFUSION_ON) && getVal(PARAM_IFN_G_UPTAKE)>0;

}

//! update from QSP parameters
void Param::update_from_qsp(void){

	// cell
	_paramFloatInternal[PARAM_CELL] = QP(9)*AVOGADROS;	

	// minimum tumor volume
	_paramFloatInternal[PARAM_VT_MIN] = QP(13);

	// volume of cancer cells
	_paramFloatInternal[PARAM_V_CELL] = QP(11);

	// volume of T cells
	_paramFloatInternal[PARAM_V_TCELL] = QP(12);

	// maximum concentration of cancer cells
	_paramFloatInternal[PARAM_CMAX] = QP(15);

	// maximum concentration of Tregs per volume in the tumor
	_paramFloatInternal[PARAM_TREGMAX] = QP(178);		

	// maximum concentration of MDSC per volume
	_paramFloatInternal[PARAM_MDSCMAX] = QP(177);

	// MDSC base recruitment
	_paramFloatInternal[PARAM_K_BASE_REC] = QP(180);

	// MDSC recruitment by CCL2
	_paramFloatInternal[PARAM_K_REC_BY_CCL2] = QP(162);		

	// half maximal effective concentration of CCL2 on recruitment of MDSC into the tumor (mol/m^3)
	_paramFloatInternal[PARAM_EC50_CCL2_REC] = QP(175);	

	// half maximal inhibitory concentration of NO on inhibition of CD8+ T cell cytotoxic activity (ng/ml)
	//_paramFloatInternal[PARAM_IC50_NO_CTL] = QP(174)* _paramFloat[PARAM_NO_MOLECULAR_WEIGHT] * 1e6;
	_paramFloatInternal[PARAM_IC50_NO_CTL] = QP(174);

	// half maximal inhibitory concentration of Arg I on inhibition of CD8+ T cell cytotoxic activity (ng/ml)
	//_paramFloatInternal[PARAM_IC50_ARGI_CTL] = QP(173)* _paramFloat[PARAM_ARGI_MOLECULAR_WEIGHT] * 1e6;	
	_paramFloatInternal[PARAM_IC50_ARGI_CTL] = QP(173);

	// half maximal effective concentration of arginase I on Treg expansion (ng/ml)
	//_paramFloatInternal[PARAM_EC50_ARGI_TREG] = QP(176) * _paramFloat[PARAM_ARGI_MOLECULAR_WEIGHT] * 1e6;	
	_paramFloatInternal[PARAM_EC50_ARGI_TREG] = QP(176);
	
	// half maximal inhibitory concentration of entinostat on anti-proliferation of tumor cell (mol/m^3)
	_paramFloatInternal[PARAM_IC50_ENT_C] = QP(164);	

	// half maximal inhibitory concentration of entinostat on inhibition of CCL2 production (mol/m^3)
	_paramFloatInternal[PARAM_IC50_ENT_CCL2] = QP(179);	

	// half maximal inhibitory concentration of entinostat on inhibition of NO production (mol/m^3)
	_paramFloatInternal[PARAM_IC50_ENT_NO] = QP(171);

	// half maximal inhibitory concentration of entinostat on inhibition of Arg I production (mol/m^3)
	_paramFloatInternal[PARAM_IC50_ENT_ARGI] = QP(181);	

	// number of PD1/PDL1 binding for half maximal inhibition 
	_paramFloatInternal[PARAM_PD1_PDL1_HALF] = QP(142);

	// total number of PD1 per synapse
	_paramFloatInternal[PARAM_PD1_SYN] = QP(146)*QP(71)/QP(72);

	// total number of PDL1 per synapse
	_paramFloatInternal[PARAM_PDL1_SYN_MAX] = QP(150)*QP(71)/QP(73);
	
	// k1 for PDL1-PD1 calculation
	_paramFloatInternal[PARAM_PDL1_K1] = 1 / (QP(129)/QP(120)) / QP(71);

	// k2 for PDL1-PD1 calculation
	_paramFloatInternal[PARAM_PDL1_K2] = 1 / (QP(131)/QP(121)) / (QP(100)/2);

	// k3 for PDL1-PD1 calculation
	_paramFloatInternal[PARAM_PDL1_K3] = QP(139) * QP(121) / (2 * QP(131));
	// hill coefficient
	_paramFloatInternal[PARAM_N_PD1_PDL1] = QP(143);

	/*
	std::cout << "k1, k2, k3, T1, PDL1_tot" 
		<< ": " << _paramFloatInternal[PARAM_PDL1_K1]
		<< ", " << _paramFloatInternal[PARAM_PDL1_K2]
		<< ", " << _paramFloatInternal[PARAM_PDL1_K3]
		<< ", " << _paramFloatInternal[PARAM_PD1_SYN]
		<< ", " << _paramFloatInternal[PARAM_PDL1_SYN_MAX]
		<< std::endl;
	std::cout << "k50: " << _paramFloatInternal[PARAM_PD1_PDL1_HALF]<< std::endl;
	*/

	// Parameters calculated from QSP parameter values
	double t_step_sec = _paramFloat[PARAM_SEC_PER_TIME_SLICE];

	// T cell killing of Cancer cell
	// QP(40): k_C_death_by_T (day^-1, sec^-1 internal)
	_paramFloatInternal[PARAM_ESCAPE_BASE] = std::exp(-t_step_sec * QP(49));
	
	// T cell exhaustion from PDL1
	_paramFloatInternal[PARAM_EXHUAST_BASE_PDL1] = std::exp(-t_step_sec * QP(48));

	// T cell exhaustion from Treg inhibition 
	_paramFloatInternal[PARAM_EXHUAST_BASE_TREG] = std::exp(-t_step_sec * QP(37));

	// time for resection
	_paramFloatInternal[PARAM_RESECT_TIME_STEP] = _paramFloat[PARAM_QSP_T_RESECTION] * SEC_PER_DAY / t_step_sec;

	// Recruitment

	//The number of adhesion site per voxel is:
	double site_per_voxel = _paramFloat[PARAM_ADH_SITE_DENSITY] * std::pow(double(_paramInt[PARAM_VOXEL_SIZE]) / 1e6, 3) * AVOGADROS;
	//The number of adhesion site represented by each ABM recruitment port
	double site_per_port = _paramFloat[PARAM_REC_SITE_FACTOR];
	//percentage of voxels to be assigned as ports
	_paramFloatInternal[PARAM_REC_PORT_PROB] = 1;

	/*When calculating recruitment probability:
	*/
	double  w = _paramFloat[PARAM_WEIGHT_QSP];
	// Teff -> k (1/mol) // p = k (1/mol) * Cent.T (mol)
	_paramFloatInternal[PARAM_TEFF_RECRUIT_K] = QP(46) * t_step_sec; 
	// Treg -> k (1/mol) // p = k (1/mol) * Cent.T (mol)
	_paramFloatInternal[PARAM_TREG_RECRUIT_K] = QP(27) * t_step_sec;

	std::cout << "keff" << _paramFloatInternal[PARAM_TEFF_RECRUIT_K] << std::endl;
	std::cout << "kreg" << _paramFloatInternal[PARAM_TREG_RECRUIT_K] << std::endl;
	// MDSC -> k (m^3/mol) // p = k (m^3/mol) * (Tum.MDSCmax * Tum.Vol - Tum.MDSC) (mol) / Tum.Vol(m^3)
	_paramFloatInternal[PARAM_MDSC_RECRUIT_K] = QP(180) * t_step_sec;
	_paramFloatInternal[PARAM_MDSC_RECRUIT_BY_CCL2_K] = QP(162) * t_step_sec;


	// mean life of Treg, unit: time step
	_paramFloatInternal[PARAM_TREG_LIFE_MEAN] = 1 / (QP(24) * t_step_sec);
	// mean life of MDSC, unit: time step
	_paramFloatInternal[PARAM_MDSC_LIFE_MEAN] = 1 / (QP(163) * t_step_sec);
	//std::cout << "Internal param: " << QP(144)  <<  ", "  << getVal(PARAM_TREG_LIFE_MEAN) << std::endl;
	/*
	std::cout
	<< t_step_sec << ", " << QP(40) << ", " << QP(148) << "\n"
	<< "PARAM_ESCAPE_BASE, " << _paramFloatInternal[PARAM_ESCAPE_BASE] << "\n"
	<< "PARAM_EXHUAST_BASE_PDL1, " << _paramFloatInternal[PARAM_EXHUAST_BASE_PDL1] << "\n"
	<< "PARAM_EXHUAST_BASE_TREG, " << _paramFloatInternal[PARAM_EXHUAST_BASE_TREG] << "\n"
	<< std::endl;
	*/

	/*Treg expansion*/
	_paramFloatInternal[PARAM_TREG_EXP_INTERVAL_SLICE]
		= std::log(2) / QP(172) / t_step_sec;

	/*Cancer cell dynamics parameters*/

	// stem cell division rate is calculated from QSP parameter
	// unit: s^-1
	double rs = QP(14) / (1 - _paramFloat[PARAM_CANCER_STEM_ASYMMETRIC_DIV_PROB]);
	// unit: day^-1
	_paramFloatInternal[PARAM_CSC_GROWTH_RATE] = rs * SEC_PER_DAY;

	//Entire tumor (tumor growth 0.01)
	//_paramFloatInternal[PARAM_CANCER_PROG_GROWTH_RATE] = 0.005;

	//Partial tumor
	_paramFloatInternal[PARAM_CANCER_PROG_GROWTH_RATE] = QP(14);

	_paramFloatInternal[PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE]
		= ((1 - _paramFloat[PARAM_CANCER_STEM_ASYMMETRIC_DIV_PROB])*std::log(2)+_paramFloat[PARAM_CANCER_STEM_ASYMMETRIC_DIV_PROB])/rs / getVal(PARAM_SEC_PER_TIME_SLICE);

		std::cout << "stem" << _paramFloatInternal[PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE] << std::endl;

	_paramFloatInternal[PARAM_CANCER_SENESCENT_MEAN_LIFE] =
		1 / _paramFloat[PARAM_CANCER_SENESCENT_DEATH_RATE]
		/ _paramFloat[PARAM_SEC_PER_TIME_SLICE] * SEC_PER_DAY;
		

	_paramFloatInternal[PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE]
		= std::log(2)/ _paramFloatInternal[PARAM_CANCER_PROG_GROWTH_RATE]
		* SEC_PER_DAY / getVal(PARAM_SEC_PER_TIME_SLICE);

		std::cout << "prog" << _paramFloatInternal[PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE] << std::endl;

	/*
	std::cout << getVal(PARAM_FLOAT_CANCER_SENESCENT_DEATH_PROB)
		<< ", " << getVal(PARAM_INT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE) 
		<< ", " << rs << ": " << getVal(PARAM_INT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE) << std::endl;
	std::cout << getVal(PARAM_INT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE)
		<< ", "<< getVal(PARAM_CANCER_SENESCENT_MEAN_LIFE)
		<< ", "<< getVal(PARAM_INT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE)
		<< std::endl;
	*/

	return;
}

};
};