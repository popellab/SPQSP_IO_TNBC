#include "LymphCentral.h"
#include "Tumor.h"
#include "../../core/GlobalUtilities.h"

// shorthands
// get raw value (original units)
#define GET_PARAM_RAW(x) _QSP_model.getSystem()->getParameterVal(x, true)
#define SET_PARAM_RAW(x, y) _QSP_model.getSystem()->setParameterVal(x, y, true)
#define GET_VAR_RAW(x) _QSP_model.getSystem()->getSpeciesVar(x, true)
#define SET_VAR_RAW(x, y) _QSP_model.getSystem()->setSpeciesVar(x, y, true)
// get value (SI units)
#define GET_PARAM(x) _QSP_model.getSystem()->getParameterVal(x, false)
#define SET_PARAM(x, y) _QSP_model.getSystem()->setParameterVal(x, y, false)
#define GET_VAR(x) _QSP_model.getSystem()->getSpeciesVar(x, false)
#define SET_VAR(x, y) _QSP_model.getSystem()->setSpeciesVar(x, y, false)
// parameter (SI units)
#define QSP_CONST(x) LymphBloodQSP::get_class_param(x)

// indices of parameter/variables in their vectors
// y
#define QSP_ID_TUM_C1 20
#define QSP_ID_CENT_TEFF 1
#define QSP_ID_CENT_TREG 0
#define QSP_ID_TUM_TEFF 22
#define QSP_ID_TUM_TREG 21
#define QSP_ID_CENT_NIVO 2
#define QSP_ID_CENT_DURV 3
#define QSP_ID_CENT_IPI 4
#define QSP_ID_CENT_ENT 5
#define QSP_ID_TUM_MDSC 32
#define QSP_ID_TUM_NIVO 26
#define QSP_ID_TUM_ENT 36
#define QSP_ID_PD1_PDL1 62
#define QSP_ID_P0 46
#define QSP_ID_C  25
#define QSP_ID_P1 47
#define QSP_ID_TUM_CX 18
#define QSP_ID_TUM_TEXH 19
#define QSP_ID_TUM_CCL2 33
#define QSP_ID_TUM_ARGI 35
#define QSP_ID_TUM_NO 34

// class_param
#define QSP_C_MAX 15
#define QSP_CELL 9
#define QSP_P1_C1 85
#define QSP_P0_C1 70
#define QSP_DAMPS 59
#define QSP_N_T0_CLONES 19
#define QSP_N_T1_CLONES 38
#define QSP_VT_MIN 13
#define QSP_INIT_TUM_DIAM 18
#define QSP_VOL_CELL 11
#define QSP_VOL_TCELL 12

// constants
#define AVOGADROS 6.022140857E23 
#define SEC_PER_DAY 86400
#define PI 3.1416

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

LymphCentral::LymphCentral()
	: _QSP_model()
	, _var_qsp_to_abm()
{
	// Cent.Teff, Cent.Treg, Cent.Nivo
	_var_qsp_to_abm = std::vector<double>(QSPEX_VAR_NUM, 0);
}

LymphCentral::~LymphCentral()
{
}

void LymphCentral::setup_param(LymphBloodParam& p){

	bool steadystate = true;
	unsigned int n = _QSP_model.getSystem()->get_num_variables();
	std::vector<double> ss_val(n, 0);

	if (steadystate)
	{
		QSP ss;
		LymphBloodQSP::_QSP_weight = 1;
		LymphBloodQSP::use_steady_state = true;
		LymphBloodQSP::use_resection = false;
		LymphBloodQSP::setup_class_parameters(p);
		ss.getSystem()->setup_instance_tolerance(p);
		ss.getSystem()->setup_instance_varaibles(p);
		ss.getSystem()->eval_init_assignment();

		// run to steady state or until volume condition is met
		double tss = params.getVal(PARAM_QSP_STEADYSTATE) * SEC_PER_DAY;
		double tt = 0;
		double deltatt = params.getVal(PARAM_SEC_PER_TIME_SLICE);
		double tumor_volume = QSP_CONST(QSP_VT_MIN) + (QSP_CONST(QSP_VOL_CELL) * (ss_val[20] + ss_val[18]) + QSP_CONST(QSP_VOL_TCELL) * (ss_val[19] + ss_val[21] + ss_val[22])) / AVOGADROS;
		double tumor_volume_ref = (PI * std::pow(QSP_CONST(QSP_INIT_TUM_DIAM),3)) / 6;

		//Total tumor (tumor growth): params.getVal(PARAM_QSP_INIT_CANCER_CELL) = 1

		//Total tumor: params.getVal(PARAM_QSP_INIT_CANCER_CELL) = 8.75e7

		//Partial tumor
		while (ss_val[20]<params.getVal(PARAM_QSP_INIT_CANCER_CELL)) 		

		//QSP solutions
		//while (tt < tss && tumor_volume < tumor_volume_ref)
		{
			ss.solve(tt, deltatt);
			for (size_t i = 0; i < n; i++)
			{
				ss_val[i] = ss.getSystem()->getSpeciesVar(i);
			}
			tumor_volume = QSP_CONST(QSP_VT_MIN) + (QSP_CONST(QSP_VOL_CELL) * (ss_val[20] + ss_val[18]) + QSP_CONST(QSP_VOL_TCELL) * (ss_val[19] + ss_val[21] + ss_val[22])) / AVOGADROS;
			tt += deltatt;
		}
		//if (tumor_volume < tumor_volume_ref)
		//{		
		//	std::cout<<"tumor volume condition is not met"<<std::endl;
		//	exit(0);
		//}	
	}

	// setup

	LymphBloodQSP::_QSP_weight = params.getVal(PARAM_WEIGHT_QSP);
	LymphBloodQSP::use_steady_state = false;
	LymphBloodQSP::setup_class_parameters(p);
	_QSP_model.getSystem()->setup_instance_tolerance(p);
	_QSP_model.getSystem()->setup_instance_varaibles(p);

	// load steady state
	if (steadystate)
	{
		for (size_t i = 0; i < n; i++)
		{
			_QSP_model.getSystem()->setSpeciesVar(i, ss_val[i]);
		}
	}

	_QSP_model.getSystem()->eval_init_assignment();
	_QSP_model.getSystem()->updateVar();
}

/*! solve QSP from t to t + dt
*/
void LymphCentral::time_step(double t, double dt){

	// Pharmacokinetics
		
	double tumor_volume = QSP_CONST(QSP_VT_MIN) + (QSP_CONST(QSP_VOL_CELL) * (GET_VAR_RAW(QSP_ID_TUM_C1) + GET_VAR_RAW(QSP_ID_TUM_CX)) + QSP_CONST(QSP_VOL_TCELL) * (GET_VAR_RAW(QSP_ID_TUM_TEFF)+GET_VAR_RAW(QSP_ID_TUM_TREG)+GET_VAR_RAW(QSP_ID_TUM_TEXH))) / AVOGADROS;
	double tumor_volume_ref = (PI * std::pow(QSP_CONST(QSP_INIT_TUM_DIAM),3)) / 6;
	double cent_nivo = GET_VAR(QSP_ID_CENT_NIVO);

	//Nivo
	if (params.getVal(PARAM_NIVO_ON) != 0 && (tumor_volume > tumor_volume_ref || cent_nivo))
	{
		double week = t / (SEC_PER_DAY * params.getVal(PARAM_NIVO_DOSE_INTERVAL_TIME));
		int week_int = floor(week);
		double nivo_dose = params.getVal(PARAM_NIVO_DOSE);

		if (week == week_int)
		{
			cent_nivo += nivo_dose;
			SET_VAR(QSP_ID_CENT_NIVO, cent_nivo);
		}
	}

	//Durv
		if (params.getVal(PARAM_DURV_ON) != 0)
	{
		double week = t / (SEC_PER_DAY * params.getVal(PARAM_DURV_DOSE_INTERVAL_TIME));
		int week_int = floor(week);
		double durv_dose = params.getVal(PARAM_DURV_DOSE);
		double cent_durv = GET_VAR(QSP_ID_CENT_DURV);

		if (week == week_int)
		{
			cent_durv += durv_dose;
			SET_VAR(QSP_ID_CENT_DURV, cent_durv);
		}
	}

	//Ent
		if (params.getVal(PARAM_ENT_ON) != 0)
	{
		double week = t / (SEC_PER_DAY * params.getVal(PARAM_ENT_DOSE_INTERVAL_TIME));
		int week_int = floor(week);
		double ent_dose = params.getVal(PARAM_ENT_DOSE);
		double cent_ent = GET_VAR(QSP_ID_CENT_ENT);
		std::cout<<"ENT1: "<< week_int <<std::endl;
		std::cout<<"ENT2: "<< ent_dose <<std::endl;
		std::cout<<"ENT3: "<< cent_ent <<std::endl;

		if (week == week_int)
		{
			cent_ent += ent_dose;
			std::cout<<"ENT: "<< cent_ent <<std::endl;
			SET_VAR(QSP_ID_CENT_ENT, cent_ent);
		}
	}	

	//Ipi
		if (params.getVal(PARAM_IPI_ON) != 0)
	{
		double week = t / (SEC_PER_DAY * params.getVal(PARAM_IPI_DOSE_INTERVAL_TIME));
		int week_int = floor(week);
		double ipi_dose = params.getVal(PARAM_IPI_DOSE);
		double cent_ipi = GET_VAR(QSP_ID_CENT_IPI);

		if (week == week_int)
		{
			cent_ipi += ipi_dose;
			SET_VAR(QSP_ID_CENT_IPI, cent_ipi);
		}
	}	

	// solve QSP for dt
	_QSP_model.solve(t, dt);

	return;
}

/*! Get QSP variables for ABM.
	
	# Tum.C1 (unit: cell)
    # Cent.Teff (unit: convert from cell to mole)
	# Cent.Treg (unit: convert from cell to mole)
	# Tum.MDSC (unit: convert from cell to mole)
	# Tum.Nivo (unit: convert from 1e-6 mole/m^3 to mole/m^3)
	# Tum.ENT (unit: convert from 1e-6 mole/m^3 to mole/m^3)
    # Tum.Cx (unit: convert from cell to mole)
	# Tum.Texh (unit: convert from cell to mole)
	# Tum.CCL2 (unit: convert from 1e-6 mole/m^3 to mole/m^3)	
*/
const std::vector<double>& LymphCentral::get_var_exchange(void){

	// need to be cell count for calculating ABM scalor
	_var_qsp_to_abm[QSPEX_TUM_C] = GET_VAR_RAW(QSP_ID_TUM_C1);
	// internal SI unit for calculating rates and probabilities.
	_var_qsp_to_abm[QSPEX_CENT_TEFF] = GET_VAR(QSP_ID_CENT_TEFF);
	_var_qsp_to_abm[QSPEX_CENT_TREG] = GET_VAR(QSP_ID_CENT_TREG);
	_var_qsp_to_abm[QSPEX_TUM_TEFF] = GET_VAR(QSP_ID_TUM_TEFF);
	_var_qsp_to_abm[QSPEX_TUM_TREG] = GET_VAR(QSP_ID_TUM_TREG);
	_var_qsp_to_abm[QSPEX_TUM_MDSC] = GET_VAR(QSP_ID_TUM_MDSC);
	_var_qsp_to_abm[QSPEX_TUM_NIVO] = GET_VAR(QSP_ID_TUM_NIVO);
	_var_qsp_to_abm[QSPEX_TUM_CX] = GET_VAR(QSP_ID_TUM_CX);
	_var_qsp_to_abm[QSPEX_TUM_TEXH] = GET_VAR(QSP_ID_TUM_TEXH);		
	_var_qsp_to_abm[QSPEX_PD1_PDL1] = GET_VAR(QSP_ID_PD1_PDL1);
	_var_qsp_to_abm[QSPEX_TUM_CCL2] = GET_VAR(QSP_ID_TUM_CCL2);	
	_var_qsp_to_abm[QSPEX_TUM_ARGI] = GET_VAR(QSP_ID_TUM_ARGI);
	_var_qsp_to_abm[QSPEX_TUM_NO] = GET_VAR(QSP_ID_TUM_NO);
	_var_qsp_to_abm[QSPEX_TUM_VOL] = params.getVal(PARAM_VT_MIN) + params.getVal(PARAM_V_CELL) * (GET_VAR_RAW(QSP_ID_TUM_C1)/AVOGADROS + GET_VAR(QSP_ID_TUM_CX)) + params.getVal(PARAM_V_TCELL) * (GET_VAR(QSP_ID_TUM_TEFF) + GET_VAR(QSP_ID_TUM_TREG) + GET_VAR(QSP_ID_TUM_TEXH));
	return _var_qsp_to_abm;
}

/*! update QSP module with output from ABM
	unit convert: item to mole
    # cancer cell death (total)
	# cancer cell death (Teff kill)
	# Teff recruitment
	# Treg recruitment
*/
void LymphCentral::update_qsp_var(const std::vector<double>& var_abm){

	// convert item to internal units
	double scalar = 1 / AVOGADROS;

	// CC death total, CC death Teff, Teff recruit, Treg recruit
	double cc_death_total = var_abm[Tumor::TUMEX_CC_DEATH] * scalar;
	double cc_death_Teff = var_abm[Tumor::TUMEX_CC_T_KILL] * scalar;
	double Teff_recruit = var_abm[Tumor::TUMEX_TEFF_REC] * scalar;
	double Treg_recruit = var_abm[Tumor::TUMEX_TREG_REC] * scalar;
	double cc_molec = var_abm[Tumor::TUMEX_CC];
	double Teff_molec = var_abm[Tumor::TUMEX_T];
	double Treg_molec = var_abm[Tumor::TUMEX_TREG];
	double MDSC_molec = var_abm[Tumor::TUMEX_MDSC];
	double cx_molec = var_abm[Tumor::TUMEX_CX];
	double Texh_molec = var_abm[Tumor::TUMEX_TEXH];


	SET_VAR_RAW(QSP_ID_TUM_C1, cc_molec);
	SET_VAR(QSP_ID_TUM_C1, (cc_molec)*scalar);

	if (Teff_molec != 0){
	SET_VAR_RAW(QSP_ID_TUM_TEFF, Teff_molec);
	SET_VAR(QSP_ID_TUM_TEFF, (Teff_molec)*scalar);
	}

	if (Treg_molec != 0){
	SET_VAR_RAW(QSP_ID_TUM_TREG, Treg_molec);
	SET_VAR(QSP_ID_TUM_TREG, (Treg_molec)*scalar);
	}

	if (MDSC_molec != 0){
	SET_VAR_RAW(QSP_ID_TUM_MDSC, MDSC_molec);
	SET_VAR(QSP_ID_TUM_MDSC, (MDSC_molec)*scalar);
	}

	if (cx_molec != 0){
	SET_VAR_RAW(QSP_ID_TUM_CX, cx_molec);
	SET_VAR(QSP_ID_TUM_CX, (cx_molec)*scalar);	
	}

	if (Texh_molec != 0){
	SET_VAR_RAW(QSP_ID_TUM_TEXH, Texh_molec);
	SET_VAR(QSP_ID_TUM_TEXH, (Texh_molec)*scalar);	
	}

	return;
}

};
};
