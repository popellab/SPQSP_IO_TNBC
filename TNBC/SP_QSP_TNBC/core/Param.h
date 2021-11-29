#pragma once

#include "SP_QSP_shared/ABM_Base/ParamBase.h"

#include <string>

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

//! enumerator for double type parameters 
enum ParamFloat{
	// QSP
	PARAM_WEIGHT_QSP,
	PARAM_QSP_STEADYSTATE,
	PARAM_QSP_T_RESECTION,
	// ENV
	PARAM_SEC_PER_TIME_SLICE,
	PARAM_CELLS_SCALING_FACTOR,
	PARAM_SOURCES,
	PARAM_REC_SITE_FACTOR,
	PARAM_ADH_SITE_DENSITY,
	// PHARM
	PARAM_NIVO_DOSE_INTERVAL_TIME,
	PARAM_NIVO_DOSE,
	PARAM_DURV_DOSE_INTERVAL_TIME,
	PARAM_DURV_DOSE,	
	PARAM_ENT_DOSE_INTERVAL_TIME,
	PARAM_ENT_DOSE,	
	PARAM_IPI_DOSE_INTERVAL_TIME,
	PARAM_IPI_DOSE,	
	// T cell
	PARAM_T_CELL_LIFE_MEAN,
	PARAM_T_CELL_LIFE_SD,
	PARAM_T_CELL_MOVE_PROB,
	PARAM_IL_2_RELEASE_TIME,
	PARAM_IL_2_PROLIF_TH,
	PARAM_IFN_RELEASE_TIME,
	// Treg
	PARAM_TREG_MOVE_PROB,
	// MDSC
	PARAM_MDSC_MOVE_PROB,
	// Cancer cell
	PARAM_CANCER_STEM_ASYMMETRIC_DIV_PROB,
	PARAM_CANCER_SENESCENT_DEATH_RATE,
	PARAM_CANCER_STEM_MOVE_PROB,
	PARAM_CANCER_CELL_MOVE_PROB,
	PARAM_C1_MIN,
	PARAM_IFN_G_UPTAKE,
	// Agent cytokine 
	PARAM_PDL1_HIGH_TH,
	PARAM_IFN_G_PDL1_HALF,
	PARAM_IFN_G_PDL1_N,
	PARAM_PDL1_DECAY_DAY,
	// diffusion grid
	PARAM_IFN_G_DIFFUSIVITY,
	PARAM_IFN_G_RELEASE,
	PARAM_IFN_G_DECAY_RATE,
	PARAM_IL_2_DIFFUSIVITY,
	PARAM_IL_2_RELEASE,
	PARAM_IL_2_DECAY_RATE,
	PARAM_CCL2_DIFFUSIVITY,
	PARAM_CCL2_RELEASE,
	PARAM_CCL2_DECAY_RATE,	
	PARAM_CCL2_MOLECULAR_WEIGHT,
	PARAM_ARGI_DIFFUSIVITY,
	PARAM_ARGI_RELEASE,
	PARAM_ARGI_DECAY_RATE,
	PARAM_ARGI_MOLECULAR_WEIGHT,
	PARAM_NO_DIFFUSIVITY,
	PARAM_NO_RELEASE,
	PARAM_NO_DECAY_RATE,
	PARAM_NO_MOLECULAR_WEIGHT,
	// dummy
	PARAM_FLOAT_COUNT // dummy for count
};

//! enumerator for int type parameters 
enum ParamInt{
	PARAM_TUMOR_X,
	PARAM_TUMOR_Y,
	PARAM_TUMOR_Z,
	PARAM_VOXEL_SIZE,
	PARAM_N_T_VOXEL,
	PARAM_N_T_VOXEL_C,
	PARAM_SHUFFLE_CELL_VEC_INTERVAL,
	PARAM_SHIFTGRID_INTERVAL,
	PARAM_T_DIV_INTERVAL,
	PARAM_T_DIV_LIMIT,	
	PARAM_TREG_EXP_DIV_LIMIT,
	PARAM_CANCER_CELL_PROGENITOR_DIV_MAX,
	// lymphatic compartment
	// molecular
	PARAM_MOLECULAR_STEP_PER_SLICE,
	// dummy
	PARAM_INT_COUNT // dummy for count
};

//! enumerator for boolean type parameters 
enum ParamBool{
	PARAM_QSP_RESECTION,
	PARAM_ALL_MOLECULAR_OFF,
	PARAM_DIFFUSION_OFF,
	PARAM_NIVO_ON,
	PARAM_DURV_ON,
	PARAM_ENT_ON,
	PARAM_IPI_ON,	
	PARAM_BOOL_COUNT
};

//! parameters not directly from paramter file (float type, place holder)
enum ParamFloatInternal{
	PARAM_VOXEL_SIZE_CM,
	PARAM_CANCER_SENESCENT_MEAN_LIFE,
	PARAM_T_CELL_LIFE_MEAN_SLICE,
	PARAM_T_CELL_LIFE_SD_SLICE,
	PARAM_PDL1_DECAY_SLICE,
	// parameters calculated from QSP param
	PARAM_AVOGADROS,
	PARAM_CELL,	
	PARAM_VT_MIN,
	PARAM_V_CELL,
	PARAM_V_TCELL,
	PARAM_CMAX,
	PARAM_TREGMAX,
	PARAM_MDSCMAX,
	PARAM_K_BASE_REC,
	PARAM_K_REC_BY_CCL2,	
	PARAM_EC50_CCL2_REC,
	PARAM_IC50_NO_CTL,
	PARAM_IC50_ARGI_CTL,
	PARAM_EC50_ARGI_TREG,	
	PARAM_IC50_ENT_C,	
	PARAM_IC50_ENT_CCL2,	
	PARAM_IC50_ENT_NO,
	PARAM_IC50_ENT_ARGI,
	PARAM_PD1_PDL1_HALF,
	PARAM_PD1_SYN,
	PARAM_PDL1_SYN_MAX,
	PARAM_PDL1_K1,
	PARAM_PDL1_K2,
	PARAM_PDL1_K3,
	PARAM_N_PD1_PDL1,
	PARAM_ESCAPE_BASE,
	PARAM_EXHUAST_BASE_PDL1,
	PARAM_EXHUAST_BASE_TREG,
	PARAM_RESECT_TIME_STEP,
	PARAM_REC_PORT_PROB,
	PARAM_TEFF_RECRUIT_K,
	PARAM_TREG_RECRUIT_K,
	PARAM_MDSC_RECRUIT_K,
	PARAM_MDSC_RECRUIT_BY_CCL2_K,
	PARAM_TUM_MAX_C,
	PARAM_TREG_LIFE_MEAN,
	PARAM_MDSC_LIFE_MEAN,
	PARAM_TREG_EXP_INTERVAL_SLICE,
	PARAM_CSC_GROWTH_RATE,
	PARAM_CANCER_PROG_GROWTH_RATE,	
	PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE,
	PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE,
	PARAM_FLOAT_INTERNAL_COUNT
};

//! parameters not directly from paramter file (int type, place holder)
enum ParamIntInternal{
	PARAM_INT_INTERNAL_COUNT
};

//! parameters not directly from paramter file (boolean type)
enum ParamBoolInternal{
	PARAM_MOLECULAR_MODULES_ON,
	PARAM_DIFFUSION_ON,
	PARAM_IFN_SINK_ON,
	PARAM_BOOL_INTERNAL_COUNT
};
//! Model prameters
class Param: public ParamBase
{
public:
	Param();
	~Param(){};
	//! get parameter value (float)
	inline double getVal(ParamFloat n) const { return _paramFloat[n];};
	//! get parameter value (int)
	inline int getVal(ParamInt n) const { return _paramInt[n]; };
	//! get parameter value (bool)
	inline bool getVal(ParamBool n) const { return _paramBool[n]; };

	inline double getVal(ParamFloatInternal n) const { return _paramFloatInternal[n];};
	inline int getVal(ParamIntInternal n) const { return _paramIntInternal[n]; };
	inline bool getVal(ParamBoolInternal n) const { return _paramBoolInternal[n]; };

	//! update from QSP parameters
	void update_from_qsp(void);

private:

	//! setup content of _paramDesc
	virtual void setupParam();
	//! process all internal parameters
	virtual void processInternalParams();
};

};
};

