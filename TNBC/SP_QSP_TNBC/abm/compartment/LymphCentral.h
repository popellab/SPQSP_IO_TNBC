#pragma once

#include <boost/serialization/nvp.hpp>

#include <string>
#include <vector>

#include "SP_QSP_shared/Numerical_Adaptor/CVODE/MolecularModelCVode.h"
#include "TNBC/SP_QSP_TNBC/ode/ODE_system.h"


namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

class LymphCentral
{
typedef CancerVCT::ODE_system LymphBloodQSP;
typedef CancerVCT::Param LymphBloodParam;
typedef MolecularModelCVode<LymphBloodQSP> QSP;

public:
enum QSPExVar{
	QSPEX_TUM_C,
	QSPEX_CENT_TEFF,
	QSPEX_CENT_TREG,
	QSPEX_TUM_TEFF,
	QSPEX_TUM_TREG,
	QSPEX_TUM_MDSC,
	QSPEX_TUM_NIVO,
	QSPEX_TUM_CX,
	QSPEX_TUM_TEXH,
	QSPEX_PD1_PDL1,
	QSPEX_TUM_CCL2,
	QSPEX_TUM_ARGI,
	QSPEX_TUM_NO,	
	QSPEX_TUM_VOL,
	QSPEX_VAR_NUM
};
public:
	LymphCentral();
	~LymphCentral();

	//! setup parameters of ODE (one time)
	void setup_param(LymphBloodParam& p);

	//! time step 
	void time_step(double t, double dt);
	//! return variable values from QSP module that are needed for ABM
	const std::vector<double>& get_var_exchange(void);
	//! update QSP module with output from ABM
	void update_qsp_var(const std::vector<double>&);
	//! QSP headers
	std::string getQSPHeaders(void)const { return LymphBloodQSP::getHeader();};
	//! write QSP variables 
	friend std::ostream & operator<<(std::ostream &os, const LymphCentral& l);

private:

	//! boost serialization
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);


	//! QSP model excluding tumor dynamics.
	QSP _QSP_model;

	//! variable values passed from QSP to ABM
	std::vector<double> _var_qsp_to_abm;
};

template<class Archive>
inline void LymphCentral::serialize(Archive & ar, const unsigned int  version) {
	//ar & BOOST_SERIALIZATION_NVP(_cancer_debris);// no need to serialize
	ar & BOOST_SERIALIZATION_NVP(_QSP_model);
	LymphBloodQSP::classSerialize(ar, version);
}

inline std::ostream & operator<<(std::ostream &os, const LymphCentral& l){
	os << l._QSP_model;
	return os;
}

};
};

