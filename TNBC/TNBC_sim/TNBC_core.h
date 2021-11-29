#pragma once

#include "TNBC/SP_QSP_TNBC/abm/compartment/Tumor.h"
#include "TNBC/SP_QSP_TNBC/abm/compartment/LymphCentral.h"

#include "FileOutputHub.h"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
//! Simulation central control
/*! Host all compartments (currently only one: tumor).
	Handle interaction between compartments; process global statistics output.
*/
namespace TNBC = SP_QSP_IO::SP_QSP_TNBC;
class TNBC_Core
{
public:

public:
	TNBC_Core();
	virtual ~TNBC_Core();
	//! setup qsp module
	void setup_qsp(CancerVCT::Param& p);
	//! Initialize all compartments: random population 
	void initializeSimulation(void);
	//! Proceed on slice in time for all compartments
	virtual void timeSlice(const long slice);
	//! Write ABM stats header to stats file
	void write_stats_header(void) const;
	//! Write ABM stats of current time slice to stats file
	void write_stats_slice(unsigned long slice) const;
	//! write QSP content to file
	void write_QSP(unsigned long slice, bool header)const;
	//! Write grid snapshots of current time slice to file 
	virtual void writeGrids(unsigned long slice, unsigned int option);
	//! Print brief summary of current slice
	virtual void briefStats(unsigned long slice);

		int cc_margin;
		int cc_core; 
		int teff_margin;
		int teff_core;		
		int treg_margin;
		int treg_core;
		int mdsc_margin;
		int mdsc_core;
		int cx_margin;
		int cx_core;			
		int texh_margin;
		int texh_core;

private:
	friend class boost::serialization::access;
	//! boost serialization 
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);

	//! tumor compartment instance
	std::vector<TNBC::Tumor*> _ROI_core;
	std::vector<TNBC::Tumor*> _ROI_margin;
	SP_QSP_IO::SP_QSP_TNBC::LymphCentral _lymph;
	SP_QSP_IO::SP_QSP_TNBC::Tumor _tumor;

};

template<class Archive>
inline void TNBC_Core::serialize(Archive & ar, const unsigned int version){
	ar & BOOST_SERIALIZATION_NVP(_ROI_core);
	ar & BOOST_SERIALIZATION_NVP(_ROI_margin);
	ar & BOOST_SERIALIZATION_NVP(_lymph);
	ar & BOOST_SERIALIZATION_NVP(_tumor);
	SP_QSP_IO::CellAgent::classSerialize(ar, version);
}