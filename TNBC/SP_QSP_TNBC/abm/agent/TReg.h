#pragma once

#include "Cell_Tumor.h"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

#include "../../pde/DiffuseGrid.h"

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

class Treg:
	public Cell_Tumor
{
public:
	Treg(){};
	Treg(SpatialCompartment* c );
	Treg(const Treg& c);
	virtual ~Treg();

	virtual CellAgent* createCellCopy() const { return new Treg(*this); };

	//! print cancer cell information
	virtual std::string toString() const;

	static int getTregLife();

	//! step function for cancer cell
	bool agent_movement_step(double t, double dt, Coord& c);
	bool agent_state_step(double t, double dt, Coord& c, int cc_p, int teff_p, int treg_p, int mdsc_p, int texh_p);

	//! move sources (IL10)
	void move_all_source_sink(void)const;

	//! get cell agent type, in this case cancer cell.
	virtual AgentType getType() const { return AgentTypeEnum::CELL_TYPE_TREG; };

private:
	friend class boost::serialization::access;
	//! boost serialization
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);

	//! cool down timer for division
	int _divide_cd_Treg_exp;
	//! number of remaining division
	int _divide_limit_Treg_exp;
	
	//BioFVMSinkSource* _source_IL_10;
};

//BOOST_CLASS_EXPORT_KEY(Treg)

template<class Archive>
inline void Treg::serialize(Archive & ar, const unsigned int /* version */){
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Cell_Tumor);
	ar & BOOST_SERIALIZATION_NVP(_divide_cd_Treg_exp);
	ar & BOOST_SERIALIZATION_NVP(_divide_limit_Treg_exp);	
	//ar & BOOST_SERIALIZATION_NVP(_source_IL_10);
}

};
};
