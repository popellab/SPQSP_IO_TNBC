#pragma once

#include "SP_QSP_shared/ABM_Base/StatsBase.h"
#include "AgentEnum.h"

#include <boost/serialization/base_object.hpp>

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

enum StatsEventShared
{
	STATS_EVENT_REC,
	STATS_EVENT_PROLIF,
	STATS_EVENT_DIE,
	STATS_EVENT_MOVE,
	STATS_EVENT_DROP_IN,
	STATS_EVENT_DROP_OUT,
	STATS_EVENT_SHARED_COUNT,
};

enum StatsEventSpecial
{
	STATS_EVENT_SPECIAL_COUNT,
};

enum StatsMisc {
	STATS_MISC_PDL1_POS,
	STATS_MISC_COUNT,
};

//! record and output statistics 
class Stats : public StatsBase
{
public:
	Stats();
	~Stats();

	//! increment recruitment count
	void incRecruit(AgentType, AgentState);
	//! increment proliferation count
	void incProlif(AgentType, AgentState);
	//! increment death count
	void incDeath(AgentType, AgentState);
	//! increment movement count
	void incMove(AgentType, AgentState);
	//! increment cell generated otherwise 
	void incDropIn(AgentType, AgentState);
	//! increment cell dropping out of grid 
	void incDropOut(AgentType, AgentState);

	//! get total number of eff/cyt T cells
	int getTCell() const;
	//! get total number of cancer cells
	int getCancerCell() const;
	//! get total number of Treg cell
	int getTreg() const;
	//! get total number of MDSC cell
	int getMDSC() const;
	//! get total number of exh T cells
	int getTexh() const;		

private:
	
	friend class boost::serialization::access;

	virtual void initStats();

	// interval stats accumulator
	template<class Archive>
	//! boost serialization 
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StatsBase);
	}

};

};
};
