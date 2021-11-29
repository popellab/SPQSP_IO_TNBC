#include "Stats.h"

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

using namespace std;

struct cellTypeState{
	const char* name;
	AgentTypeEnum type;
	AgentStateEnum state;
};

const cellTypeState keygenCellTypeState[] = 
{
	{ "CD8.effector", AgentTypeEnum::CELL_TYPE_T, AgentStateEnum::T_CELL_EFF },
	{ "CD8.cytotoxic", AgentTypeEnum::CELL_TYPE_T, AgentStateEnum::T_CELL_CYT },
	{ "CD8.suppressed", AgentTypeEnum::CELL_TYPE_T, AgentStateEnum::T_CELL_SUPP },
	{ "Treg.default", AgentTypeEnum::CELL_TYPE_TREG, AgentStateEnum::DEFAULT_STATE},
	{ "MDSC.default", AgentTypeEnum::CELL_TYPE_MDSC, AgentStateEnum::DEFAULT_STATE},	
	{ "cancerCell.Stem", AgentTypeEnum::CELL_TYPE_CANCER, AgentStateEnum::CANCER_STEM},
	{ "cancerCell.Progenitor", AgentTypeEnum::CELL_TYPE_CANCER, AgentStateEnum::CANCER_PROGENITOR},
	{ "cancerCell.Senescent", AgentTypeEnum::CELL_TYPE_CANCER, AgentStateEnum::CANCER_SENESCENT},
	{ "mac.default", AgentTypeEnum::CELL_TYPE_MAC, AgentStateEnum::DEFAULT_STATE},
	{ "fib.default", AgentTypeEnum::CELL_TYPE_FIB, AgentStateEnum::DEFAULT_STATE}
};

const char* eventSharedName[] = {
	"recruit",
	"prolif",
	"death",
	"move",
	"drop_in",
	"drop_out"
};

const char* eventSpecialName[] = {
	 "dummy" 
};

const char* miscStatsName[] = {
	 "PDL1_frac" 
};


Stats::Stats()
: StatsBase()
{
	initStats();
}

Stats::~Stats()
{
}

void Stats::incRecruit(AgentType e, AgentState s){
	incEventShared(STATS_EVENT_REC, e, s);
}

void Stats::incProlif(AgentType e, AgentState s){
	incEventShared(STATS_EVENT_PROLIF, e, s);
}

void Stats::incDeath(AgentType e, AgentState s){
	incEventShared(STATS_EVENT_DIE, e, s);
}

void Stats::incMove(AgentType e, AgentState s){
	incEventShared(STATS_EVENT_MOVE, e, s);
}

//! increment cell generated otherwise 
void Stats::incDropIn(AgentType e, AgentState s){
	incEventShared(STATS_EVENT_DROP_IN, e, s);
}
//! increment cell dropping out of grid 
void Stats::incDropOut(AgentType e, AgentState s){
	incEventShared(STATS_EVENT_DROP_OUT, e, s);
}

int Stats::getTCell() const {
	int n = 0;
	n += getCountTypeState(AgentTypeEnum::CELL_TYPE_T, AgentStateEnum::T_CELL_EFF);
	n += getCountTypeState(AgentTypeEnum::CELL_TYPE_T, AgentStateEnum::T_CELL_CYT);
	return n;
}
int Stats::getTexh() const {
	int n = 0;
	n += getCountTypeState(AgentTypeEnum::CELL_TYPE_T, AgentStateEnum::T_CELL_SUPP);
	return n;
}
int Stats::getCancerCell() const {
	int n = 0;
	n += getCountTypeState(AgentTypeEnum::CELL_TYPE_CANCER, AgentStateEnum::CANCER_STEM);
	n += getCountTypeState(AgentTypeEnum::CELL_TYPE_CANCER, AgentStateEnum::CANCER_PROGENITOR);
	n += getCountTypeState(AgentTypeEnum::CELL_TYPE_CANCER, AgentStateEnum::CANCER_SENESCENT);
	return n;
}

int Stats::getTreg()const{
	return getCountTypeState(AgentTypeEnum::CELL_TYPE_TREG, AgentStateEnum::DEFAULT_STATE);
}

int Stats::getMDSC()const{
	return getCountTypeState(AgentTypeEnum::CELL_TYPE_MDSC, AgentStateEnum::DEFAULT_STATE);
}


/*! Initialize stats member variables, including:
	# Headers:
		# _typeStateHeader
		# _eventHeader
		# _eventHeaderSpecial
	# Counters:
		# _agCount
		# eventTypeState
		# eventSpecial
*/
void Stats::initStats(){

	size_t tsSize = sizeof(keygenCellTypeState) / sizeof(keygenCellTypeState[0]);
	// headers for agent type/stats; _tsMap
	for (size_t i = 0; i < tsSize; i++)
	{
		_typeStateHeader.push_back(keygenCellTypeState[i].name);
		setTypeStateIndex(keygenCellTypeState[i].type, keygenCellTypeState[i].state, i);
	}
	_agCount = std::vector<int>(tsSize, 0);

	// headers for shared events
	for (size_t i = 0; i < STATS_EVENT_SHARED_COUNT; i++)
	{
		_eventHeader.push_back(eventSharedName[i]);
	}
	_eventTypeState = std::vector<std::vector<long>>(STATS_EVENT_SHARED_COUNT, 
		std::vector<long>(tsSize, 0));

	// headers for special events
	for (size_t i = 0; i < STATS_EVENT_SPECIAL_COUNT; i++)
	{
		_eventHeaderSpecial.push_back(eventSpecialName[i]);
	}
	_eventSpecial = std::vector<long>(STATS_EVENT_SPECIAL_COUNT, 0);

	for (size_t i = 0; i < STATS_MISC_COUNT; i++)
	{
		_miscHeader.push_back(miscStatsName[i]);
	}
	_misc = std::vector<double>(STATS_MISC_COUNT, 0);
}

};
};