#pragma once

#include "SP_QSP_shared/ABM_Base/SpatialCompartment.h"
#include "TumorGridVoxel.h"
#include "VoxelContentGen.h"
#include "../agent/TCell.h"
#include "../agent/CancerCell.h"
#include "../agent/Mac.h"
#include "../agent/Fib.h"
#include "../agent/TReg.h"
#include "../agent/MDSC.h"
#include "../../pde/DiffuseGrid.h"
#include "../../core/Stats.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
//#include <boost/serialization/export.hpp>

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

typedef std::vector<Coord3D> CellSource;
typedef boost::property_tree::ptree icProperty;

//! Generic tumor compartment
class Tumor : public SpatialCompartment
{
public:
enum TumExVar{
	TUMEX_CC,
	TUMEX_T,
	TUMEX_TEFF_EXH,
	TUMEX_TREG,
	TUMEX_MDSC,
	TUMEX_CX,
	TUMEX_TEXH,
	TUMEX_CC_DEATH,
	TUMEX_CC_T_KILL,
	TUMEX_TEFF_REC,
	TUMEX_TREG_REC,
	TUMEX_CC_D,
	TUMEX_TEFF_D,		
	TUMEX_TREG_D,
	TUMEX_MDSC_D,
	TUMEX_CX_D,
	TUMEX_TEXH_D,
	TUMEX_VAR_NUM
};

public:
	//! Need default constructor for serialization to work
	Tumor() {};
	Tumor(int x, int y, int z);
	virtual ~Tumor();

	// add one source of entry
	void add_lymphocyte_source(const Coord& c);

	// add one source of entry
	void add_mdsc_source(const Coord& c);	

	//! simulate for one slice
	virtual void timeSlice(unsigned long slice, int cc, int teff, int treg, int mdsc, int texh);

	//! recruit T Cells from invasive front 
	void recruitTCellsFront(double min, double max);

	//! get stats
	const Stats& get_stats(void) const{ return _stats; };

	//! default header for extra remark column when writing cell grid to file
	virtual std::string getExtraRemarkHeader() const;
	//! print grid snapshot to file
	std::string printGridToFile() const;
	//! print cell ODE stats to file
	void printCellOdeToFile(unsigned long slice) const;
	//! print grid snapshot to screen
	void printGridToScreen(unsigned long slice) const;
	//! get chemical grid
	DiffuseGrid & get_chem_grid(void) { return _chem; };
	//! get concentration of chemokine i.
	double get_chem(const Coord3D&c, chem_ID i)const; 

	//! add to abm var exchange counter
	void inc_abm_var_exchange(TumExVar v){ _var_abm_to_qsp[v] += 1.0;};
	//! return variables needed for QSP module 
	const std::vector<double>& get_var_exchange(void);

	//! update ABM module with variables from QSP 
	void update_abm_with_qsp(const std::vector<double>&);

	//! tumor volume
	double get_Tum_Vol(void)const{ return _tumor_volume; };
	//! teff concentration in the tumor
	double get_T(void)const{ return _concentration_t_eff_tum; };	
	//! treg concentration in the tumor
	double get_Treg(void)const{ return _concentration_t_reg_tum; };
	//! nivo concentration
	double get_Nivo(void)const{ return _concentration_nivo; };
	//! cc
	double get_cc(void)const{ return _concentration_cc; };	
	//! texh concentration in the tumor
	double get_Texh(void)const{ return _concentration_t_exh; };
	//! pd1pdl1
	double get_pd1pdl1(void)const{ return _concentration_pd1pdl1; };	
	//! argi
	double get_argi(void)const{ return _concentration_argi; };	
	//! no
	double get_no(void)const{ return _concentration_no; };	

	//! create one initial cell and configure to explicity initial condicitons
	void TEST_AddOneCell(AgentType type, AgentState state, const Coord3D& c) {
		createOneInitCell(type, state, c);
	};

	//! set allow shift
	void set_allow_shift(bool allowed){ _allow_shift_grid = allowed; };


	VoxelContentGen _voxel_ic;

protected:
	//! initialize agent grid 
	virtual bool initAgentGrid();	

private:

	friend class boost::serialization::access;
	template<class Archive>
	//! boost serialization
	void serialize(Archive & ar, const unsigned int /*version*/);

	void time_slice_recruitment(double t, int cc_p, int teff_p, int reg_p, int mdsc_p, int texh_p);
	void time_slice_movement(double t, double dt);
	void time_slice_state_change(double t, double dt, int cc_p, int teff_p, int treg_p, int mdsc_p, int texh_p);
	//! scan all agents for the last round in a slice
	void time_slice_final_scan(void);
	void time_slice_molecular(double t);

	//! initialize compartment: setup environmen��t
	void initEnvironment();
	//! initialize compartment: setup initial cells
	void initCell(std::string filename);
	//! initial cell: center cancer cell
	void init_cell_single_center(void);
	//! initial cell: fill grid 
	void init_cell_fill_grid(void);
	//! create one random cell
	bool populate_voxel_random(const Coord3D&);
	//! create a cluster of initial cells centered around provided coordinate
	void createClusterInitCell(icProperty &ic);
	//! create one initial cell and configure to explicity initial condicitons
	CellAgent* createOneInitCell(AgentType type, AgentState state, const Coord3D& crd);
	//! Teff/Treg recruitment probability at each entry point
	double get_T_recruitment_prob(double c, double base) const;
	//! MDSc recruitment probability at each entry point
	double get_MDSC_recruitment_prob(double c, double base) const;	
  
	//! adjust camera center by shifting grid
	void shift_adjust_center(void);
	//! get the vector of camera shift
	Coord3D get_cam_shift(void);
	//! shift contents of the grid by crd
	void shift_grid(Coord3D& crd);

	// stats
	Stats _stats;
	// Diffusion grid
	DiffuseGrid _chem;
	//! list of T cell sources
	CellSource _t_source;
	//! list of T cell sources
	CellSource _mdsc_source;
	//! Dummy T cell for initialization
	TCell * _tInitDummy;
	TCell * _tInitDummy0;	
	//! Dummy MDSC for initialization
	MDSC * _MDSCInitDummy;	
	MDSC * _MDSCInitDummy0;
	//! Dummy Cancer cell for initialization
	CancerCell * _cInitDummy;
	Mac * _macInitDummy;
	Fib * _fibInitDummy;
	Treg * _TregInitDummy;
	Treg * _TregInitDummy0;
	//! voxel dimension, no need to serialize
	//double _voxel_size;

	/* variables used in grid shifting */

	//! allow shifting
	bool _allow_shift_grid;
	//! mass center anchor
	Coord3D _center_target;
	//! temporary AgentGrid
	AgentGrid _agGrid_temp;

	/* following are variables communicated between tumor and blood */

	//! variable values passed from QSP to ABM (no need to serialize)
	std::vector<double> _var_abm_to_qsp;
	//! concentration of cancer cells(no need to serialize).
	double _concentration_cc;
	//! concentration of cytotoxic t cells in blood(Unit: SI; no need to serialize).
	double _concentration_t_cyt;
	//! concentration of regulatory t cells in blood(Unit: SI; no need to serialize).
	double _concentration_t_reg;
	//! concentration of cytotoxic t cells in the tumor(Unit: SI; no need to serialize).
	double _concentration_t_eff_tum;	
	//! concentration of regulatory t cells in the tumor(Unit: SI; no need to serialize).
	double _concentration_t_reg_tum;	
	//! concentration of mdsc(Unit: SI; no need to serialize).
	double _concentration_mdsc;
	//! concentration of nivo in Tumor(Unit: SI; no need to serialize).
	double _concentration_nivo;
	//! concentration of dead cancer cells in Tumor(Unit: SI; no need to serialize).
	double _concentration_cx;
	//! concentration of exhausted t cells in Tumor(Unit: SI; no need to serialize).
	double _concentration_t_exh;
	//! concentration of ent in Tumor(Unit: SI; no need to serialize).
	//double _concentration_ent;	
	//! pd1-pdl1(Unit: SI; no need to serialize).
	double _concentration_pd1pdl1;					
	//! concentration of ccl2 in Tumor(Unit: SI; no need to serialize).	
	double _concentration_ccl2;	
	//! concentration of argi in Tumor(Unit: SI; no need to serialize).	
	double _concentration_argi;			
	//! concentration of no in Tumor(Unit: SI; no need to serialize).	
	double _concentration_no;		
	//! tumor volume(Unit: SI; no need to serialize).
	double _tumor_volume;

};

//BOOST_CLASS_EXPORT_KEY(Tumor);

template<class Archive>
inline void Tumor::serialize(Archive & ar, const unsigned int /* version */) {
	ar.template register_type<CancerCell>();
	ar.template register_type<TCell>();
	ar.template register_type<MDSC>();
	ar.template register_type<Mac>();
	ar.template register_type<Fib>();
	ar.template register_type<Treg>();
	ar.template register_type<TumorGridVoxel>();
	//ar.template register_type<BioFVMSinkSource>();
	// this is needed because cell._compartment is of type SpatialCompartment*
	ar.template register_type<Tumor>();
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpatialCompartment);
	ar & BOOST_SERIALIZATION_NVP(_voxel_ic);
	ar & BOOST_SERIALIZATION_NVP(_stats);
	ar & BOOST_SERIALIZATION_NVP(_t_source);
	ar & BOOST_SERIALIZATION_NVP(_mdsc_source);
	ar & BOOST_SERIALIZATION_NVP(_chem);
	ar & BOOST_SERIALIZATION_NVP(_allow_shift_grid);
	ar & BOOST_SERIALIZATION_NVP(_center_target);
	/*
	ar & BOOST_SERIALIZATION_NVP(_agGrid_temp);
	ar & BOOST_SERIALIZATION_NVP(_tInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_cInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_macInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_fibInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_TregInitDummy);
	*/
}

};
};

