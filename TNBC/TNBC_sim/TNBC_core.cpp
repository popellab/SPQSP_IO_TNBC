#include "TNBC_core.h"

#include "TNBC/SP_QSP_TNBC/core/GlobalUtilities.h"
#include "InitialCondition.h"

#include <algorithm>    // std::max

#define AVOGADROS 6.022140857E23 
#define PI 3.1416

extern FileOutputHub output_hub;

extern RNG rng;

//extern SP_QSP_IO::Param params;

namespace SP_QSP_IO {
	namespace SP_QSP_TNBC {
		extern Param params;
	}
};
static auto& params = TNBC::params;

extern InitialCondition ic;
extern std::string initialCellFileName_core;
extern std::string initialCellFileName_margin;

typedef SP_QSP_IO::Coord Coord;

TNBC_Core::TNBC_Core()
: _ROI_core()
, _ROI_margin()
, _lymph()
, _tumor()
{
}

TNBC_Core::~TNBC_Core()
{
	for (auto& ptumor : _ROI_core) {
		delete ptumor;
	}
}

/*! Setup QSP module.
*/
void TNBC_Core::setup_qsp(CancerVCT::Param& p){
	_lymph.setup_param(p);
	params.update_from_qsp();
}
/*! initialize compartments: randomly populate voxels
	This function is called only when creating the model for the first time,
	and not reading from saved state. Objects created in this function should
	be already in the serialization and can be loaded directly.
*/
void TNBC_Core::initializeSimulation(void){

	for (int i = 0; i < ic.getVal(IC_NUM_ROI_core); i++)
	{
		auto pTumor = new TNBC::Tumor(params.getVal(TNBC::PARAM_TUMOR_X),
			params.getVal(TNBC::PARAM_TUMOR_Y),
			params.getVal(TNBC::PARAM_TUMOR_Z));
		_ROI_core.push_back(pTumor);
	}

	for (int i = 0; i < ic.getVal(IC_NUM_ROI_margin); i++)
	{
		auto pTumor = new TNBC::Tumor(params.getVal(TNBC::PARAM_TUMOR_X),
			params.getVal(TNBC::PARAM_TUMOR_Y),
			params.getVal(TNBC::PARAM_TUMOR_Z));
		_ROI_margin.push_back(pTumor);
	}

	std::string s;
	// rule to randomly populate voxel during initlaization or grid shifting
	for (auto& ptumor : _ROI_core) {
		auto& tumor = *ptumor;
		tumor.set_allow_shift(ic.getVal(IC_CORE_GRID_SHIFT));
		tumor._voxel_ic.setup(ic.getVal(IC_CORE_STATIONARY),
			ic.getVal(IC_DENSITY_CORE_CANCER),
			params.getVal(TNBC::PARAM_TUMOR_X),
			params.getVal(TNBC::PARAM_TUMOR_Y),
			params.getVal(TNBC::PARAM_TUMOR_Z));

		tumor.initCompartment(s);
		//std::cout << "core populated" << std::endl;
	}

	for (auto& ptumor : _ROI_margin) {

		auto& tumor = *ptumor;
		tumor.set_allow_shift(ic.getVal(IC_MARGIN_GRID_SHIFT));
		tumor._voxel_ic.setup(ic.getVal(IC_MARGIN_STATIONARY),
			ic.getVal(IC_DENSITY_MARGIN_CANCER),
			params.getVal(TNBC::PARAM_TUMOR_X),
			params.getVal(TNBC::PARAM_TUMOR_Y),
			params.getVal(TNBC::PARAM_TUMOR_Z));

		tumor.initCompartment(s);
		//std::cout << "margin populated" << std::endl;
	}

	//t cell sources
	// core:
	for (auto& ptumor : _ROI_core) {
		auto& tumor = *ptumor;
		std::vector<Coord> c_tumor;
		//std::cout << "num_source: " << c_tumor.size() << std::endl;
		unsigned int nr_source_tumor;
		tumor.for_each_grid_coord(true, true, true, [&](Coord&c){
			c_tumor.push_back(c);
		});
		//std::cout << "num_source: " << c_tumor.size() << std::endl;
		nr_source_tumor = int(params.getVal(TNBC::PARAM_SOURCES)*ic.getVal(IC_CORE_TUMOR_VAS_FOLD));
		rng.shuffle_first_k(c_tumor, nr_source_tumor);
		for (size_t i = 0; i < nr_source_tumor; i++)
		{
			tumor.add_lymphocyte_source(c_tumor[i]);
			tumor.add_mdsc_source(c_tumor[i]);
		}
		std::cout << "core nr sources: tumor: " << nr_source_tumor << std::endl;
	}

	// margin:
	for (auto& ptumor : _ROI_margin) {
		auto& tumor = *ptumor;
		std::vector<Coord> c_normal, c_tumor;
		unsigned int nr_source_normal, nr_source_tumor;
		tumor.for_each_grid_coord(true, true, true, [&](Coord&c){
			if (c.x >= 0 && c.x < params.getVal(TNBC::PARAM_TUMOR_X)
			&& c.y >= 0 && c.y < params.getVal(TNBC::PARAM_TUMOR_Y)
			&& c.z >= 0 && c.z < params.getVal(TNBC::PARAM_TUMOR_Z)){
				c_tumor.push_back(c);
			}
			else{
				c_normal.push_back(c);				
			}
		});
		nr_source_normal = int(0*params.getVal(TNBC::PARAM_SOURCES)*ic.getVal(IC_MARGIN_NORMAL_VAS_FOLD));
		rng.shuffle_first_k(c_normal, nr_source_normal);
		nr_source_tumor = int(1*params.getVal(TNBC::PARAM_SOURCES)*ic.getVal(IC_MARGIN_TUMOR_VAS_FOLD));
		rng.shuffle_first_k(c_tumor, nr_source_tumor);

		for (size_t i = 0; i < nr_source_normal; i++)
		{
			tumor.add_lymphocyte_source(c_normal[i]);
			tumor.add_mdsc_source(c_normal[i]);
		}
		for (size_t i = 0; i < nr_source_tumor; i++)
		{
			tumor.add_lymphocyte_source(c_tumor[i]);
			tumor.add_mdsc_source(c_tumor[i]);
		}
		std::cout << "margin nr sources: tumor: " << nr_source_tumor 
				<< ", normal: " << nr_source_normal << std::endl;
	}
	//std::cout << "sources generated" << std::endl;
}


void TNBC_Core::timeSlice(const long slice){

	const double dt = params.getVal(TNBC::PARAM_SEC_PER_TIME_SLICE);
	const double t0 = slice * dt;
	// std::cout << "RNG check (" << slice << ") START : " << rng.get_unif_01() << std::endl;

	/* update cancer number and blood concentration */
	auto& qsp_var = _lymph.get_var_exchange();
	double lymphCC = qsp_var[TNBC::LymphCentral::QSPEX_TUM_C];

	/* if QSP halted, skip*/
	std::cout << "lymph CC: " << lymphCC << std::endl;
	double abm_min_cc = params.getVal(TNBC::PARAM_C1_MIN);

if (lymphCC > abm_min_cc){
	for (auto& ptumor : _ROI_core) {		
		auto& abm_var_tum = ptumor->get_var_exchange();
		size_t abm_var_len = abm_var_tum.size();

		for (auto& ptumor : _ROI_margin) {
			ptumor->update_abm_with_qsp(qsp_var);
		}
		for (auto& ptumor : _ROI_core) {
			ptumor->update_abm_with_qsp(qsp_var);
		}

		std::cout << "nivo: " << qsp_var[6] << std::endl;

		/*
		for (auto& v : qsp_var)
		{
		std::cout << v << ", ";
		}
		std::cout << std::endl;
		*/	

		/* ABM time step */	

		for (auto& ptumor : _ROI_margin) {
			const auto& stats = _ROI_margin[0]->get_stats();
			cc_margin = stats.getCancerCell();
			teff_margin = stats.getTCell();
			treg_margin = stats.getTreg();
			mdsc_margin = stats.getMDSC();
			texh_margin = stats.getTexh();
			int ccm = cc_margin;
			int tem = teff_margin;
			int trm = treg_margin;
			int mdscm = mdsc_margin;
			int texhm = mdsc_margin;
			ptumor->timeSlice(slice, ccm, tem, trm, mdscm, texhm);				
		}
			for (auto& ptumor : _ROI_core) {
				const auto& stats = _ROI_core[0]->get_stats();
				cc_core = stats.getCancerCell();
				teff_core = stats.getTCell();
				treg_core = stats.getTreg();
				mdsc_core = stats.getMDSC();
				texh_core = stats.getTexh();
				int ccc = cc_core;
				int tec = teff_core;
				int trc = treg_core;
				int mdscc = mdsc_core;
				int texhc = texh_core;
				ptumor->timeSlice(slice, ccc, tec, trc, mdscc, texhc);
			}


		//std::cout << "RNG check (" << slice << ") CORE: " << rng.get_unif_01() << std::endl;
		//std::cout << "RNG check (" << slice << ") MARGI : " << rng.get_unif_01() << std::endl;

		/* update QSP variables */
		auto abm_var = std::vector<double>(abm_var_len, 0);
		auto abm_var_core = std::vector<double>(abm_var_len, 0);
		auto abm_var_margin = std::vector<double>(abm_var_len, 0);;
		for (auto& ptumor : _ROI_core) {
			auto& _abm_var_core = ptumor->get_var_exchange();
			for (size_t i = 0; i < abm_var_len; i++)
			{
				abm_var_core[i] += _abm_var_core[i]; 
			}
		}
		for (auto& ptumor : _ROI_margin) {
			auto& _abm_var_margin = ptumor->get_var_exchange();
			for (size_t i = 0; i < abm_var_len; i++)
			{
				abm_var_margin[i] += _abm_var_margin[i]; 
			}
		}
		double w = params.getVal(TNBC::PARAM_WEIGHT_QSP);

		double tumCC_core = abm_var_core[TNBC::Tumor::TUMEX_CC];
		double tumCC_margin = abm_var_margin[TNBC::Tumor::TUMEX_CC];
		double mm1;

		//Entire and partial tumor
		// example: SC =  50000
		//Entire tumor (tumor growth)
		// SC = 1

		double SC = params.getVal(TNBC::PARAM_CELLS_SCALING_FACTOR);	

		std::cout << "SC:" << SC << std::endl;	

		for (size_t i = 0; i < abm_var_len; i++){

		if (i==TNBC::Tumor::TUMEX_CC)
		{
			if (cc_core == 0){
				abm_var[i] = 0;
			}
			else{
				//double mm1 =  (abm_var_core[i]-cc_core-abm_var_core[TNBC::Tumor::TUMEX_CC_D])/cc_core;
			//abm_var[i] = qsp_var[TNBC::LymphCentral::QSPEX_TUM_C] * (1 + mm1);
			abm_var[i] =  SC*abm_var_core[i];
			//std::cout << "DROPOUT:" << abm_var_core[TNBC::Tumor::TUMEX_CC_D] << std::endl;	
			}
		}

		else if (i==TNBC::Tumor::TUMEX_T){	
			if (teff_core == 0){
				abm_var[i] = 0;
			}			
			else{
				double mm1 =  (abm_var_core[i]-teff_core)/teff_core - abm_var_core[TNBC::Tumor::TUMEX_TEFF_D]/teff_core;
			//abm_var[i] = qsp_var[TNBC::LymphCentral::QSPEX_TUM_TEFF]*AVOGADROS * (1 + mm1);	
			abm_var[i] = abm_var_core[i] * abm_var[TNBC::Tumor::TUMEX_CC]/abm_var_core[TNBC::Tumor::TUMEX_CC];
			std::cout << "TEFF_QSP:" << abm_var[i] << std::endl;
		}
		}
		else if (i==TNBC::Tumor::TUMEX_TREG){			
			if (treg_core == 0){
				abm_var[i] = 0;
			}			
			else{
				double mm1 =  (abm_var_core[i]-treg_core)/treg_core - abm_var_core[TNBC::Tumor::TUMEX_TREG_D]/treg_core;
			//abm_var[i] = qsp_var[TNBC::LymphCentral::QSPEX_TUM_TREG]*AVOGADROS * (1 + mm1);
			abm_var[i] = abm_var_core[i] * abm_var[TNBC::Tumor::TUMEX_CC]/abm_var_core[TNBC::Tumor::TUMEX_CC];
			std::cout << "TREG_QSP:" << abm_var[i] << std::endl;			
			}
		}
		else if (i==TNBC::Tumor::TUMEX_MDSC){	
			if (mdsc_core == 0){
				abm_var[i] = 0;
			}			
			else{
				double mm1 =  (abm_var_core[i]-mdsc_core)/mdsc_core - abm_var_core[TNBC::Tumor::TUMEX_MDSC_D]/mdsc_core;
			//abm_var[i] = qsp_var[TNBC::LymphCentral::QSPEX_TUM_MDSC]*AVOGADROS * (1 + mm1);
			abm_var[i] = abm_var_core[i] * abm_var[TNBC::Tumor::TUMEX_CC]/abm_var_core[TNBC::Tumor::TUMEX_CC];
			std::cout << "MDSC_QSP:" << abm_var[i] << std::endl;			
			}
		}
		else if (i==TNBC::Tumor::TUMEX_CX){	
			if (cx_core == 0){
				abm_var[i] = 0;
			}			
			else{
				double mm1 =  (abm_var_core[i]-cx_core)/cx_core - abm_var_core[TNBC::Tumor::TUMEX_CX_D]/cx_core;
			//abm_var[i] = qsp_var[TNBC::LymphCentral::QSPEX_TUM_CX]*AVOGADROS * (1 + mm1);
			abm_var[i] = abm_var_core[i] * abm_var[TNBC::Tumor::TUMEX_CC]/abm_var_core[TNBC::Tumor::TUMEX_CC];
			std::cout << "CX_QSP:" << abm_var[i] << std::endl;
			}
		}
		else if (i==TNBC::Tumor::TUMEX_TEXH){

			if (texh_core == 0){
				abm_var[i] = 0;
			}	
			else if (abm_var_core[i] == 1){		
				abm_var[i] = abm_var[TNBC::Tumor::TUMEX_CC]/abm_var_core[TNBC::Tumor::TUMEX_CC];
			}
			else{
				double mm1 =  (abm_var_core[i]-texh_core)/texh_core - abm_var_core[TNBC::Tumor::TUMEX_TEXH_D]/texh_core;
			//abm_var[i] = qsp_var[TNBC::LymphCentral::QSPEX_TUM_TEXH]*AVOGADROS * (1 + mm1);
			abm_var[i] = abm_var_core[i] * abm_var[TNBC::Tumor::TUMEX_CC]/abm_var_core[TNBC::Tumor::TUMEX_CC];			
			std::cout << "TEXH_QSP:" << abm_var[i] << std::endl;
			}
		}

		}


		/* QSP time step */
		_lymph.time_step(t0, dt);
		//std::cout << "RNG check (" << slice << ") QSP: " << rng.get_unif_01() << std::endl;

		_lymph.update_qsp_var(abm_var);
	}

}
	else{
		_lymph.time_step(t0, dt);
	}

	return;

}

void TNBC_Core::write_stats_header(void) const {

	for (int i = 0; i < ic.getVal(IC_NUM_ROI_core); i++)
	{
		//std::cout << "header: core: " << i << std::endl;
		auto& statsStream = output_hub.getStatsFstream(true, i);
		statsStream << _ROI_core[i]->get_stats().writeHeader();
		statsStream.flush();
	}
	for (int i = 0; i < ic.getVal(IC_NUM_ROI_margin); i++)
	{
		//std::cout << "header: margin: " << i << std::endl;
		auto& statsStream = output_hub.getStatsFstream(false, i);
		statsStream << _ROI_margin[i]->get_stats().writeHeader();
		statsStream.flush();
	}
	return;
}

void TNBC_Core::write_stats_slice(unsigned long slice)const{
	for (int i = 0; i < ic.getVal(IC_NUM_ROI_core); i++)
	{
		auto& statsStream = output_hub.getStatsFstream(true, i);
		statsStream << _ROI_core[i]->get_stats().writeSlice(slice);
		statsStream.flush();
	}
	for (int i = 0; i < ic.getVal(IC_NUM_ROI_margin); i++)
	{
		auto& statsStream = output_hub.getStatsFstream(false, i);
		statsStream << _ROI_margin[i]->get_stats().writeSlice(slice);
		statsStream.flush();
	}
	return;
}


void TNBC_Core::write_QSP(unsigned long slice, bool header)const{
	auto& stream = output_hub.get_lymph_blood_QSP_stream();
	if (header){
		stream << "time" << _lymph.getQSPHeaders() << std::endl;
	}
	else{
		stream << slice << _lymph << std::endl;
	}
	return;
}


void TNBC_Core::briefStats(unsigned long slice){
	std::cout << "Time: " << slice << std::endl;
	for (int i = 0; i < ic.getVal(IC_NUM_ROI_core); i++)
	{
		const auto& stats = _ROI_core[i]->get_stats();
		std::cout << "Core_" << i << ": " << "nrCell: " << _ROI_core[i]->getNrCell()
			<< ", CD8 (eff+cyt): " << stats.getTCell()
			<< ", CD8 (supp): " << stats.getTexh()
			<< ", Treg: " << stats.getTreg()
			<< ", MDSC: " << stats.getMDSC()
			<< ", Cancer cell:" << stats.getCancerCell() << std::endl;
	}
	for (int i = 0; i < ic.getVal(IC_NUM_ROI_margin); i++)
	{
		const auto& stats = _ROI_margin[i]->get_stats();
		std::cout << "Margin_" << i << ": " << "nrCell: " << _ROI_margin[i]->getNrCell()
			<< ", CD8 (eff+cyt): " << stats.getTCell()
			<< ", CD8 (supp): " << stats.getTexh()
			<< ", Treg: " << stats.getTreg()
			<< ", MDSC: " << stats.getMDSC()
			<< ", Cancer cell:" << stats.getCancerCell() << std::endl;
	}
}

/*! Print grid info to file.
    \param [in] slice
	\param [in] option: 1. only cellular scale; 2. only molecular scale; 3. both scales
*/
void TNBC_Core::writeGrids(unsigned long slice, unsigned int option){
	if (option == 1 || option == 3)
	{
		for (int i = 0; i < ic.getVal(IC_NUM_ROI_core); i++)
		{
			std::string prefix = "cell_core_" + std::to_string(i) + "_";
			std::ofstream& snap = output_hub.getNewGridToSnapshotStream(slice, prefix);
			snap << _ROI_core[i]->compartment_cells_to_string();
			snap.close();
		}
		for (int i = 0; i < ic.getVal(IC_NUM_ROI_margin); i++)
		{
			std::string prefix = "cell_margin" + std::to_string(i) + "_";
			std::ofstream& snap = output_hub.getNewGridToSnapshotStream(slice, prefix);
			snap << _ROI_margin[i]->compartment_cells_to_string();
			snap.close();
		}
	}
	if (option == 2 || option == 3)
	{
		for (int i = 0; i < ic.getVal(IC_NUM_ROI_core); i++)
		{
			std::string prefix = "grid_core_" + std::to_string(i) + "_";
			std::ofstream&  snap = output_hub.getNewGridToSnapshotStream(slice, prefix);
			snap << _ROI_core[i]->printGridToFile();
			snap.close();

		}
		for (int i = 0; i < ic.getVal(IC_NUM_ROI_margin); i++)
		{
			std::string prefix = "grid_margin" + std::to_string(i) + "_";
			std::ofstream&  snap = output_hub.getNewGridToSnapshotStream(slice, prefix);
			snap << _ROI_margin[i]->printGridToFile();
			snap.close();
		}
	}
}

