#include "VoxelContentGen.h"
#include "../../core/GlobalUtilities.h"

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

VoxelContentGen::VoxelContentGen()
	: _stationary(true)
	, _x_min(0)
	, _y_min(0)
	, _z_min(0)
	, _x_max(0)
	, _y_max(0)
	, _z_max(0)
	, _celltype_cdf()
{
	int dmax = params.getVal(PARAM_CANCER_CELL_PROGENITOR_DIV_MAX);
	// 0: stem; 1-dmax: progenitor; dmax+1: senescent; dmax+2: empty
	_celltype_cdf = std::vector<double>(dmax + 3, 0.0);
	return;
}


VoxelContentGen::~VoxelContentGen()
{
}

/*! return true if any cell is to be created
*/
bool VoxelContentGen::get_type_state(const Coord3D& c, RNG& rng,
	AgentType& type, AgentState& state, int&div)const{

	bool create_cell = false;
	int dmax = params.getVal(PARAM_CANCER_CELL_PROGENITOR_DIV_MAX);
	type = AgentTypeEnum::CELL_TYPE_CANCER;
	state = AgentStateEnum::CANCER_PROGENITOR;
	

	//Total tumor (tumor growth)
	//params.getVal(PARAM_MEAN_INIT_NORM_DIST) = 0, params.getVal(PARAM_SD_INIT_NORM_DIST) = 0.99
	//double mean_dist = params.getVal(PARAM_MEAN_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//double sdx_dist = params.getVal(PARAM_SD_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//double sdy_dist = params.getVal(PARAM_SD_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//double sdz_dist = params.getVal(PARAM_SD_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;

	//Total tumor

	//double mean_dist = params.getVal(PARAM_MEAN_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//double sdx_dist = params.getVal(PARAM_SD_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//double sdy_dist = params.getVal(PARAM_SD_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//double sdz_dist = params.getVal(PARAM_SD_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;


	//Partial tumor (3D)
	//double mean_dist = params.getVal(PARAM_MEAN_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//double sd_dist = params.getVal(PARAM_SD_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//double dx2cc = std::pow(c.x,2);
	//double dy2cc = std::pow(c.y,2);
	//double dz2cc = std::pow(c.z,2);
	//double dist_sourcecc = std::pow(dx2cc+dy2cc+dz2cc, 0.5)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//if (dist_sourcecc <= rng.get_normal(mean_dist,sd_dist))	

	//Partial tumor (thin slice)
	double mean_dist = params.getVal(PARAM_MEAN_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	double sd_dist = params.getVal(PARAM_SD_INIT_NORM_DIST)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	double dx2cc = std::pow(c.x,2);
	double dy2cc = std::pow(c.y,2);
	double dist_sourcecc = std::pow(dx2cc+dy2cc, 0.5)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	if ((dist_sourcecc <= rng.get_normal(mean_dist,sd_dist)) && (c.z >= 0 && c.z < params.getVal(PARAM_TUMOR_Z)))	


	//Entire tumor
	//double dx2cc = std::pow(c.x,2);
	//double dy2cc = std::pow(c.y,2);	
	//double dz2cc = std::pow(c.z,2);	
	//double dx2cc = std::pow(c.x-params.getVal(PARAM_TUMOR_X)/2,2);
	//double dy2cc = std::pow(c.y-params.getVal(PARAM_TUMOR_Y)/2,2);	
	//double dz2cc = std::pow(c.z-params.getVal(PARAM_TUMOR_Z)/2,2);		
	//double dist_sourceccx = std::pow(dx2cc, 0.5)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//double dist_sourceccy = std::pow(dy2cc, 0.5)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//double dist_sourceccz = std::pow(dz2cc, 0.5)*params.getVal(PARAM_VOXEL_SIZE)/1e6;
	//if ((dist_sourceccx <= rng.get_normal(mean_dist,sdx_dist))
	//&& (dist_sourceccy <= rng.get_normal(mean_dist,sdy_dist))
	//&& (dist_sourceccz <= rng.get_normal(mean_dist,sdz_dist)))	
		
	{
		int i = rng.sample_cdf(_celltype_cdf);
		if (i <= dmax+1)
		{
			create_cell = true;
			if (i==0)
			{
				state = AgentStateEnum::CANCER_STEM;
			}
			else if (i==dmax+1)
			{
				state = AgentStateEnum::CANCER_SENESCENT;
			}
			else{
				state = AgentStateEnum::CANCER_PROGENITOR;
				div = dmax + 1 - i;
			}
		}
	}
	return create_cell;
}

void VoxelContentGen::setup(bool stationary, double cancer_prob,
	int xlim, int ylim, int zlim){
	_stationary = stationary;

	_x_min= 0;
	_y_min= 0;
	_z_min= 0;
	_x_max= _x_min + xlim;
	_y_max= _y_min + ylim;
	_z_max= _z_min + zlim;
	unsigned int dmax = params.getVal(PARAM_CANCER_CELL_PROGENITOR_DIV_MAX);
	double k, r, rs, rp, mu, l0, l1, l2;
	k = params.getVal(PARAM_CANCER_STEM_ASYMMETRIC_DIV_PROB);
	rs = params.getVal(PARAM_CSC_GROWTH_RATE);
	rp = params.getVal(PARAM_CANCER_PROG_GROWTH_RATE);
	mu = params.getVal(PARAM_CANCER_SENESCENT_DEATH_RATE);
	r = rs * (1 - k);
	l0 = k*rs / (r + rp);
	l1 = 2 * rp / (r + rp);
	l2 = 2 * rp / (r + mu);
	double C;
	if (l1==1){
		C = 1 + l0 + l0*l2*std::pow(l1,(dmax - 1));
	}
	else{
		C = 1 + l0*(std::pow(l1, dmax) - 1) / (l1 - 1) + l0*l2*std::pow(l1,(dmax - 1));
	}
	double p ;
	_celltype_cdf[0] = p = 1 / C * cancer_prob; // joint P
	p *= l0;
	_celltype_cdf[1] = _celltype_cdf[0] + p;

	for (size_t i = 2; i <= dmax; i++)
	{
		p *= l1;
		_celltype_cdf[i] = _celltype_cdf[i - 1] + p;	
	}
	p *= l2;
	_celltype_cdf[dmax + 1] = _celltype_cdf[dmax] + p;
	_celltype_cdf[dmax + 2] = 1.0;
	/*
	std::cout << "cancer cell cdf:" << std::endl;
	for (size_t i = 0; i < _celltype_cdf.size(); i++)
	{
		std::cout << i << ", " << _celltype_cdf[i] << std::endl;
	}
	*/
	return;
}

/* Fill (0,0,0), (xlim, ylim, zlim) with cancer cells
*/
void VoxelContentGen::setup( double pstem, int xlim, int ylim, int zlim,
	int x0 = 0, int y0=0, int z0 = 0) {

	_stationary = false;

	_x_min= x0;
	_y_min= y0;
	_z_min= z0;
	_x_max= _x_min + xlim;
	_y_max= _y_min + ylim;
	_z_max= _z_min + zlim;

	unsigned int dmax = params.getVal(PARAM_CANCER_CELL_PROGENITOR_DIV_MAX);

	_celltype_cdf[0] = pstem;
	_celltype_cdf[1] = 1.0;
	for (size_t i = 2; i <= dmax+2; i++)
	{
		_celltype_cdf[i] = 1.0;
	}
	/*
	std::cout << "cancer cell cdf:" << std::endl;
	for (size_t i = 0; i < _celltype_cdf.size(); i++)
	{
		std::cout << i << ", " << _celltype_cdf[i] << std::endl;
	}
	*/
	return;
}
};// end of namespace
};