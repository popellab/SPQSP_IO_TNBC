//#include <boost/serialization/export.hpp>
#include "TReg.h" 

//BOOST_CLASS_EXPORT_IMPLEMENT(Treg)

#include <iostream>
#include <sstream>

#include "../../core/GlobalUtilities.h"
#include "../compartment/Tumor.h"
#include "TCell.h"

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

Treg::Treg(SpatialCompartment* c)
	:Cell_Tumor(c)
	, _divide_cd_Treg_exp(0)
	, _divide_limit_Treg_exp(params.getVal(PARAM_TREG_EXP_DIV_LIMIT))				
	//, _source_IL_10(NULL)
{
	_life = getTregLife();
}

Treg::Treg(const Treg& c)
	:Cell_Tumor(c)
	, _divide_cd_Treg_exp(c._divide_cd_Treg_exp)
	, _divide_limit_Treg_exp(c._divide_limit_Treg_exp)			
	//, _source_IL_10(NULL)
{
	_life = getTregLife();
	//setup_chem_source(_source_IL_10, CHEM_IL_10, params.getVal(PARAM_IL_10_RELEASE));
}

Treg::~Treg()
{
}

std::string Treg::toString()const{
	std::stringstream ss;
	ss << Cell_Tumor::toString();
	return ss.str();
}

bool Treg::agent_movement_step(double t, double dt, Coord& c){
	bool move = false;
	/**/
	if (rng.get_unif_01() < params.getVal(PARAM_TREG_MOVE_PROB))
	{
		// move
		int idx;
		const auto shape = getCellShape();
		if (_compartment->getOneOpenVoxel(shape->getMoveDestinationVoxels(), 
			shape->getMoveDirectionAnchor(), _coord, getType(), idx, rng))
		{
			c = getCellShape()->getMoveDirectionAnchor()[idx] + _coord;				
				move = true;
		}
	}
	return move;
}

bool Treg::agent_state_step(double t, double dt, Coord& c, int cc_p, int teff_p, int treg_p, int mdsc_p, int texh_p){
	bool divide = false;
	if (!isDead())
	{
		_life--;
		if (_life == 0)
		{			
			setDead();
			// remove source when cell die
			return divide;
		}
	}

	const auto shape = getCellShape();
	Cell_Tumor::agent_state_step(t, dt, c, cc_p, teff_p, treg_p, mdsc_p, texh_p);

		auto tumor = dynamic_cast<Tumor*>(_compartment);
		double tumvol = tumor->get_Tum_Vol();
		double tregul = tumor->get_Treg();
		double ArgI = tumor->get_argi();
		
	//double ArgI = get_tumor().get_chem(c, CHEM_ARGI);

	if (ArgI > 0)
	{
		if (_divide_cd_Treg_exp > 0)
		{
			_divide_cd_Treg_exp--;
		}

		double prob_exp = (ArgI / (ArgI + params.getVal(PARAM_EC50_ARGI_TREG)) * (1-tregul/(tumvol*params.getVal(PARAM_TREGMAX))))/params.getVal(PARAM_TREG_EXP_INTERVAL_SLICE);
		double p;
		if (prob_exp > 0){
			p = (prob_exp < 1 ? prob_exp : 1);
		}
		else{
			p = 0;
		}


		if (_divide_limit_Treg_exp > 0 && _divide_cd_Treg_exp == 0 && rng.get_unif_01() < p)
		{
			int idx;
			if (_compartment->getOneOpenVoxel(shape->getProlifDestinationVoxels(), 
				shape->getProlifDestinationAnchor(), _coord, getType(), idx, rng))
		{
			c = getCellShape()->getProlifDestinationAnchor()[idx] + _coord;				
				divide = true;
				_divide_limit_Treg_exp -= 1;
				_divide_cd_Treg_exp = int(params.getVal(PARAM_TREG_EXP_INTERVAL_SLICE)/(ArgI / (ArgI + params.getVal(PARAM_EC50_ARGI_TREG)) * (1-tregul/(tumvol*params.getVal(PARAM_TREGMAX)))) + .5);	
							
		}	
		}

	}
	return divide;
}

void Treg::move_all_source_sink(void) const
{
	//move_source_sink(_source_IL_10);
}

int Treg::getTregLife(){

	double lifeMean = params.getVal(PARAM_T_CELL_LIFE_MEAN_SLICE);
	//double tLifeD = lifeMean;
	double rn = rng.get_unif_01();
	double tLifeD = lifeMean*std::log(1/rn);

	int tLife = int(tLifeD + 0.5);
	tLife = tLife > 0 ? tLife : 0;
	//std::cout << "random Treg life: " << tLife << std::endl;
	return tLife;
}

};
};