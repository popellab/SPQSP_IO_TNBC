//#include <boost/serialization/export.hpp>
#include "MDSC.h" 

//BOOST_CLASS_EXPORT_IMPLEMENT(MDSC)

#include <iostream>
#include <sstream>

#include "../../core/GlobalUtilities.h"
#include "../compartment/Tumor.h"
#include "../compartment/LymphCentral.h"
//#include "TCell.h"

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

MDSC::MDSC(SpatialCompartment* c)
	:Cell_Tumor(c)
	//, _source_IL_10(NULL)
	, _source_ArgI(NULL)
	, _source_NO(NULL)    
{
	_life = getMDSCLife();
}

MDSC::MDSC(const MDSC& c)
	:Cell_Tumor(c)
	//, _source_IL_10(NULL)
	, _source_ArgI(NULL)
	, _source_NO(NULL)      
{
	_life = getMDSCLife();
	setup_chem_source(_source_ArgI, CHEM_ARGI, params.getVal(PARAM_ARGI_RELEASE));
	setup_chem_source(_source_NO, CHEM_NO, params.getVal(PARAM_NO_RELEASE));	
	//setup_chem_source(_source_IL_10, CHEM_IL_10, params.getVal(PARAM_IL_10_RELEASE));
}

MDSC::~MDSC()
{
}

std::string MDSC::toString()const{
	std::stringstream ss;
	ss << Cell_Tumor::toString();
	return ss.str();
}


bool MDSC::agent_movement_step(double t, double dt, Coord& c){
	bool move = false;
	/**/
	if (rng.get_unif_01() < params.getVal(PARAM_MDSC_MOVE_PROB))
	{
		// move
		int idx;
		const auto shape = getCellShape();
		if (_compartment->getOneOpenVoxel(shape->getMoveDestinationVoxels(), 
			shape->getMoveDirectionAnchor(), _coord, getType(), idx, rng))
		{

			double c1 = 0;
			double c2 = 0;
			double c3 = 0;
			double alpha1 = (c.x-0)*(c.x-params.getVal(PARAM_TUMOR_X)+1);
			double alpha2 = (c.y-0)*(c.y-params.getVal(PARAM_TUMOR_Y)+1);
			double alpha3 = (c.z-0)*(c.z-params.getVal(PARAM_TUMOR_Z)+1);

			if (alpha1+alpha2+alpha3 == 0 && rng.get_unif_01()<c1){
				move = false;
			}
			else if(alpha1*alpha2+alpha1*alpha3+alpha2*alpha3 == 0 && rng.get_unif_01()<c2){
				move = false;
			}
			else if(alpha1*alpha2*alpha3 == 0 && rng.get_unif_01()<c3){
				move = false;
			}
			else{
			c = getCellShape()->getMoveDirectionAnchor()[idx] + _coord;				
				move = true;
			}
		}
	}
	return move;
}

bool MDSC::agent_state_step(double t, double dt, Coord& c, int cc_p, int teff_p, int treg_p, int mdsc_p, int texh_p){
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

	return divide;
}


int MDSC::getMDSCLife(){

	double lifeMean = params.getVal(PARAM_MDSC_LIFE_MEAN);
	//double tLifeD = rng.get_exponential(lifeMean);
	double rn = rng.get_unif_01();
	double tLifeD = lifeMean*std::log(1/rn);
	int tLife = int(tLifeD + 0.5);
	tLife = tLife > 0 ? tLife : 0;
	//std::cout << "random MDSC life: " << tLife << std::endl;
	return tLife;
}

//! move sources (NO and ArgI)
void MDSC::move_all_source_sink(void)const{
	//std::cout << "moving sources: " << _coord << std::endl;
	move_source_sink(_source_ArgI);
	move_source_sink(_source_NO);  
	return;
}

//! remove sources (NO and ArgI)
void MDSC::remove_all_source_sink(void){
	remove_source_sink(_source_ArgI);
	remove_source_sink(_source_NO); 
	return;
}

};
};