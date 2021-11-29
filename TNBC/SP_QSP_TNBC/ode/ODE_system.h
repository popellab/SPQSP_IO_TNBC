#ifndef __CancerVCT_ODE__
#define __CancerVCT_ODE__

#include "SP_QSP_shared/Numerical_Adaptor/CVODE/CVODEBase.h"
#include "Param.h"

namespace CancerVCT{

class ODE_system :
	public CVODEBase
{
public:
    //! ODE right hand side
    static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
    //! Root finding (for events)
    static int g(realtype t, N_Vector y, realtype *gout, void *user_data);
    static std::string getHeader();
    static void setup_class_parameters(Param& param);

    // accessing class parameters
    static double get_class_param(unsigned int i);
    static void set_class_param(unsigned int i, double v);

    template<class Archive>
    static void classSerialize(Archive & ar, const unsigned int  version);

	static bool use_steady_state;
	static bool use_resection;
    static double _QSP_weight;
private:
    // parameters shared by all instances of this class
    static state_type _class_parameter;
public:
    ODE_system();
    ODE_system(const ODE_system& c);
    ~ODE_system();
   
   // read tolerance from parameter file
   void setup_instance_tolerance(Param& param);
   // read variable values from parameter file
   void setup_instance_varaibles(Param& param);
   // evaluate intial assignment
   void eval_init_assignment(void);

   unsigned int get_num_variables(void)const { return _species_var.size(); };   

protected:
    void setupVariables(void);
    void setupEvents(void);
    void initSolver(realtype t0);
    void update_y_other(void);
    bool triggerComponentEvaluate(int i, realtype t, bool curr);
    //! evaluate one event trigger
    bool eventEvaluate(int i);
    //! execute one event
    bool eventExecution(int i, bool delay, realtype& dt);
    //! unit conversion factor for species (y and non-y)
    realtype get_unit_conversion_species(int i) const;
    //! unit conversion foactor for non-species variables
    realtype get_unit_conversion_nspvar(int i) const;
    bool allow_negative(int i)const{ return false; };
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/);
};

template<class Archive>
inline void ODE_system::serialize(Archive & ar, const unsigned int /* version */){
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CVODEBase);
}

template<class Archive>
void ODE_system::classSerialize(Archive & ar, const unsigned int /* version */){
    ar & BOOST_SERIALIZATION_NVP(_class_parameter);
    ar & BOOST_SERIALIZATION_NVP(use_steady_state);
    ar & BOOST_SERIALIZATION_NVP(use_resection);
    ar & BOOST_SERIALIZATION_NVP(_QSP_weight);
}


inline double ODE_system::get_class_param(unsigned int i){
    if (i < _class_parameter.size())
        return _class_parameter[i];
    else
        throw std::invalid_argument("Accessing ODE class parameter: out of range");
}
inline void ODE_system::set_class_param(unsigned int i, double v){
    if (i < _class_parameter.size())
        _class_parameter[i] = v;
    else
        throw std::invalid_argument("Assigning ODE class parameter: out of range");
}

};
#endif
