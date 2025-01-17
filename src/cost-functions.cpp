#include <cost-functions.h>

/**Cost abstract class**/
Cost::Cost()
{};

Vector Cost::jacobian(Vector const& state,
                      Vector const& param)
{
        f_.new_dynamic(param);
        return f_.Jacobian(state);
}

Vector Cost::jacobian_cg(Vector const& state,
                            Vector const& param
                            )
{
        Vector X(state.rows() + param.rows());
        X.head(state.rows()) = state;
        X.tail(param.rows()) = param;
        // std::cout << "m : " << cg_model_jacobian_->Range() << std::endl;
        return cg_model_jacobian_->ForwardZero<Vector>(X);
}

Vector Cost::jacobian_cg(Vector const& state,
                        Vector const& param,
                        int const& dep_vars
                        )
{
        Vector X(state.rows() + param.rows());
        X.head(state.rows()) = state;
        X.tail(param.rows()) = param;
        Vector dep(Vector::Zero(dep_vars));
        cg_model_jacobian_->ForwardZero(X, dep);
        return dep;
}


Vector Cost::forward(Vector const& state,
                        Vector const& param
                        )
{       
        f_.new_dynamic(param);
        return f_.Forward(0, state);
}

Vector Cost::forward_cg(Vector const& state,
                        Vector const& param
                        )
{
        Vector X(state.rows() + param.rows());
        X.head(state.rows()) = state;
        X.tail(param.rows()) = param;
        return cg_model_forward_->ForwardZero<Vector>(X);
}

Vector Cost::forward_cg(Vector const& state,
                        Vector const& param,
                        int const& dep_vars
                        )
{
        Vector X(state.rows() + param.rows());
        X.head(state.rows()) = state;
        X.tail(param.rows()) = param;
        Vector dep(Vector::Zero(dep_vars));
        cg_model_forward_->ForwardZero(X, dep);
        return dep;
}

Vector Cost::hessian(Vector const& state,
                Vector const& param,
                Vector const& lambda)
{
        f_.new_dynamic(param);
        return f_.Hessian(state, lambda);
}

Vector Cost::hessian(Vector const& state,
                    Vector const& param,
                    int const& sigma)
{
    f_.new_dynamic(param);
    return f_.Hessian(state, sigma);
}

/**constr_DOC class with human to box constraints **/
Constr_HT_filip_DOC::Constr_HT_filip_DOC() {};
Constr_HT_filip_DOC::Constr_HT_filip_DOC(MechanicalModel &mech_model) 
{
    // initialize_AD_Function(mech_model); 
    compile_CG_lib(mech_model);  
};

Constr_HT_filip_DOC::~Constr_HT_filip_DOC() {};

void Constr_HT_filip_DOC::initialize_AD_Function(MechanicalModel & mech_model)
{
};
void Constr_HT_filip_DOC::compile_CG_lib(MechanicalModel & mech_model)
{
    if(!system::isFile("../cg_libs/constr_cg_func.so"))
    {
        std::cout << "No cg lib found for constraints function " << std::endl;
    }else 
    {
        // std::cout << "Existing ./constr_cg_func found" << std::endl;
        load_CG_lib();
    }

};

void Constr_HT_filip_DOC::load_CG_lib()
{
        cg_lib_ = std::make_unique<LinuxDynamicLib<Scalar>>("../cg_libs/constr_cg_func.so");
        cg_model_forward_ = cg_lib_->modelFunctor("forward_fun");
        cg_model_jacobian_ = cg_lib_->modelFunctor("jac_fun");
};



/**Cost_DOC class **/
Cost_DOC_filip::Cost_DOC_filip() {};
Cost_DOC_filip::Cost_DOC_filip(MechanicalModel &mech_model) 
{
    // initialize_AD_Function(mech_model);   
    compile_CG_lib(mech_model);
};

Cost_DOC_filip::~Cost_DOC_filip() {};

void Cost_DOC_filip::initialize_AD_Function(MechanicalModel & mech_model)
{
};

void Cost_DOC_filip::compile_CG_lib(MechanicalModel & mech_model)
{
    if(!system::isFile("../cg_libs/cost_cg_func.so"))
    {
        std::cout << "No cg lib found for cost function " << std::endl;
    }else 
    {
        // std::cout << "Existing ./cost_cg_func found" << std::endl;
        load_CG_lib();
    }
};

void Cost_DOC_filip::load_CG_lib()
{
        cg_lib_ = std::make_unique<LinuxDynamicLib<Scalar>>("../cg_libs/cost_cg_func.so");
        cg_model_forward_ =  cg_lib_->modelFunctor("forward_fun");
        cg_model_jacobian_ =  cg_lib_->modelFunctor("jac_fun");
};
