#pragma once
#include <mechanical_model.h>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>
#include <polynomial_splines.h>
#include <cppad/cg/model/functor_generic_model.hpp>

class Cost 
{
    protected:
        ADFun<Scalar> f_;
        ADFun<CGScalar> cg_f_;
        std::unique_ptr<LinuxDynamicLib<Scalar>> cg_lib_;
        std::unique_ptr<FunctorGenericModel<double>> cg_model_forward_;
        std::unique_ptr<FunctorGenericModel<double>> cg_model_jacobian_;

    public:
        Cost();
        virtual ~Cost() = default;

        /**
         * This method initializes the AD function with correct calculations sequence
        */
        virtual void initialize_AD_Function(MechanicalModel & mech_model) = 0;
        /**
         * This method returns the Jacobian (stored in a vector) given a point of the state space 
         * and a vector of parameters
         * @param state State vector
         * @param param Vector of parameters
        */
        Vector jacobian(Vector const& state,
                        Vector const& param
                        );

        /**
         * This method returns the Jacobian (stored in a vector) given a point of the state space 
         * and a vector of parameters
         * @param state State vector
         * @param param Vector of parameters
         * @param lambda Vector of lambda parameters of hessian
        */
        Vector hessian(Vector const& state,
                       Vector const& param,
                       Vector const& lambda
                        );

        /**
         * This method returns the Jacobian (stored in a vector) given a point of the state space 
         * and a vector of parameters
         * @param state State vector
         * @param param Vector of parameters
         * @param sigma Value of the i-th hessian desired
        */
        Vector hessian(Vector const& state,
                       Vector const& param,
                       int const& sigma
                        );

        /**
         * This method returns the Jacobian (stored in a vector) given a point of the state space 
         * and a vector of parameters. It uses codegen
         * @param state State vector
         * @param param Vector of parameters
        */
        virtual Vector jacobian_cg(Vector const& state,
                                    Vector const& param
                                    );

        virtual Vector jacobian_cg(Vector const& state,
                                    Vector const& param,
                                    int const& dep_vars
                                    );

        /**
         * Returns the residual given a point of the state space and a parameters vector
         * @param state State vector
         * @param param Vector of parameters
        */
        Vector forward(Vector const& state,
                        Vector const& param
                        );

        /**
         * Returns the residual given a point of the state space and a parameters vector.
         * It uses codegen
         * @param state State vector
         * @param param Vector of parameters
        */
        virtual Vector forward_cg(Vector const& state,
                            Vector const& param
                            );

        virtual Vector forward_cg(Vector const& state,
                            Vector const& param,
                            int const& dep_vars
                            );

        virtual void compile_CG_lib(MechanicalModel & mech_model) = 0;
        virtual void load_CG_lib() = 0;

        /**
         * Returns the cppad function by reference
        */
        ADFun<Scalar>& f();
        // void set_Mechanical_Model(MechanicalModel* model);
        // MechanicalModel* mechanical_Model();
};

class Constr_HT_filip_DOC : public Cost
{
    public:
    Constr_HT_filip_DOC();
    Constr_HT_filip_DOC(MechanicalModel &mech_model);

    virtual ~Constr_HT_filip_DOC();

    void initialize_AD_Function(MechanicalModel & mech_model) override;
    void compile_CG_lib(MechanicalModel & mech_model) override;
    virtual void load_CG_lib();
};

class Cost_DOC_filip : public Cost
{
    public:
    Cost_DOC_filip();
    Cost_DOC_filip(MechanicalModel &mech_model);

    virtual ~Cost_DOC_filip();

    void initialize_AD_Function(MechanicalModel & mech_model) override;
    void compile_CG_lib(MechanicalModel & mech_model) override;
    virtual void load_CG_lib();
};