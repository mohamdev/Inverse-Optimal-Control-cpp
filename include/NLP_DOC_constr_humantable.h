#pragma once
#include <cost-functions.h>
#include <IpTNLP.hpp>

using namespace Ipopt;

class DOC_NLP_HTconstr: public TNLP
{   
    private:

    /**@name Methods to block default compiler methods.
        *
        * The compiler automatically generates the following three methods.
        *  Since the default compiler implementation is generally not what
        *  you want (for all but the most simple classes), we usually
        *  put the declarations of these methods in the private section
        *  and never implement them. This prevents the compiler from
        *  implementing an incorrect "default" behavior without us
        *  knowing. (See Scott Meyers book, "Effective C++")
        */
    //@{
    DOC_NLP_HTconstr(
        const DOC_NLP_HTconstr&
    );

    DOC_NLP_HTconstr& operator=(
        const DOC_NLP_HTconstr&
    );
    //@}


    public:
    Cost * cost_func_ptr_; //Pointer to Cost Function
    Cost * constr_func_ptr_; //Pointer to Cost Function
    MechanicalModel mech_model_; //Mechanical model
    NlpParam nlp_param_; //Parameters of the non-linear problem
    Vector q_initial_;
    Vector dq_initial_;
    Vector cost_params_;
    Vector constr_params_;
    Vector cp_est_;
    Vector h_sph_constr_val_;
    double t_sph_constr_val_;
    double cop_l_constr_val_;
    double cop_u_constr_val_;
    Vector weights;

    /** Default constructor */
    DOC_NLP_HTconstr();

    /** Constructor with DOC_NLP_HTconstr parameters and urdf model path */
    DOC_NLP_HTconstr(NlpParam const& nlp_param, Vector const& weights_par);

    /** Default destructor */
    virtual ~DOC_NLP_HTconstr();

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the DOC_NLP_HTconstr */
    virtual bool get_nlp_info(
        Index&          n,
        Index&          m,
        Index&          nnz_jac_g,
        Index&          nnz_h_lag,
        IndexStyleEnum& index_style
    );

    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(
        Index   n,
        Number* x_l,
        Number* x_u,
        Index   m,
        Number* g_l,
        Number* g_u
    );

    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(
        Index   n,
        bool    init_x,
        Number* x,
        bool    init_z,
        Number* z_L,
        Number* z_U,
        Index   m,
        bool    init_lambda,
        Number* lambda
    );

    /** Method to return the objective value */
    virtual bool eval_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number&       obj_value
    );

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number*       grad_f
    );

    /** Method to return the constraint residuals */
    virtual bool eval_g(
        Index         n,
        const Number* x,
        bool          new_x,
        Index         m,
        Number*       g
    );

    /** Method to return:
        *   1) The structure of the jacobian (if "values" is NULL)
        *   2) The values of the jacobian (if "values" is not NULL)
        */
    virtual bool eval_jac_g(
        Index         n,
        const Number* x,
        bool          new_x,
        Index         m,
        Index         nele_jac,
        Index*        iRow,
        Index*        jCol,
        Number*       values
    );

    /** Method to return:
        *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
        *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
        */
    virtual bool eval_h(
        Index         n,
        const Number* x,
        bool          new_x,
        Number        obj_factor,
        Index         m,
        const Number* lambda,
        bool          new_lambda,
        Index         nele_hess,
        Index*        iRow,
        Index*        jCol,
        Number*       values
    );

    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(
        SolverReturn               status,
        Index                      n,
        const Number*              x,
        const Number*              z_L,
        const Number*              z_U,
        Index                      m,
        const Number*              g,
        const Number*              lambda,
        Number                     obj_value,
        const IpoptData*           ip_data,
        IpoptCalculatedQuantities* ip_cq
    );
    //@}

    /**
     * Returns the cost function pointer
    */
    Cost * cost_Func();

    /**
     * Returns the Non Linear Porblem parameters by reference
    */
    NlpParam& nlp_Param();

    /**
     * Returns the mechanical model by reference
    */
    MechanicalModel& mechanical_Model();

    /**
     * Returns a pointer to the mechanical model
    */
    MechanicalModel* mechanical_Model_Ptr();


    /**
     * Sets the non linear problem parameters
    */
    void set_Nlp_Param(NlpParam const& param);
    
};