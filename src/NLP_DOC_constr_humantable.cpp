#include <NLP_DOC_constr_humantable.h>
#include <cassert>
#include <iostream>


/**
 * **********************************************************************************************
 * **********************************************************************************************
 * ---------------------------------------- CONSTRUCTORS ----------------------------------------
 * **********************************************************************************************
 * **********************************************************************************************
*/
// constructor
DOC_NLP_HTconstr::DOC_NLP_HTconstr() 
//   cost_func_ptr_(nullptr)
{ }

DOC_NLP_HTconstr::DOC_NLP_HTconstr(NlpParam const& nlp_param,
                                   Vector const& weights_par) :
                  nlp_param_(nlp_param),
                  weights(weights_par),
                  mech_model_(),
                  cost_func_ptr_(nullptr)
             
{
    /*Build the mechanical model and cost/constr functions*/
    mech_model_.build_Pinocchio_Model();
    cost_func_ptr_ = new Cost_DOC_filip(mech_model_);
    constr_func_ptr_ = new Constr_HT_filip_DOC(mech_model_);

    /** Get the radius of segments spheres*/
    double shank_sph_rad = 0.4281073446*0.2;
    double thigh_sph_rad = 0.4271186441*0.2;
    double pelv_sph_rad = 0.2412429379*0.2;
    double torso_sph_rad = 0.3302259887*0.2;
    double upa_sph_rad = 0.2738700565*0.2;

    /**CoP constraints */
    cop_l_constr_val_ = -0.03525247161 - 0.3*0.1374293785; //Heel_X - 0.3*foot_l
    cop_u_constr_val_ = 0.1857351821 + 0.3*0.1374293785; //Toe_X + 0.3*foot_l

    double box_mass = 0.0;
    double box_width = 0.0;
    double table_width = 0.0;
    double table_height = 0.0;
    Vector box_pos(Vector::Zero(3,1));
    Vector table_center(Vector::Zero(2,1));
    Vector wrist_final_pos(Vector::Zero(2,1));
    std::vector<double> t_cp;
    std::vector<Vector> q_cp, dq_cp;
    std::string cp_path = "../data/spline_traj_s2_trial_1.txt";
    std::string box_params_path = "../data/s2_trial_1_lifting_environment.txt";
    read_control_points(t_cp, q_cp, dq_cp, cp_path);
    read_box_params(box_mass, box_pos, box_width, box_params_path);
    read_final_hand_pos(wrist_final_pos, box_params_path);
    read_table_params(table_width, table_height, table_center, box_params_path);
    std::vector<Eigen::VectorXd> dq_ref, ddq_ref;
    dq_ref = read_q6_data("../data/dq_ref.txt");
    ddq_ref = read_q6_data("../data/ddq_ref.txt");

    Eigen::VectorXd dq_in = dq_ref[0]; Eigen::VectorXd ddq_in = ddq_ref[0];
    Eigen::VectorXd dq_fin = dq_ref[dq_ref.size()-1]; Eigen::VectorXd ddq_fin = ddq_ref[ddq_ref.size()-1];

   //  std::cout << "Box mass : " << box_mass << std::endl;
   //  std::cout << "Box pos : " << box_pos.transpose() << std::endl;
   //  std::cout << "Box width : " << box_width << std::endl;
   //  std::cout << "wrist_final_pos : " << wrist_final_pos.transpose() << std::endl;
   //  std::cout << "q init : " << q_cp[0].transpose() << std::endl;
   //  std::cout << "dq init : " << dq_cp[0].transpose() << std::endl;
    double box_sph_rad = box_width*0.2;

    /* Vector of constraints params*/
    constr_params_ = Vector::Zero(10+1+3+2+2+2+6+6, 1); //Vector of constr_params
    for(int i = 0; i < t_cp.size(); i++)
    {
        constr_params_(i) = t_cp[i];
    }
    constr_params_(10) = box_mass;
    constr_params_.segment(11,3) = box_pos;
    constr_params_.segment(14,2) = wrist_final_pos;
    constr_params_(16) = table_width;
    constr_params_(17) = table_height;
    constr_params_.segment(18,2) = table_center;
    constr_params_.segment(20,6) = dq_in;
    constr_params_.segment(26,6) = ddq_in;
   //  constr_params_.segment(32,6) = dq_fin;
   //  constr_params_.segment(38,6) = ddq_fin;
   //  constr_params_.segment(14,6) = q_cp[0];
   //  constr_params_.segment(20,6) = dq_cp[0];
   //  constr_params_.segment(26,2) = wrist_final_pos;

    /*Human spheres constraints values*/
    h_sph_constr_val_ = Vector::Zero(5*10);
    for(int i = 0; i < h_sph_constr_val_.rows(); i+=5)
    {
        h_sph_constr_val_(i) = shank_sph_rad + box_sph_rad;
        h_sph_constr_val_(i+1) = thigh_sph_rad + box_sph_rad;
        h_sph_constr_val_(i+2) = pelv_sph_rad + box_sph_rad;
        h_sph_constr_val_(i+3) = torso_sph_rad + box_sph_rad;
        h_sph_constr_val_(i+4) = upa_sph_rad + box_sph_rad;
    }

    /*Table spheres constraints values*/
    t_sph_constr_val_ = table_width/10.0 + box_sph_rad;

    // std::cout << h_sph_constr_val_.transpose() << std::endl;
    /*Vector of const func params*/
    cost_params_ = Vector::Zero(10+7*3+1+3+6+6, 1); //Vector of cost_params
    for(int i = 0; i < t_cp.size(); i++)
    {
        cost_params_(i) = t_cp[i];
    }
    cost_params_(10) = box_mass;
    cost_params_.segment(11,3) = box_pos;
    // cost_params_.tail(4) = Vector::Ones(4,1)*0.1;
   //  Vector weights(Vector::Zero(7*3,1));
   //  weights(0) = 1e-3;
   //  weights(1) = 100.0;
   //  weights(2) = 50.0;
   //  weights(3) = 1e-1;
   //  weights(4) = 0.0;
   //  weights(5) = 100.0;
   //  weights(6) = 150.0;
   //  weights(7) = 1e-3;
   //  weights(8) = 100.0;
   //  weights(9) = 50.0;
   //  weights(10) = 1e-1;
   //  weights(11) = 0.0;
   //  weights(12) = 100.0;
   //  weights(13) = 150.0;
   //  weights(14) = 1e-3;
   //  weights(15) = 100.0;
   //  weights(16) = 50.0;
   //  weights(17) = 1e-1;
   //  weights(18) = 0.0;
   //  weights(19) = 100.0;
   //  weights(20) = 150.0;
    cost_params_.segment(14,7*3) = weights;
    cost_params_.segment(35,6) = dq_in;
    cost_params_.segment(41,6) = ddq_in;
   //  cost_params_.segment(30,6) = dq_in;
   //  cost_params_.segment(36,6) = ddq_in;

   // std::cout << "weights size in : " << weights.rows() << std::endl;
   // std::cout << "weights in : " << weights.transpose() << std::endl;
    q_initial_ = Vector::Zero(60,1);
    q_initial_.segment(0,6) = q_cp[0];
    q_initial_.segment(6,6) = q_cp[1];
    q_initial_.segment(12,6) = q_cp[2];
    q_initial_.segment(18,6) = q_cp[3];
    q_initial_.segment(24,6) = q_cp[4];
    q_initial_.segment(30,6) = q_cp[5];
    q_initial_.segment(36,6) = q_cp[6];
    q_initial_.segment(42,6) = q_cp[7];
    q_initial_.segment(48,6) = q_cp[8];
    q_initial_.segment(54,6) = q_cp[9];
    dq_initial_ = Vector::Zero(60,1);
    dq_initial_.segment(0,6) = dq_cp[0];
    dq_initial_.segment(6,6) = dq_cp[1];
    dq_initial_.segment(12,6) = dq_cp[2];
    dq_initial_.segment(18,6) = dq_cp[3];
    dq_initial_.segment(24,6) = dq_cp[4];
    dq_initial_.segment(30,6) = dq_cp[5];
    dq_initial_.segment(36,6) = dq_cp[6];
    dq_initial_.segment(42,6) = dq_cp[7];
    dq_initial_.segment(48,6) = dq_cp[8];
    dq_initial_.segment(54,6) = dq_cp[9];

   //  std::cout << "q initial : " << q_initial_.transpose() << std::endl;
   //  std::cout << "dq initial : " << dq_initial_.transpose() << std::endl;
}

// destructor
DOC_NLP_HTconstr::~DOC_NLP_HTconstr()
{ 
   // delete cost_func_ptr_;
}

/**
 * **********************************************************************************************
 * **********************************************************************************************
 * ---------------------------------------- IPOPT METHODS ---------------------------------------
 * **********************************************************************************************
 * **********************************************************************************************
*/

// [TNLP_get_nlp_info]
// returns the size of the problem
bool DOC_NLP_HTconstr::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
   // std::cout << "get_nlp_info" << std::endl;
   n = nlp_param_.n; //Number of variables of the problem

   m = nlp_param_.m; //Number of constraints

   nnz_jac_g = nlp_param_.nnz_jac_g;

   nnz_h_lag = nlp_param_.nnz_hes_lag;

   // use the C style indexing (0-based)
   index_style = TNLP::C_STYLE;

   return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool DOC_NLP_HTconstr::get_bounds_info(
   Index   n,
   Number* x_l,
   Number* x_u,
   Index   m,
   Number* g_l,
   Number* g_u
)
{
   // std::cout << "get_bounds_info" << std::endl;
   assert(n == nlp_param_.n); //Number of variables of the problem

   assert(m == nlp_param_.m); //Number of constraints

   for(int i = 0; i<n; i++)
   {
      x_l[i] = nlp_param_.x_l(i);
      x_u[i] = nlp_param_.x_u(i);
   }

   for(int i = 0;  i<m-5*10-10-5*10; i++)
   {
      g_l[i] = g_u[i] = 0.0;
   }
   for(int i = m-5*10-10-5*10; i < m-10-5*10; i++)
   {
       g_l[i] = h_sph_constr_val_(i-(m-5*10-10-5*10));
       g_u[i] = 2e19;
    //    std::cout << g_l[i] << std::endl;
   }
   for(int i = m-10-5*10; i < m-5*10; i++)
   {
       g_l[i] = cop_l_constr_val_;
       g_u[i] = cop_u_constr_val_;
   }
   for(int i = m-5*10; i < m; i++)
   {
       g_l[i] = t_sph_constr_val_;
       g_u[i] = 2e19;
    //    std::cout << g_l[i] << std::endl;
   }
//    for(int i = 0; i < m; i++)
//    {
//        std::cout << "g_l : " << g_l[i] << std::endl;
//        std::cout << "g_u : " << g_u[i] << std::endl;
//    }
   return true;
}
// [TNLP_get_bounds_info]

// [TNLP_get_starting_point]
// returns the initial point for the problem
bool DOC_NLP_HTconstr::get_starting_point(
   Index   n,
   bool    init_x,
   Number* x,
   bool    init_z,
   Number* z_L,
   Number* z_U,
   Index   m,
   bool    init_lambda,
   Number* lambda
)
{
   // std::cout << "get_starting_point" << std::endl;
   assert(n == nlp_param_.n); //Number of variables of the problem

   assert(init_x == true);
   assert(init_z == false);
   assert(init_lambda == false);

   Vector init = Vector::Zero(n);
   init.head(60) = q_initial_;
   for (int i = 0 ; i<n; i++)
   {
      x[i] = init(i);
   }
   return true;
}

// [TNLP_get_starting_point]

// [TNLP_eval_f]
// returns the value of the objective function
bool DOC_NLP_HTconstr::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
      // std::cout << "eval_f" << std::endl;
   assert(n == nlp_param_.n);

   //Store evaluation variables in eval vector
   Vector eval = Vector::Zero(n);
   for (int i = 0; i<n; i++)
   {
      eval(i) = x[i];
   }

   obj_value = (cost_func_ptr_->forward_cg(eval, cost_params_))(0);
   
   return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool DOC_NLP_HTconstr::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{

      // std::cout << "eval_grad_f" << std::endl;
   assert(n == nlp_param_.n);

   //Store evaluation variables in eval vector
   Vector eval = Vector::Zero(n);
   for (int i = 0; i<n; i++)
   {
      eval(i) = x[i];
   }

   Vector jac = cost_func_ptr_->jacobian_cg(eval, cost_params_);
   
   for (int i = 0; i<nlp_param_.n; i++)
   {
      grad_f[i] = jac(i);
   }

   return true;
}
// [TNLP_eval_g]
/** Method to return the constraint residuals */
bool DOC_NLP_HTconstr::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
   assert(n == nlp_param_.n);
   assert(m == nlp_param_.m);
   //Store evaluation variables in eval vector
   Vector eval = Vector::Zero(n);
   for (int i = 0; i<n; i++)
   {
       eval(i) = x[i];
   }

   Vector constr_val = constr_func_ptr_->forward_cg(eval, constr_params_);
   // std::cout << "constraints value : " << constr_val.transpose() << std::endl;
   for(int i = 0; i < m; i++)
   {
      g[i] = constr_val(i);
   }


   return true;
}

// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool DOC_NLP_HTconstr::eval_jac_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Index         nele_jac,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
   assert(n == nlp_param_.n);
   assert(m == nlp_param_.m);

   if( values == NULL )
   {
      int k_th_element=0;
      for(int i = 0; i < m; i++)
      {
         for(int j = 0; j < n; j++)
         {
            iRow[k_th_element] = i;
            jCol[k_th_element] = j;
            k_th_element++;
         }
      }
   }
   else
   {
      Vector eval = Vector::Zero(n);
      for (int i = 0; i<n; i++)
      {
         eval(i) = x[i];
      }
      
      Vector constr_jac = constr_func_ptr_->jacobian_cg(eval, constr_params_);
      for(int i = 0; i < m*n; i++)
      {
         values[i] = constr_jac[i];
      }
   }
   // std::cout << "eval_jac_g" << std::endl;
   return true;
}
// [TNLP_eval_jac_g]

// [TNLP_eval_h]
//return the structure or values of the Hessian
bool DOC_NLP_HTconstr::eval_h(
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
)
{
   assert(n == nlp_param_.n);
   assert(m == nlp_param_.m);

   // if( values == NULL )
   // {
   //    int k_th_element=0;
   //    for(int i = 0; i < n; i++)
   //    {
   //       for(int j = 0; j < n; j++)
   //       {
   //          iRow[k_th_element] = i;
   //          jCol[k_th_element] = j;
   //          k_th_element++;
   //       }
   //    }
   // }
   // else
   // {
   //    Vector eval = Vector::Zero(n);
   //    Vector lambdas = Vector::Zero(m);
   //    for (int i = 0; i<n; i++)
   //    {
   //       eval(i) = x[i];
   //    }
   //    for(int i = 0; i < m; i++)
   //    {
   //       lambdas(i) = lambda[i];
   //    }
   //    Vector constr_hes = constr_func_ptr_->hessian(eval, constr_params_, lambdas);
   //    Vector cost_hes = cost_func_ptr_->hessian(eval, cost_params_, 0);
   //    // std::cout << "constr hes size : " << constr_hes.rows() << std::endl;
   //    // std::cout << "cost hes size : " << cost_hes.rows() << std::endl;

   //    for(int i = 0; i < n*n; i++)
   //    {
   //       values[i] = obj_factor*cost_hes(i) + constr_hes(i);
   //    }
   // }
      // std::cout << "eval_h" << std::endl;
   return true;
}
// [TNLP_eval_h]

// [TNLP_finalize_solution]
void DOC_NLP_HTconstr::finalize_solution(
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
)
{
   assert(n == nlp_param_.n);
   assert(m == nlp_param_.m);

   cp_est_ = Vector::Zero(nlp_param_.n,1);


   std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
   for( int i = 0; i < n; i++ )
   {
      cp_est_(i) = x[i];
      // std::cout << "x[" << i << "] = " << x[i] << std::endl;
   }
 
   // std::cout << std::endl << std::endl << "Objective value" << std::endl;
   // std::cout << "f(x*) = " << obj_value << std::endl;
}

/**
 * **********************************************************************************************
 * **********************************************************************************************
 * ---------------------------------------- ACCESSORS -------------------------------------------
 * **********************************************************************************************
 * **********************************************************************************************
*/

Cost* DOC_NLP_HTconstr::cost_Func()
{
    /**
     * Returns the cost function pointer
    */
    return cost_func_ptr_;
}

NlpParam& DOC_NLP_HTconstr::nlp_Param()
{
    /**
     * Returns the Non Linear Porblem parameters by reference
    */
    return nlp_param_;
}


MechanicalModel& DOC_NLP_HTconstr::mechanical_Model()
{
    /**
     * Returns the mechanical model by reference
    */
    return mech_model_;
}

MechanicalModel* DOC_NLP_HTconstr::mechanical_Model_Ptr()
{
    /**
     * Returns a pointer to the mechanical model
    */
    return &mech_model_;
}


void DOC_NLP_HTconstr::set_Nlp_Param(NlpParam const& param)
{
    /**
     * Sets the non linear problem parameters
    */
    nlp_param_ = param;
}
