#include <bayesopt/bayesopt.hpp>
#include <NLP_DOC_constr_humantable.h>
#include "IpIpoptApplication.hpp"

#define NB_DOC_ITER 10

class MyOptimization: public bayesopt::ContinuousModel
{
 public:

  NlpParam nlp_param;
  std::vector<double> t_cp;
  std::vector<Eigen::VectorXd> q_cp, dq_cp;
  Ipopt::SmartPtr<DOC_NLP_HTconstr> myNlp;
  std::vector<Eigen::VectorXd> q_ref, dq_ref, ddq_ref;
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app;
  Ipopt::ApplicationReturnStatus status;
  Eigen::VectorXd cp_time;
  int iter_no;
  MyOptimization(bayesopt::Parameters param,
                 int input_dimension):
    bayesopt::ContinuousModel(input_dimension, param),
    iter_no(0)
  {
        std::string cp_path = "../data/spline_traj_s2_trial_1.txt";
        std::string box_params_path = "../data/s2_trial_1_lifting_environment.txt";
        read_control_points(t_cp, q_cp, dq_cp, cp_path);
        cp_time = Eigen::VectorXd::Zero(10,1);
        for(int i = 0; i < t_cp.size(); i++)
        {
            cp_time(i) = t_cp[i];
        }
        q_ref = read_q6_data("../data/q_ref.txt");
        dq_ref = read_q6_data("../data/dq_ref.txt");
        ddq_ref = read_q6_data("../data/ddq_ref.txt");

        nlp_param.n = 60+6+6;
        nlp_param.m = 2 + 5*10 + 10 + 5*10;
        nlp_param.nnz_jac_g = nlp_param.m*nlp_param.n;
        nlp_param.nnz_hes_lag = nlp_param.n*nlp_param.n;

        Eigen::VectorXd x_l(nlp_param.n);
        x_l(0) = q_cp[0](0);
        x_l(1) = q_cp[0](1);
        x_l(2) = q_cp[0](2);
        x_l(3) = q_cp[0](3);
        x_l(4) = q_cp[0](4);
        x_l(5) = q_cp[0](5);
        for (int i = 6; i < x_l.rows()-(6+6); i+=6)
        {
            x_l(i) = 0;
            x_l(i+1) = -0.1745329252;
            x_l(i+2) = -M_PI;
            x_l(i+3) = -1.256637061;
            x_l(i+4) = -4.71238898;
            x_l(i+5) = -0.1745329252;
        }

        Eigen::VectorXd x_u(nlp_param.n);
        x_u(0) = q_cp[0](0);
        x_u(1) = q_cp[0](1);
        x_u(2) = q_cp[0](2);
        x_u(3) = q_cp[0](3);
        x_u(4) = q_cp[0](4);
        x_u(5) = q_cp[0](5);
        for (int i = 6; i < x_u.rows()-(6+6); i+=6)
        {
            x_u(i) = 2.35619449;
            x_u(i+1) = 3.141592654;
            x_u(i+2) = 0.3141592654;
            x_u(i+3) = 0.7853981634;
            x_u(i+4) = 0.2617993878;
            x_u(i+5) = 3.141592654;
        }
        for(int i = 60; i < 66; i++)
        {
            x_l(i) = -0.01;
            x_u(i) = 0.01;
        }
        for(int i = 66; i < 72; i++)
        {
            x_u(i) = 1.5;
            x_l(i) = -1.5;
        }
        nlp_param.x_l = x_l;
        nlp_param.x_u = x_u;

        app = new Ipopt::IpoptApplication();

        // Change some options
        // Note: The following choices are only examples, they might not be
        //       suitable for your optimization problem.
        app->Options()->SetNumericValue("tol", 1e-3);
        //    app->Options()->SetNumericValue("acceptable_tol", 5e-15);
        //    app->Options()->SetIntegerValue("acceptable_iter", 1);
        app->Options()->SetIntegerValue("print_level", 0);
        //    app->Options()->SetStringValue("mehrotra_algorithm", "no");
        //    app->Options()->SetStringValue("mu_strategy", "adaptive");
        app->Options()->SetIntegerValue("max_iter", NB_DOC_ITER);
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
        //    app->Options()->SetStringValue("derivative_test", "first-order");
        //    app->Options()->SetStringValue("derivative_test_print_all", "yes");

        // Initialize the IpoptApplication and process the options
        status = app->Initialize();
        if( status != Ipopt::Solve_Succeeded)
        {
            std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        }
        
  }

  double evaluateSample( const boost::numeric::ublas::vector<double> &query ) 
  {
        
        Eigen::VectorXd weights(Eigen::VectorXd::Zero(7*3,1));
        // for(int i = 0; i < weights.rows(); i++)
        // {
        //     weights(i) = query(i);
        // }

        weights(0) = query(0);
        weights(3) = query(1);
        weights(4) = query(2);
        weights(7) = query(3);
        weights(10) = query(4);
        weights(11) = query(5);
        weights(14) = query(6);
        weights(17) = query(7);
        weights(18) = query(8);
        // std::cout << "weights size out : " << weights.rows() << std::endl;
        // std::cout << "weights out : " << weights.transpose() << std::endl;
        myNlp = new DOC_NLP_HTconstr(nlp_param, weights);

        std::vector<Eigen::VectorXd> q_est, dq_est, ddq_est;
        status = app->OptimizeTNLP(myNlp);

        Eigen::VectorXd state = myNlp->cp_est_;
        Eigen::VectorXd dq_in = dq_ref[0]; Eigen::VectorXd ddq_in = ddq_ref[0];
        Eigen::VectorXd dq_fin = dq_ref[dq_ref.size()-1]; Eigen::VectorXd ddq_fin = ddq_ref[ddq_ref.size()-1];

        Eigen::MatrixXd in_cond(Eigen::MatrixXd::Zero(6,2)),
                        fin_cond(Eigen::MatrixXd::Zero(6,2));
        in_cond.col(0) = dq_in; in_cond.col(1) = ddq_in;
        fin_cond.col(0) = state.segment(60,6); fin_cond.col(1) = state.segment(66,6);

        Eigen::MatrixXd q_cp_eig(Eigen::MatrixXd::Zero(6,10));
        int col_idx(0);
        for(int i = 0; i < state.rows()-(6+6); i+=6)
        {
            q_cp_eig.col(col_idx) = state.segment(i,6);
            col_idx++;
        }

        Eigen::MatrixXd pol_coefs = compute_coefs<double>(5, cp_time, q_cp_eig, in_cond, fin_cond);
        Eigen::MatrixXd q_est_mat(Eigen::MatrixXd::Zero(6,100)),
                        dq_est_mat(Eigen::MatrixXd::Zero(6,100)),
                        ddq_est_mat(Eigen::MatrixXd::Zero(6,100));

        for(int i = 0; i < q_est_mat.rows(); i++)
        {
            q_est_mat.row(i) = eval_polynomial_piecewise<double>(pol_coefs.col(i), cp_time, 2.23, 100, 5);
            dq_est_mat.row(i) = eval_Dpolynomial_piecewise<double>(pol_coefs.col(i), cp_time, 2.23, 100, 5);
            ddq_est_mat.row(i) = eval_DDpolynomial_piecewise<double>(pol_coefs.col(i), cp_time, 2.23, 100, 5);
        }

        for(int i = 0; i < q_est_mat.cols(); i++)
        {
            q_est.push_back(q_est_mat.col(i));
            dq_est.push_back(dq_est_mat.col(i));
            ddq_est.push_back(ddq_est_mat.col(i));
        }

        double mse = 0.0;
        for(int i = 0; i < q_est.size(); i++)
        {
            mse += (q_est[i] - q_ref[i]).transpose()*(q_est[i] - q_ref[i]);
        }
        std::cout << "ITER : " << iter_no << std::endl;
        // std::cout << "Solution : " << weights.transpose() << std::endl;
        std::cout << "Rmse : " << (std::sqrt(mse/(double)q_est.size())*180.0/M_PI) << " deg" << std::endl;
        iter_no++;
        return mse;
  };

  bool checkReachability( const boost::numeric::ublas::vector<double> &query )
  { 
     return true;
  };
};