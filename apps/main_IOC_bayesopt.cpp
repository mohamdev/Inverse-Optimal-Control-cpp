#include <bayesopt_ioc.h>
#include <bayesopt/parameters.hpp>

#define NB_IOC_ITER 30
#define NB_IOC_SEEDS 5

int main()
{
    int nb_bayesopt_vars = 3*3;
    bopt_params param;
    param.l_type = L_MCMC;

    bayesopt::Parameters b_params;

    b_params = initialize_parameters_to_default();

    b_params.n_iterations = NB_IOC_ITER;
    b_params.n_init_samples = NB_IOC_SEEDS;

    b_params.l_type = L_MCMC;
    b_params.sc_type = SC_MAP;

    MyOptimization optimizer(b_params, nb_bayesopt_vars);
    //Define bounds and prepare result.
    boost::numeric::ublas::vector<double> bestPoint(nb_bayesopt_vars);
    boost::numeric::ublas::vector<double> lowerBound(nb_bayesopt_vars);
    boost::numeric::ublas::vector<double> upperBound(nb_bayesopt_vars);

    for(int i = 0; i < nb_bayesopt_vars; i++)
    {
        lowerBound(i) = 1e-5;
        upperBound(i) = 100.0;
    }

    //Set the bounds. This is optional. Default is [0,1]
    //Only required because we are doing continuous optimization
    optimizer.setBoundingBox(lowerBound,upperBound);
    //Collect the result in bestPoint
    optimizer.optimize(bestPoint);

    std::cout << "bayesopt optim ok " << std::endl;

    Eigen::VectorXd weights_sol(Eigen::VectorXd::Zero(7*3,1));
    Eigen::VectorXd weights_write(Eigen::VectorXd::Zero(nb_bayesopt_vars,1));
    for(int i = 0; i < weights_write.rows(); i++)
    {
        weights_write(i) = bestPoint(i);
    }
    weights_sol(0) = bestPoint(0);
    weights_sol(3) = bestPoint(1);
    weights_sol(4) = bestPoint(2);
    weights_sol(7) = bestPoint(3);
    weights_sol(10) = bestPoint(4);
    weights_sol(11) = bestPoint(5);
    weights_sol(14) = bestPoint(6);
    weights_sol(17) = bestPoint(7);
    weights_sol(18) = bestPoint(8);

    std::cout << "Best weights : " << weights_sol.transpose() << std::endl;


    std::vector<double> t_cp;
    Eigen::VectorXd box_pos(Eigen::VectorXd::Zero(3,1));
    Eigen::VectorXd table_center(Eigen::VectorXd::Zero(3,1));
    double box_mass, table_width, table_height;
    std::vector<Eigen::VectorXd> q_cp, dq_cp;
    std::string cp_path = "../data/spline_traj_s2_trial_1.txt";
    std::string box_params_path = "../data/s2_trial_1_lifting_environment.txt";
    read_control_points(t_cp, q_cp, dq_cp, cp_path);
    Eigen::VectorXd cp_time(Eigen::VectorXd::Zero(10,1));
    for(int i = 0; i < t_cp.size(); i++)
    {
        cp_time(i) = t_cp[i];
    }
    Eigen::VectorXd wrist_final_pos(Eigen::VectorXd::Zero(2,1));
    read_final_hand_pos(wrist_final_pos, box_params_path);
    read_box_params(box_mass, box_pos, box_params_path);
    read_table_params(table_width, table_height, table_center, box_params_path);

    std::vector<Eigen::VectorXd> q_ref, dq_ref, ddq_ref;
    q_ref = read_q6_data("../data/q_ref.txt");
    dq_ref = read_q6_data("../data/dq_ref.txt");
    ddq_ref = read_q6_data("../data/ddq_ref.txt");

    Eigen::VectorXd pos_wrist_final(Eigen::VectorXd::Zero(3,1));
    pos_wrist_final.head(2) = wrist_final_pos;

    NlpParam nlp_param;
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

    Ipopt::SmartPtr<DOC_NLP_HTconstr> myNlp = new DOC_NLP_HTconstr(nlp_param, weights_sol);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = new Ipopt::IpoptApplication();

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
    Ipopt::ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Ipopt::Solve_Succeeded)
    {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int) status;
    }

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


    std::string q_write_path = "../output_ioc/q_est_ioc.txt";
    std::string dq_write_path = "../output_ioc/dq_est_ioc.txt";
    std::string ddq_write_path = "../output_ioc/ddq_est_ioc.txt";
    std::string weights_path = "../output_ioc/weights_ioc_solution.txt";
    write_q_data_6dof(q_write_path, q_est);
    write_q_data_6dof(dq_write_path, dq_est);
    write_q_data_6dof(ddq_write_path, ddq_est);
    write_ioc_solution(weights_path, weights_write);

    double mse = 0.0;
    for(int i = 0; i < q_est.size(); i++)
    {
        mse += (q_est[i] - q_ref[i]).transpose()*(q_est[i] - q_ref[i]);
    }
    std::cout << "Solution : " << weights_sol.transpose() << std::endl;
    std::cout << "Rmse : " << (std::sqrt(mse/(double)q_est.size())*180.0)/M_PI << " deg" << std::endl;
    return 0;
}