#pragma once
#include <utils.h>

int cutoff_fact(int const& n,
                int const& k);

template<typename Base>
Base eval_polynomial_atomic(Eigen::Matrix<Base,
                                          Eigen::Dynamic,
                                          1> const & coefs,
                            Base const& t);

template<typename Base>
Base eval_polynomial_atomic(Eigen::Matrix<Base,
                                          Eigen::Dynamic,
                                          1> const & coefs,
                            Base const& t)
{
    Base res(0.0);
    Base tpow(1.0);
    for(int i = 0; i < coefs.rows(); i++)
    {
        res += coefs(i)*tpow;
        tpow *= t;
    }
    return res;
};

template<typename Base>
Base eval_Dpolynomial_atomic(Eigen::Matrix<Base,
                                          Eigen::Dynamic,
                                          1> const & coefs,
                            Base const& t);

template<typename Base>
Base eval_Dpolynomial_atomic(Eigen::Matrix<Base,
                                          Eigen::Dynamic,
                                          1> const & coefs,
                             Base const& t)
{
    Base res(0.0);
    Base tpow(1.0);
    for(int i = 1; i < coefs.rows(); i++)
    {
        res += coefs(i)*tpow*Base(i);
        tpow *= t;
    }
    return res;
};

template<typename Base>
Base eval_DDpolynomial_atomic(Eigen::Matrix<Base,
                                          Eigen::Dynamic,
                                          1> const & coefs,
                            Base const& t);

template<typename Base>
Base eval_DDpolynomial_atomic(Eigen::Matrix<Base,
                                          Eigen::Dynamic,
                                          1> const & coefs,
                             Base const& t)
{
    Base res(0.0);
    Base tpow(1.0);
    for(int i = 2; i < coefs.rows(); i++)
    {
        res += coefs(i)*tpow*Base(i*(i-1));
        tpow *= t;
    }
    return res;
};

template<typename Base>
Base eval_DDDpolynomial_atomic(Eigen::Matrix<Base,
                                          Eigen::Dynamic,
                                          1> const & coefs,
                            Base const& t);

template<typename Base>
Base eval_DDDpolynomial_atomic(Eigen::Matrix<Base,
                                          Eigen::Dynamic,
                                          1> const & coefs,
                             Base const& t)
{
    Base res(0.0);
    Base tpow(1.0);
    for(int i = 3; i < coefs.rows(); i++)
    {
        res += coefs(i)*tpow*Base(i*(i-1)*(i-2));
        tpow *= t;
    }
    return res;
};

template<typename Base>
Eigen::Matrix<Base, Eigen::Dynamic, 1> eval_polynomial_piecewise(Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const & piecewise_coefs,
                                                                 Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const& t_cp,
                                                                 double const& traj_duration,
                                                                 int const& size_traj,
                                                                 int const& deg);

template<typename Base>
Eigen::Matrix<Base, Eigen::Dynamic, 1> eval_polynomial_piecewise(Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const & piecewise_coefs,
                                                                 Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const& t_cp,
                                                                 double const& traj_duration,
                                                                 int const& size_traj,
                                                                 int const& deg)
{
    int nb_splines = int(piecewise_coefs.rows()/(deg+1));
    // std::cout << "nb splines : " << nb_splines << std::endl;
    double nb_points_per_spline = (double)size_traj/(double)nb_splines;

    Base dt = Base(traj_duration/(size_traj-1));
    Eigen::Matrix<Base, Eigen::Dynamic, 1> traj_val(Eigen::Matrix<Base, Eigen::Dynamic, 1>::Zero(size_traj,1));
    int coefs_idx(0);


    for(int i = 0; i < size_traj; i++)
    {
        // std::cout << " i : " << i << std::endl;
        int curr_spline_order = std::floor((double)i/nb_points_per_spline);
        coefs_idx = (deg+1)*curr_spline_order;
        // std::cout << "coefs_idx : " << coefs_idx << std::endl;
        Base t_curr = t_cp(0) + i*dt;
        // std::cout << "t_curr : " << t_curr << std::endl; 
        Eigen::Matrix<Base, Eigen::Dynamic, 1> coefs = piecewise_coefs.segment(coefs_idx, deg+1);
        traj_val(i) = eval_polynomial_atomic<Base>(coefs, t_curr);
    }

    return traj_val;
};


template<typename Base>
Eigen::Matrix<Base, Eigen::Dynamic, 1> eval_Dpolynomial_piecewise(Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const & piecewise_coefs,
                                                                 Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const& t_cp,
                                                                 double const& traj_duration,
                                                                 int const& size_traj,
                                                                 int const& deg);

template<typename Base>
Eigen::Matrix<Base, Eigen::Dynamic, 1> eval_Dpolynomial_piecewise(Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const & piecewise_coefs,
                                                                 Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const& t_cp,
                                                                 double const& traj_duration,
                                                                 int const& size_traj,
                                                                 int const& deg)
{
    int nb_splines = int(piecewise_coefs.rows()/(deg+1));
    // std::cout << "nb splines : " << nb_splines << std::endl;
    double nb_points_per_spline = (double)size_traj/(double)nb_splines;

    Base dt = Base(traj_duration/(size_traj-1));
    Eigen::Matrix<Base, Eigen::Dynamic, 1> traj_val(Eigen::Matrix<Base, Eigen::Dynamic, 1>::Zero(size_traj,1));
    int coefs_idx(0);


    for(int i = 0; i < size_traj; i++)
    {
        // std::cout << " i : " << i << std::endl;
        int curr_spline_order = std::floor((double)i/nb_points_per_spline);
        coefs_idx = (deg+1)*curr_spline_order;
        // std::cout << "coefs_idx : " << coefs_idx << std::endl;
        Base t_curr = t_cp(0) + i*dt;
        // std::cout << "t_curr : " << t_curr << std::endl; 
        Eigen::Matrix<Base, Eigen::Dynamic, 1> coefs = piecewise_coefs.segment(coefs_idx, deg+1);
        traj_val(i) = eval_Dpolynomial_atomic<Base>(coefs, t_curr);
    }

    return traj_val;
};

template<typename Base>
Eigen::Matrix<Base, Eigen::Dynamic, 1> eval_DDpolynomial_piecewise(Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const & piecewise_coefs,
                                                                 Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const& t_cp,
                                                                 double const& traj_duration,
                                                                 int const& size_traj,
                                                                 int const& deg);

template<typename Base>
Eigen::Matrix<Base, Eigen::Dynamic, 1> eval_DDpolynomial_piecewise(Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const & piecewise_coefs,
                                                                 Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const& t_cp,
                                                                 double const& traj_duration,
                                                                 int const& size_traj,
                                                                 int const& deg)
{
    int nb_splines = int(piecewise_coefs.rows()/(deg+1));
    // std::cout << "nb splines : " << nb_splines << std::endl;
    double nb_points_per_spline = (double)size_traj/(double)nb_splines;

    Base dt = Base(traj_duration/(size_traj-1));
    Eigen::Matrix<Base, Eigen::Dynamic, 1> traj_val(Eigen::Matrix<Base, Eigen::Dynamic, 1>::Zero(size_traj,1));
    int coefs_idx(0);


    for(int i = 0; i < size_traj; i++)
    {
        // std::cout << " i : " << i << std::endl;
        int curr_spline_order = std::floor((double)i/nb_points_per_spline);
        coefs_idx = (deg+1)*curr_spline_order;
        // std::cout << "coefs_idx : " << coefs_idx << std::endl;
        Base t_curr = t_cp(0) + i*dt;
        // std::cout << "t_curr : " << t_curr << std::endl; 
        Eigen::Matrix<Base, Eigen::Dynamic, 1> coefs = piecewise_coefs.segment(coefs_idx, deg+1);
        traj_val(i) = eval_DDpolynomial_atomic<Base>(coefs, t_curr);
    }

    return traj_val;
};

template<typename Base>
Eigen::Matrix<Base, Eigen::Dynamic, 1> eval_DDDpolynomial_piecewise(Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const & piecewise_coefs,
                                                                 Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const& t_cp,
                                                                 double const& traj_duration,
                                                                 int const& size_traj,
                                                                 int const& deg);

template<typename Base>
Eigen::Matrix<Base, Eigen::Dynamic, 1> eval_DDDpolynomial_piecewise(Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const & piecewise_coefs,
                                                                 Eigen::Matrix<Base,
                                                                               Eigen::Dynamic,
                                                                               1> const& t_cp,
                                                                 double const& traj_duration,
                                                                 int const& size_traj,
                                                                 int const& deg)
{
    int nb_splines = int(piecewise_coefs.rows()/(deg+1));
    // std::cout << "nb splines : " << nb_splines << std::endl;
    double nb_points_per_spline = (double)size_traj/(double)nb_splines;

    Base dt = Base(traj_duration/(size_traj-1));
    Eigen::Matrix<Base, Eigen::Dynamic, 1> traj_val(Eigen::Matrix<Base, Eigen::Dynamic, 1>::Zero(size_traj,1));
    int coefs_idx(0);


    for(int i = 0; i < size_traj; i++)
    {
        // std::cout << " i : " << i << std::endl;
        int curr_spline_order = std::floor((double)i/nb_points_per_spline);
        coefs_idx = (deg+1)*curr_spline_order;
        // std::cout << "coefs_idx : " << coefs_idx << std::endl;
        Base t_curr = t_cp(0) + i*dt;
        // std::cout << "t_curr : " << t_curr << std::endl; 
        Eigen::Matrix<Base, Eigen::Dynamic, 1> coefs = piecewise_coefs.segment(coefs_idx, deg+1);
        traj_val(i) = eval_DDDpolynomial_atomic<Base>(coefs, t_curr);
    }
    return traj_val;
};

template<typename Base>
Eigen::Matrix<Base,Eigen::Dynamic,Eigen::Dynamic> compute_coefs(int deg,
                                                                Eigen::Matrix<Base,
                                                                                Eigen::Dynamic,
                                                                                1> const & knotTimes,
                                                                Eigen::Matrix<Base,
                                                                                Eigen::Dynamic,
                                                                                Eigen::Dynamic> const & knotValues,
                                                                Eigen::Matrix<Base,
                                                                            Eigen::Dynamic,
                                                                            Eigen::Dynamic> const & in_cond,
                                                                Eigen::Matrix<Base,
                                                                             Eigen::Dynamic,
                                                                             Eigen::Dynamic> const & fin_cond);

template<typename Base>
Eigen::Matrix<Base,Eigen::Dynamic,Eigen::Dynamic> compute_coefs(int deg,
                                                                Eigen::Matrix<Base,
                                                                                Eigen::Dynamic,
                                                                                1> const & knotTimes,
                                                                Eigen::Matrix<Base,
                                                                                Eigen::Dynamic,
                                                                                Eigen::Dynamic> const & knotValues,
                                                                Eigen::Matrix<Base,
                                                                            Eigen::Dynamic,
                                                                            Eigen::Dynamic> const & in_cond,
                                                                Eigen::Matrix<Base,
                                                                            Eigen::Dynamic,
                                                                            Eigen::Dynamic> const & fin_cond)
{
    /**
     * This function computes the coefficients polynomials of degrees deg
     * @param deg degrees of the polnomials between the control points
     * @param n_bnd_cnd number of boundary conditions
     * @param knotTimes vector containing the time values of the control points
     * @param knotValues matrix containing the values of the control poitns [cp_0, ..., cp_T]
    */

    /*Number of polynomials*/
    int n = knotTimes.rows() - 1;

    /*Number of boundary conditions*/
    int n_bnd_cnd = deg - 1;

    /* A */
    Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic> A = Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic>::
                                                            Zero((deg+1)*n, (deg+1)*n);

    /* b */
    Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic> b = Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic>::
                                                            Zero((deg+1)*n, knotValues.rows());

    /* For each polynomial up to the next to last one add interior and exterior 
       knot-point equalities*/
    for(int i = 0; i < n-1; i++)
    {
        for(int j = 0; j < (deg+1); j++)
        {
            A((deg+1)*i, (deg+1)*i + j) = CppAD::pow<Base>(knotTimes(i), j);
            A((deg+1)*i + 1, (deg+1)*i + j) = CppAD::pow<Base>(knotTimes(i+1), j);
            A((deg+1)*i + 1, (deg+1)*(i+1) + j) = -CppAD::pow<Base>(knotTimes(i+1), j);
        }
        b.row((deg+1)*i) = knotValues.col(i).transpose();
    }

    for(int j = 0; j < (deg+1); j++)
    {
        A((deg+1)*(n-1), (deg+1)*(n-1) + j) = CppAD::pow<Base>(knotTimes(n-1), j);
        A((deg+1)*(n-1) + 1, (deg+1)*(n-1) + j) = CppAD::pow<Base>(knotTimes(n), j);
    }
    b.row((deg+1)*(n-1)) = knotValues.col(n-1).transpose();
    b.row((deg+1)*(n-1) + 1) = knotValues.col(n).transpose();

    for(int i = 0; i < n-1; i++)
    {
        for(int d = 1; d < deg; d++)
        {
            for(int j = d; j < deg+1; j++)
            {
                A((deg+1)*i + 1 + d, (deg+1)*i + j) = cutoff_fact(j, d)
                                                      *CppAD::pow<Base>(knotTimes(i+1), j-d);
                A((deg+1)*i + 1 + d, (deg+1)*i + deg + 1 + j) = -cutoff_fact(j, d)
                                                                *CppAD::pow<Base>(knotTimes(i+1), j-d);
            }
        }
    }

    for(int num_con = 1; num_con <= n_bnd_cnd/2; num_con++)
    {
        int d = num_con;
        int num_knot = 0;
        int num_pol = num_knot;
        Base x_c = knotTimes(num_knot);

        for(int j = d; j < deg+1; j++)
        {
            A((deg+1)*(n-1) + 1 + num_con, (deg+1)*(num_pol) + j) = cutoff_fact(j, d)
                                                                      *CppAD::pow<Base>(x_c, j-d);
        }
        // std::cout << "NIK OMOK " << n << " " << deg << std::endl;
        b.row((deg+1)*(n-1) + 1 + num_con) = in_cond.col(num_con-1).transpose();
        // std::cout << b << std::endl;
    }
    

    for(int num_con = 1; num_con <= n_bnd_cnd/2; num_con++)
    {
        int d = num_con;
        int num_knot = n;
        int num_pol = num_knot-1;
        Base x_c = knotTimes(num_knot);

        for(int j = d; j < deg+1; j++)
        {
            A((deg+1)*(n-1) + 1 + num_con + n_bnd_cnd/2, (deg+1)*(num_pol) + j) = cutoff_fact(j, d)
                                                                                    *CppAD::pow<Base>(x_c, j-d);
        }

        b.row((deg+1)*(n-1) + 1 + num_con + n_bnd_cnd/2) = fin_cond.col(num_con-1).transpose();
    }

    Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic> invA = A.inverse();

    return invA*b;
};