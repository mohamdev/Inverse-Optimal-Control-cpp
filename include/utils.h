#pragma once
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <netdb.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <stdlib.h>
#include <random>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <limits>
#include <fstream>
#include <time.h>
#include <ctime>
#include <math.h>
#include <mutex>
// #include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/codegen/cppadcg.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/math/quaternion.hpp"
#include <pinocchio/algorithm/frames.hpp>
#include "pinocchio/autodiff/cppad.hpp"
#include <cppad/cppad.hpp> // the CppAD package
// #include "cppad/cg/cppadcg.hpp"

namespace pin = pinocchio;
using namespace CppAD::cg;
using namespace CppAD;

typedef double Scalar;
typedef pin::Model Model;
typedef pin::Model::Data Data;
typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> Vector;
typedef AD<Scalar> ADScalar;
typedef pin::ModelTpl<ADScalar> ADModel;
typedef ADModel::Data ADData;
typedef Eigen::Matrix<ADScalar,Eigen::Dynamic,1> ADVector;
typedef Eigen::Matrix<ADScalar,Eigen::Dynamic,Eigen::Dynamic> ADMatrix;
typedef ADModel::ConfigVectorType ADConfigVectorType;
typedef ADModel::TangentVectorType ADTangentVectorType;

typedef CppAD::cg::CG<Scalar> CGScalar;
typedef CppAD::AD<CGScalar> ADCGScalar;
typedef pinocchio::ModelTpl<ADCGScalar> CGModel;
typedef CGModel::Data CGData;
typedef CGModel::ConfigVectorType CGConfigVectorType;
typedef CGModel::TangentVectorType CGTangentVectorType;
typedef Eigen::Matrix<ADCGScalar,Eigen::Dynamic,Eigen::Dynamic> CGMatrix;
typedef Eigen::Matrix<ADCGScalar,Eigen::Dynamic,1> CGVector;

/**
 * This function is used to convert a std::vector of Eigen column vectors
 * into an Eigen matrix
 * @param vector_mat The vector of Vectors to convert
 * @see Vector
*/
Matrix vector_To_Eigen(std::vector<Vector> & vector_mat);


std::vector<std::vector<double>> eigen_To_Vector(Matrix const& eigen_mat);


template<typename Base>
Eigen::Matrix<Base,Eigen::Dynamic,1> rot2quat(Eigen::Matrix<Base,Eigen::Dynamic,Eigen::Dynamic> const & m);

template<typename Base>
Eigen::Matrix<Base,Eigen::Dynamic,1> rot2quat(Eigen::Matrix<Base,Eigen::Dynamic,Eigen::Dynamic> const & m){
    /**
     * This function is used to convert a rotation matrix into a quaternion using a given
     * base (generaly Scalar which is based on double, or ADScalar ..)
     * @param m 3*3 Rotation matrix
    */
    Base qw = CppAD::sqrt(1.0+m(0,0)+m(1,1)+m(2,2))/2.0;
    Base qx = (m(2,1) - m(1,2))/( 4.0 *qw);
    Base qy = (m(0,2) - m(2,0))/( 4.0 *qw);
    Base qz = (m(1,0) - m(0,1))/( 4.0 *qw);
    Eigen::Matrix<Base,Eigen::Dynamic,1> quat(4,1);
    // quat(0) = qw;
    // quat(1) = qx;
    // quat(2) = qy;
    // quat(3) = qz;
    quat(3) = qw;
    quat(0) = qx;
    quat(1) = qy;
    quat(2) = qz;
    quat.normalize();
    return quat;
};

void write_q_data_6dof(std::string const& path, std::vector<CGVector> const& q_data);

template<typename Base>
Eigen::Matrix<Base,Eigen::Dynamic,1> rot2quat_new(Eigen::Matrix<Base,Eigen::Dynamic,Eigen::Dynamic> const & m){
    /**
     * This function is used to convert a rotation matrix into a quaternion using a given
     * base (generaly Scalar which is based on double, or ADScalar ..)
     * @param m 3*3 Rotation matrix
    */
    Base qw = CppAD::sqrt(1.0+m(0,0)+m(1,1)+m(2,2))/2.0;
    Base qx = (m(2,1) - m(1,2))/( 4.0 *qw);
    Base qy = (m(0,2) - m(2,0))/( 4.0 *qw);
    Base qz = (m(1,0) - m(0,1))/( 4.0 *qw);
    Eigen::Matrix<Base,Eigen::Dynamic,1> quat(4,1);
    // quat(0) = qw;
    // quat(1) = qx;
    // quat(2) = qy;
    // quat(3) = qz;
    quat(0) = qw;
    quat(1) = qx;
    quat(2) = qy;
    quat(3) = qz;
    quat.normalize();
    return quat;
};

template<typename Base>
Eigen::Matrix<Base,Eigen::Dynamic,Eigen::Dynamic> quat2rot_new(Eigen::Matrix<Base,Eigen::Dynamic,1> const& quat)
{
    Eigen::Quaterniond q;
    q.x() = quat(1);
    q.y() = quat(2);
    q.z() = quat(3);
    q.w() = quat(0); 

    return q.toRotationMatrix();
} 


void write_anat_pos_orientation_data_2vimus(std::string const& path,
                                std::vector<Vector> const& quats_trunk,
                                std::vector<Vector> const& quats_hand,
                                std::vector<Vector> const& pos_trunk,
                                std::vector<Vector> const& pos_hand);

template<typename Base>
Eigen::Matrix<Base,Eigen::Dynamic,1> conjugate(Eigen::Matrix<Base,Eigen::Dynamic, 1> const & q)
{
    Eigen::Matrix<Base,Eigen::Dynamic,1> q_conj = q;
    q_conj(0) = -q_conj(0);
    q_conj(1) = -q_conj(1);
    q_conj(2) = -q_conj(2);

    q_conj.normalize();
    return q_conj;
}

template<typename Base>
Eigen::Matrix<Base,Eigen::Dynamic,1> inverse_quat(Eigen::Matrix<Base,Eigen::Dynamic, 1> const & q)
{
    Eigen::Matrix<Base,Eigen::Dynamic,1> q_conj = conjugate<Base>(q);
    Eigen::Matrix<Base,Eigen::Dynamic,1> q_inv = q_conj/(q.transpose()*q);
    q_inv.normalize();
    return q_inv;
}

template<typename Base>
Eigen::Matrix<Base,Eigen::Dynamic,1> multiply_quat(Eigen::Matrix<Base,Eigen::Dynamic, 1> const & q1,
                                                  Eigen::Matrix<Base,Eigen::Dynamic, 1> const & q2)
{
    Eigen::Matrix<Base,Eigen::Dynamic,1> q(Eigen::Matrix<Base,Eigen::Dynamic,1>::Zero(4,1));
    q(3) = q1(3)*q2(3) - q1(0)*q2(0) - q1(1)*q2(1) - q1(2)*q2(2);
    q(0) = q1(0)*q2(3) + q1(3)*q2(0) + q1(1)*q2(2) - q1(2)*q2(1);
    q(1) = q1(3)*q2(1) - q1(0)*q2(2) + q1(1)*q2(3) + q1(2)*q2(0);
    q(2) = q1(3)*q2(2) + q1(0)*q2(1) - q1(1)*q2(0) + q1(2)*q2(3);
    q.normalize();
    return q;
}

template<typename Base>
Eigen::Matrix<Base,Eigen::Dynamic,Eigen::Dynamic> skew_quat_4_new(Eigen::Matrix<Base,Eigen::Dynamic, 1> const & q)
{
    Eigen::Matrix<Base,Eigen::Dynamic,Eigen::Dynamic> skew(Eigen::Matrix<Base,Eigen::Dynamic,Eigen::Dynamic>::Zero(4,4));
    skew(0,0) = q(0);
    skew(1,1) = q(0);
    skew(2,2) = q(0);
    skew(3,3) = q(0);
    skew(0,1) = -q(1); 
    skew(0,2) = -q(2);
    skew(0,3) = -q(3);
    skew(1,0) = q(1);
    skew(1,2) = -q(3);
    skew(1,3) = q(2);
    skew(2,0) = q(2);
    skew(2,1) = q(3);
    skew(2,3) = -q(1);
    skew(3,0) = q(3);
    skew(3,1) = -q(2);
    skew(3,2) = q(1);
    return skew;
}

template<typename Base>
Eigen::Matrix<Base,Eigen::Dynamic,1> multiply_quat_new(Eigen::Matrix<Base,Eigen::Dynamic, 1> const & q1,
                                                       Eigen::Matrix<Base,Eigen::Dynamic, 1> const & q2)
{
    Eigen::Matrix<Base,Eigen::Dynamic,1> q = skew_quat_4_new<Base>(q1) * q2;
    q.normalize();
    return q;
}

void write_meas_data(std::string const& path, std::vector<Vector> const& meas_data);
void write_meas_data_4vimus(std::string const& path, std::vector<Vector> const& meas_data);

void write_q_data(std::string const& path, std::vector<Vector> const& q_data);
void write_q_data_6dof(std::string const& path, std::vector<Vector> const& q_data);

void write_calib_data(std::string const& path_vimu_local,
                      std::string const& path_anat_local,
                      std::map<std::string, Matrix> calib_data_map);

void write_anat_orientation_data(std::string const& path,
                                std::vector<Vector> const& quats_trunk,
                                std::vector<Vector> const& quats_upperarm,
                                std::vector<Vector> const& quats_lowerarm,
                                std::vector<Vector> const& quats_hand);

void write_arm_3D_pos_data(std::string const& path,
                            std::vector<Vector> const& pos_shoulder,
                            std::vector<Vector> const& pos_elbow,
                            std::vector<Vector> const& pos_wrist);

void write_2vimus_data(std::string const& path,
                        std::vector<Vector> const& quats_trunk,
                        std::vector<Vector> const& quats_hand,
                        std::vector<Vector> const& pos_trunk,
                        std::vector<Vector> const& pos_hand,
                        std::vector<Vector> const& gyr_trunk,
                        std::vector<Vector> const& gyr_hand,
                        std::vector<Vector> const& acc_trunk,
                        std::vector<Vector> const& acc_hand);

void write_calib_data_2vimus(std::string const& path_vimu_local,
                             std::string const& path_anat_local,
                             std::map<std::string, Matrix> calib_data_map);

void write_anat_orientation_data_2vimus(std::string const& path,
                                std::vector<Vector> const& quats_trunk,
                                std::vector<Vector> const& quats_hand);



std::vector<Vector> read_q_data(std::string const& path);
std::vector<Vector> read_q7_data(std::string const& path);
std::vector<Vector> read_q6_data(std::string const& path);

Vector compute_xyz_rmse(std::vector<Vector> est_xyz,
                        std::vector<Vector> ref_xyz);

Vector compute_xyz_covariance(std::vector<Vector> est_xyz,
                              std::vector<Vector> ref_xyz);

Vector compute_xyz_PC(std::vector<Vector> est_xyz,
                      std::vector<Vector> ref_xyz);

double calculatePC(Vector v1, Vector v2);

enum DataType{
    EST, 
    REF, 
    MEAS
};

// Calculates rotation matrix to euler angles
// The result is the same as MATLAB except the order
// of the euler angles ( x and z are swapped ).
Vector rotm_2_euler(Matrix const&R);

void challis();

struct EkfParam
{
    /**
     * Data structure containing the number of state variables, the number of 
     * degrees of freedom, the number of measurement variables and the sampling time
    */
    unsigned int nb_states; //Number of state variables
    unsigned int nb_dof; //Number of degrees of freedom
    unsigned int nb_meas_variables; //Number of measurement variables
    unsigned int nb_bias;
    double dt; //Sampling time
    bool cond_exp;
    double r_acc_gain;
    double r_gyr_gain;
    double r_pos_gain;
    double r_rot_gain;
    double q_gain;
    double dq_gain;
    double ddq_gain;
    bool custom_meas_gain;
    bool custom_model_gain;

    EkfParam() = default;
    EkfParam(EkfParam const & param) = default;
    EkfParam(EkfParam && param) = default;
    EkfParam& operator=(EkfParam const & param) = default;
    EkfParam& operator=(EkfParam && param) = default;
};

struct NlpParam
{
    /**
     * Data structore containing all the parameters needed by the Ipopt::TNLP class
     * to solve a non linear problem (NLP)
    */
    int n; //Number of variables x to estimate
    int m; //Number of constraints g(x)
    int nb_meas;
    int nnz_jac_g; //Number of non zeros of the Jacobian of g(x)
    int nnz_hes_lag; //Number of non zeros of the Hessian of the lagrangian of f
    Vector x_l; //Vector of lower bounds of x variables
    Vector x_u; //Vector of upper bounds of x variables
    Vector g_l; //Vector of lower bounds of constraints
    Vector g_u; //Vector of upper bounds of constraints

    NlpParam() = default;
    NlpParam(NlpParam const & param) = default;
    NlpParam(NlpParam && param) = default;
    NlpParam& operator=(NlpParam const & param) = default;
    NlpParam& operator=(NlpParam && param) = default;
};

struct ModelParam
{
    /**
     * Data structure containing the number of state variables, the number of 
     * degrees of freedom, the number of measurement variables and the sampling time
    */
    std::string urdf_path; //Path to the urdf model
    unsigned int nb_states; //Number of state variables
    unsigned int nb_dof; //Number of degrees of freedom
    unsigned int nb_meas_variables; //Number of measurement variables
    double dt; //Sampling time
    std::map<std::string, int> sensors_id_; //Map containing the ID of sensors

    ModelParam() = default;
    ModelParam(ModelParam const & param) = default;
    ModelParam(ModelParam && param) = default;
    ModelParam& operator=(ModelParam const & param) = default;
    ModelParam& operator=(ModelParam && param) = default;
};

struct ModelStates{
    /**
     * Data struct containing the joint variables vectors and the measurement
     * vectors at the current timestamp and previous timestamp
    */
    ModelParam param;
    DataType type_data; //Type of data
    Vector q; //Angles
    Vector dq; //Articular velocities
    Vector ddq; //Articular accelerations
    Vector state; //State vector containing q, dq and ddq and eventually other variables
    Vector meas; //Measurements vector
    //-------- A VIRER ---////
    Vector prev_q; //Previous angles
    Vector prev_dq; //Previous articular velocities
    Vector prev_ddq; //previous articular accelerations
    Vector prev_state; //previous state containing q, dq and ddq and eventually other variables
    Vector prev_meas; //Previous measurement vector

    ModelStates() = default;
    ModelStates(ModelStates const & j) = default;
    ModelStates(ModelStates && j) = default;
    ModelStates& operator=(ModelStates const & j) = default;
    ModelStates& operator=(ModelStates && j) = default;

    ModelStates(ModelParam const& model_param,
                DataType const& type)
    :
    type_data(type),
    param(model_param),
    q(Vector::Zero(model_param.nb_dof)),
    dq(Vector::Zero(model_param.nb_dof)),
    ddq(Vector::Zero(model_param.nb_dof)),
    prev_q(Vector::Zero(model_param.nb_dof)),
    prev_dq(Vector::Zero(model_param.nb_dof)),
    prev_ddq(Vector::Zero(model_param.nb_dof)),
    meas(Vector::Zero(model_param.nb_meas_variables)),
    prev_meas(Vector::Zero(model_param.nb_meas_variables))
    {};
    void initialize(EkfParam const & param)
    {
        q = Matrix::Zero(param.nb_dof, 1); //Angles
        dq = Matrix::Zero(param.nb_dof, 1); //Articular velocities
        ddq = Matrix::Zero(param.nb_dof, 1); //Articular accelerations
        state = Matrix::Zero(param.nb_states, 1); //State vector containing q, dq and ddq and eventually other variables
        prev_q = Matrix::Zero(param.nb_dof, 1); //Previous angles
        prev_dq = Matrix::Zero(param.nb_dof, 1); //Previous articular velocities
        prev_ddq = Matrix::Zero(param.nb_dof, 1); //previous articular accelerations
        prev_state = Matrix::Zero(param.nb_states, 1); //previous state containing q, dq and ddq and eventually other variables
        meas = Matrix::Zero(param.nb_meas_variables, 1); //Measurements vector
        prev_meas = Matrix::Zero(param.nb_meas_variables, 1); //Predicted measurement vector
    };
};


struct Trajectories
{
    /**
     * Data struct containing the trajectories of joint variables, and measurements
    */
    DataType type_data; //Type of data
    std::vector<Vector> q; //Angles trajectories
    std::vector<Vector> dq; //Articular velocities trajectories
    std::vector<Vector> ddq; //Articular accelerations trajectories
    std::vector<Vector> state; //State vector trajectories
    std::vector<Vector> meas; //Measurement trajectories
    
    Trajectories() = default;
    Trajectories(DataType const& type)
    : type_data(type){};
    Trajectories(Trajectories const & t) = default;
    Trajectories(Trajectories && t) = default;
    Trajectories& operator=(Trajectories const & t) = default;
    Trajectories& operator=(Trajectories && t) = default;
    void clear()
    {
        q.clear();
        dq.clear();
        ddq.clear();
        state.clear();
        meas.clear();
    }
};

struct EkfMatrices
{
    /**
     * Data struct containing the Matrices needed by the EKF
    */
    Matrix P; //Estimated state covariance matrix
    Matrix inv_P; //Inverse of P
    Matrix A; //Evolution process matrix X = A*X;
    Matrix ID_state; //Identity with dimension nb_states*nb_states
    Matrix H; //Jacobian matrix of measurement function
    Matrix S; //Covariance matrix of the measurement function
    Matrix inv_S; //Inverse of S measurement covariance matrix
    Matrix ID_meas; //Identity with dimension nb_meas_var*nb_meas_var
    Matrix Q;  //Covariance of process evolution noise
    Matrix R; //Covariance of measurement noise
    Matrix K; //Kalman gain matrix K

    EkfMatrices() = default;
    EkfMatrices(EkfMatrices const & mat) = default;
    EkfMatrices(EkfMatrices && mat) = default;
    EkfMatrices& operator=(EkfMatrices const & mat) = default;
    EkfMatrices& operator=(EkfMatrices && mat) = default;

    void initialize(EkfParam const & param)
    {
        P = Matrix::Identity(param.nb_states, param.nb_states) ; 
        inv_P = Matrix::Identity(param.nb_states, param.nb_states);
        A = Matrix::Zero(param.nb_states, param.nb_states);
        ID_state = Matrix::Identity(param.nb_states, param.nb_states);
        H = Matrix::Zero(param.nb_meas_variables, param.nb_states);
        S = Matrix::Zero(param.nb_meas_variables, param.nb_meas_variables);
        inv_S = Matrix::Zero(param.nb_meas_variables, param.nb_meas_variables);
        ID_meas = Matrix::Identity(param.nb_meas_variables, param.nb_meas_variables);
        Q = Matrix::Zero(param.nb_states, param.nb_states);
        R = Matrix::Zero(param.nb_meas_variables, param.nb_meas_variables);
        K = Matrix::Zero(param.nb_states, param.nb_meas_variables); 
    };

};

void write_ioc_solution(std::string const& path, Vector const& weights);

void read_control_points(std::vector<double> & t_cp,
                         std::vector<Vector> & q_cp,
                         std::vector<Vector> & dq_cp,
                         std::string const& cp_path);

void read_box_params(double & mass,
                     Vector & pos,
                     std::string const& box_params_path);

void read_box_params(double & mass,
                     Vector & pos,
                     double & box_width,
                     std::string const& box_params_path);

void read_final_hand_pos(Vector & pos,
                         std::string const& box_params_path);

void read_table_params(double & tw,
                       double & th,
                       Vector & pos,
                       std::string const& table_params_path);
