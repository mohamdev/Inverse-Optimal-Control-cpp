#pragma once
#include <utils.h>
#include <pinocchio/algorithm/rnea-derivatives.hpp>

class MechanicalModel
{
    public:

    Model model_; //Pinocchio model
    Data data_; //Pinocchio data linked to model

    /**** 
     * 
     * Constructors 
     * 
     * ****/
    MechanicalModel();
    MechanicalModel(MechanicalModel const& model) = default;
    MechanicalModel(MechanicalModel && model) = default;
    virtual ~MechanicalModel();

    /**** 
     * 
     * Operators 
     * 
     * ****/
    MechanicalModel& operator=(MechanicalModel const& model) = default;
    MechanicalModel& operator=(MechanicalModel && model) = default;

    /**** 
     * 
     * Accessors 
     * 
     * ****/
    
    /**
     * Returns pinocchio model by reference
    */
    Model& model();

    /**
     * Returns a pointer to pinocchio model
    */
    Model* model_Ptr();

    /**
     * Returns pinocchio Data by reference
    */
    Data& data();

    /**
     * Returns a pointer to pinocchio data object
    */
    Data* data_Ptr();

    /** Inertia in CoM to inertia in Segment*/
    Matrix inertia_in_segment(Matrix const& inertia_in_com,
                              Vector const& com,
                              double const& mass);

    /**
     * This function builds the pinocchio model_ and data_ attributes.
     * It also sets the right sensor IDs in ModelParam structure
    */
    void build_Pinocchio_Model();
    void build_Pinocchio_Model(double shank_l,
                                double thigh_l,
                                double pelv_l,
                                double torso_l,
                                double upa_l,
                                double loa_l,
                                double foot_m,
                                double shank_m,
                                double thigh_m,
                                double pelv_m,
                                double torso_m,
                                double upa_m,
                                double loa_m,
                                double hand_m,
                                double head_m,
                                double foot_in_zz,
                                double shank_in_zz,
                                double thigh_in_zz,
                                double pelv_in_zz,
                                double torso_in_zz,
                                double upa_in_zz,
                                double loa_in_zz,
                                double hand_in_zz,
                                double head_in_zz,
                                Vector foot_com,
                                Vector shank_com,
                                Vector thigh_com,
                                Vector pelv_com,
                                Vector torso_com,
                                Vector upa_com,
                                Vector loa_com,
                                Vector hand_com,
                                Vector head_com);

};

