#include <mechanical_model.h>

MechanicalModel::MechanicalModel() : model_(), data_()
{
}

MechanicalModel::~MechanicalModel(){};

/**
 * **********************************************************************************************
 * **********************************************************************************************
 * ---------------------------------------- ACCESSORS -------------------------------------------
 * **********************************************************************************************
 * **********************************************************************************************
*/

Model& MechanicalModel::model()
{
/**
 * Returns pinocchio model by reference
*/
return model_;
}

Model* MechanicalModel::model_Ptr()
{
/**
 * Returns a pointer to pinocchio model
*/
return &model_;
}


Data& MechanicalModel::data()
{
/**
 * Returns pinocchio Data by reference
*/
return data_;
}


Data* MechanicalModel::data_Ptr()
{
/**
 * Returns a pointer to pinocchio data
*/
return &data_;
}

Matrix MechanicalModel::inertia_in_segment(Matrix const& inertia_in_com,
                                           Vector const& com,
                                           double const& mass)
{

    Matrix i_in_segment(Matrix::Zero(3,3)); // Inertia expressed in segment base frame
    Matrix id_3D(Matrix::Identity(3, 3)); // 3D identity
    i_in_segment =
        inertia_in_com +
        mass *
            ((com.transpose() * com) * id_3D -
             com * com.transpose());
    return i_in_segment;
}

void MechanicalModel::build_Pinocchio_Model()
{
        double shank_l = 0.4281073446;
        double thigh_l = 0.4271186441;
        double pelv_l = 0.2412429379;
        double torso_l = 0.3302259887;
        double upa_l = 0.2738700565;
        double loa_l = 0.2798022599;
        double foot_m = 1.60864169;
        double shank_m = 6.434566761;
        double thigh_m = 16.48857732;
        double pelv_m = 11.46157204;
        double torso_m = 20.37612808;
        double upa_m = 3.21728338;
        double loa_m = 2.278909061;
        double hand_m = 0.8043208451;
        double head_m = 4.490791385;
        double foot_in_zz = 0.01585963353;
        double shank_in_zz =0.505688988;
        double thigh_in_zz =1.256514753;
        double pelv_in_zz = 0.1416091287;
        double torso_in_zz =0.7279829267;
        double upa_in_zz = 0.07785874439;
        double loa_in_zz = 0.04504209355;
        double hand_in_zz = 0.005153979839;
        double head_in_zz = 0.1274138803;
        Vector foot_com(Vector::Zero(3,1)),
        shank_com(Vector::Zero(3,1)),
        thigh_com(Vector::Zero(3,1)),
        pelv_com(Vector::Zero(3,1)),
        torso_com(Vector::Zero(3,1)),
        upa_com(Vector::Zero(3,1)),
        loa_com(Vector::Zero(3,1)),
        hand_com(Vector::Zero(3,1)),
        head_com(Vector::Zero(3,1));
        foot_com(0) = -0.02734844633;
        shank_com(0) = 0.2525833333;
        thigh_com(0) = 0.2438847458;
        pelv_com(0) = 0.0865956512;
        torso_com(0) = 0.146950565;
        upa_com(0) =  0.1320053672;
        loa_com(0) =  0.1166775424;
        hand_com(0) = 0.06636158192;
        head_com(0) = 0.1467745763;
        foot_com(1) =  -0.06898954802;
        shank_com(1) = 0.02054915254;
        thigh_com(1) = 0.01751186441;
        pelv_com(1) =  -0.004303403046;
        torso_com(1)=  0.0;
        upa_com(1) =   0.004929661017;
        loa_com(1) =   -0.003637429379;
        hand_com(1) =  0.006485875706;
        head_com(1) =  -0.005497175141;

        // double shank_l = 1.0;
        // double thigh_l = 1.0;
        // double pelv_l = 1.0;
        // double torso_l = 1.0;
        // double upa_l = 1.0;
        // double loa_l = 1.0;
        // double foot_m = 1.0;
        // double shank_m = 1.0;
        // double thigh_m = 1.0;
        // double pelv_m = 1.0;
        // double torso_m = 1.0;
        // double upa_m = 1.0;
        // double loa_m = 1.0;
        // double hand_m = 1e-15;
        // double head_m = 1e-15;
        // double foot_in_zz = 1.0/3.0;
        // double shank_in_zz = 1.0/3.0;
        // double thigh_in_zz = 1.0/3.0;
        // double pelv_in_zz = 1.0/3.0;
        // double torso_in_zz = 1.0/3.0;
        // double upa_in_zz = 1.0/3.0;
        // double loa_in_zz = 1.0/3.0;
        // double hand_in_zz = 1e-15;
        // double head_in_zz = 1e-15;
        // Vector foot_com(Vector::Zero(3,1)),
        // shank_com(Vector::Zero(3,1)),
        // thigh_com(Vector::Zero(3,1)),
        // pelv_com(Vector::Zero(3,1)),
        // torso_com(Vector::Zero(3,1)),
        // upa_com(Vector::Zero(3,1)),
        // loa_com(Vector::Zero(3,1)),
        // hand_com(Vector::Zero(3,1)),
        // head_com(Vector::Zero(3,1));
        // foot_com(0) = 1.0;
        // shank_com(0) = 1.0;
        // thigh_com(0) = 1.0;
        // pelv_com(0) = 1.0;
        // torso_com(0) = 1.0;
        // upa_com(0) =  1.0;
        // loa_com(0) =  1.0;
        // hand_com(0) = 0.0;
        // head_com(0) = 0.0;


        /**
         * This function builds the pinocchio model_ and data_ attributes.
        */
        using namespace pinocchio;

        Matrix foot_in_mat = Matrix::Identity(3,3)*1e-15;
        foot_in_mat(2,2) = foot_in_zz;
        foot_in_mat = inertia_in_segment(foot_in_mat, foot_com, foot_m);
        Matrix shank_in_mat = Matrix::Identity(3,3)*1e-15;
        shank_in_mat(2,2) = shank_in_zz;
        shank_in_mat = inertia_in_segment(shank_in_mat, shank_com, shank_m);
        Matrix thigh_in_mat = Matrix::Identity(3,3)*1e-15;
        thigh_in_mat(2,2) = thigh_in_zz;
        thigh_in_mat = inertia_in_segment(thigh_in_mat, thigh_com, thigh_m);
        Matrix pelv_in_mat = Matrix::Identity(3,3)*1e-15;
        pelv_in_mat(2,2) = pelv_in_zz;
        pelv_in_mat = inertia_in_segment(pelv_in_mat, pelv_com, pelv_m);
        Matrix torso_in_mat = Matrix::Identity(3,3)*1e-15;
        torso_in_mat(2,2) = torso_in_zz;
        torso_in_mat = inertia_in_segment(torso_in_mat, torso_com, torso_m);
        Matrix upa_in_mat = Matrix::Identity(3,3)*1e-15;
        upa_in_mat(2,2) = upa_in_zz;
        upa_in_mat = inertia_in_segment(upa_in_mat, upa_com, upa_m);
        Matrix loa_in_mat = Matrix::Identity(3,3)*1e-15;
        loa_in_mat(2,2) = loa_in_zz;
        loa_in_mat = inertia_in_segment(loa_in_mat, loa_com, loa_m);
        Matrix hand_in_mat = Matrix::Identity(3,3)*1e-15;
        hand_in_mat(2,2) = hand_in_zz;
        hand_in_mat = inertia_in_segment(hand_in_mat, hand_com, hand_m);
        Matrix head_in_mat = Matrix::Identity(3,3)*1e-15;
        head_in_mat(2,2) = head_in_zz;
        head_in_mat = inertia_in_segment(head_in_mat, head_com, head_m);

        pinocchio::Inertia foot_in(foot_m, foot_com, foot_in_mat);
        pinocchio::Inertia shank_in(shank_m, shank_com, shank_in_mat);
        pinocchio::Inertia thigh_in(thigh_m, thigh_com, thigh_in_mat);
        pinocchio::Inertia pelv_in(pelv_m, pelv_com, pelv_in_mat);
        pinocchio::Inertia torso_in(torso_m, torso_com, torso_in_mat);
        pinocchio::Inertia upa_in(upa_m, upa_com, upa_in_mat);
        pinocchio::Inertia loa_in(loa_m, loa_com, loa_in_mat);
        pinocchio::Inertia hand_in(hand_m, hand_com, hand_in_mat);
        pinocchio::Inertia head_in(head_m, head_com, head_in_mat);
        pinocchio::Inertia null_inertia(1e-15, Vector::Ones(3)*1e-15, Matrix::Identity(3,3)*1e-15);

        pinocchio::JointModelVoid fixed_joint;
        pinocchio::JointModelRevoluteTpl<double, 0, 0> revolute_jointX;
        pinocchio::JointModelRevoluteTpl<double, 0, 1> revolute_jointY;
        pinocchio::JointModelRevoluteTpl<double, 0, 2> revolute_jointZ;

        Eigen::Matrix4d id4 = Eigen::Matrix4d::Identity(4,4);
        Vector pos_thigh_in_shank(Vector::Zero(3,1)),
               pos_pelv_in_thigh(Vector::Zero(3,1)),
               pos_torso_in_pelv(Vector::Zero(3,1)),
               pos_upa_in_torso(Vector::Zero(3,1)),
               pos_loa_in_upa(Vector::Zero(3,1)),
               pos_head_in_torso(Vector::Zero(3,1)),
               pos_hand_in_loa(Vector::Zero(3,1));
        pos_thigh_in_shank(0) = shank_l;
        pos_pelv_in_thigh(0) = thigh_l;
        pos_torso_in_pelv(0) = pelv_l;
        pos_head_in_torso(0) = torso_l;
        pos_upa_in_torso(0) = torso_l;
        pos_loa_in_upa(0) = upa_l;
        pos_hand_in_loa(0) = loa_l;

        /** Initialize transformation matrices */
        Eigen::Matrix4d foot_in_base = id4;
        Eigen::Matrix4d foot_com_in_foot = id4;
        foot_com_in_foot.block(0,3,3,1) = foot_com;
        Eigen::Matrix4d shank_in_foot = id4;
        Eigen::Matrix4d shank_com_in_shank = id4;
        shank_com_in_shank.block(0,3,3,1) = shank_com;
        Eigen::Matrix4d thigh_in_shank = id4;
        thigh_in_shank.block(0,3,3,1) = pos_thigh_in_shank;
        Eigen::Matrix4d thigh_com_in_thigh = id4;
        thigh_com_in_thigh.block(0,3,3,1) = thigh_com;
        Eigen::Matrix4d pelv_in_thigh = id4;
        pelv_in_thigh.block(0,3,3,1) = pos_pelv_in_thigh;
        Eigen::Matrix4d pelv_com_in_pelv = id4;
        pelv_com_in_pelv.block(0,3,3,1) = pelv_com;
        Eigen::Matrix4d torso_in_pelv = id4;
        torso_in_pelv.block(0,3,3,1) = pos_torso_in_pelv;
        Eigen::Matrix4d torso_com_in_torso = id4;
        torso_com_in_torso.block(0,3,3,1) = torso_com;
        Eigen::Matrix4d upa_in_torso = id4;
        upa_in_torso.block(0,3,3,1) = pos_upa_in_torso;
        Eigen::Matrix4d upa_com_in_upa = id4;
        upa_com_in_upa.block(0,3,3,1) = upa_com;
        Eigen::Matrix4d head_in_torso = id4;
        head_in_torso.block(0,3,3,1) = pos_head_in_torso;
        Eigen::Matrix4d head_com_in_head = id4;
        head_com_in_head.block(0,3,3,1) = head_com + pos_head_in_torso;
        Eigen::Matrix4d loa_in_upa = id4;
        loa_in_upa.block(0,3,3,1) = pos_loa_in_upa;
        Eigen::Matrix4d loa_com_in_loa = id4;
        loa_com_in_loa.block(0,3,3,1) = loa_com;
        Eigen::Matrix4d hand_in_loa = id4;
        hand_in_loa.block(0,3,3,1) = pos_hand_in_loa;
        Eigen::Matrix4d hand_com_in_hand = id4;
        hand_com_in_hand.block(0,3,3,1) = hand_com + pos_hand_in_loa;


        JointIndex joint_id;
        Frame frame = model_.frames[0];
        Frame foot_frame = Frame("foot", 0, 0, pinocchio::SE3Tpl<double, 0>(foot_in_base), FrameType::FIXED_JOINT, foot_in);
        FrameIndex prev_frame = model_.addFrame(foot_frame, true);

        Frame foot_com_frame = Frame("foot_com", 0, prev_frame, pinocchio::SE3Tpl<double, 0>(foot_com_in_foot), FrameType::FIXED_JOINT, null_inertia);
        model_.addFrame(foot_com_frame, false);

        joint_id = model_.addJoint(0, 
                       revolute_jointZ,  pinocchio::SE3Tpl<double, 0>(shank_in_foot), 
                       "ankle_rot_Z");
        FrameIndex jointFrameId = model_.addJointFrame(joint_id, prev_frame);
        model_.appendBodyToJoint(joint_id, shank_in, SE3::Identity());
        prev_frame = model_.addBodyFrame("shank", joint_id,  SE3::Identity(), (int)jointFrameId);
        Frame shank_com_frame = Frame("shank_com", joint_id, prev_frame, pinocchio::SE3Tpl<double, 0>(shank_com_in_shank), FrameType::FIXED_JOINT, null_inertia);
        model_.addFrame(shank_com_frame, false);

        frame = model_.frames[jointFrameId];
        joint_id = model_.addJoint(joint_id, 
                       revolute_jointZ,  pinocchio::SE3Tpl<double, 0>(thigh_in_shank), 
                       "knee_rot_Z");
        jointFrameId = model_.addJointFrame(joint_id, (int)prev_frame);
        model_.appendBodyToJoint(joint_id, thigh_in, SE3::Identity());
        prev_frame = model_.addBodyFrame("thigh", joint_id,  SE3::Identity(), (int)jointFrameId);
        Frame thigh_com_frame = Frame("thigh_com", joint_id, prev_frame, pinocchio::SE3Tpl<double, 0>(thigh_com_in_thigh), FrameType::FIXED_JOINT, null_inertia);
        model_.addFrame(thigh_com_frame, false);

        joint_id = model_.addJoint(joint_id, 
                       revolute_jointZ,  pinocchio::SE3Tpl<double, 0>(pelv_in_thigh), 
                       "hip_rot_Z");
        jointFrameId = model_.addJointFrame(joint_id, (int)prev_frame);
        model_.appendBodyToJoint(joint_id, pelv_in, SE3::Identity());
        prev_frame = model_.addBodyFrame("pelvis", joint_id,  SE3::Identity(), (int)jointFrameId);
        Frame pelv_com_frame = Frame("pelvis_com", joint_id, prev_frame, pinocchio::SE3Tpl<double, 0>(pelv_com_in_pelv), FrameType::FIXED_JOINT, null_inertia);
        model_.addFrame(pelv_com_frame, false);

        joint_id = model_.addJoint(joint_id, 
                       revolute_jointZ,  pinocchio::SE3Tpl<double, 0>(torso_in_pelv), 
                       "thoracic_rot_Z");
        jointFrameId = model_.addJointFrame(joint_id, (int)prev_frame);
        model_.appendBodyToJoint(joint_id, torso_in, SE3::Identity());
        prev_frame = model_.addBodyFrame("torso", joint_id,  SE3::Identity(), (int)jointFrameId);
        Frame torso_com_frame = Frame("torso_com", joint_id, prev_frame, pinocchio::SE3Tpl<double, 0>(torso_com_in_torso), FrameType::FIXED_JOINT, null_inertia);
        model_.addFrame(torso_com_frame, false);


        Frame head_frame = Frame("head", joint_id, prev_frame, pinocchio::SE3Tpl<double, 0>(head_in_torso), FrameType::FIXED_JOINT, head_in);
        prev_frame = model_.addFrame(head_frame, true);
        Frame head_com_frame = Frame("head_com", joint_id, prev_frame, pinocchio::SE3Tpl<double, 0>(head_com_in_head), FrameType::FIXED_JOINT, null_inertia);
        model_.addFrame(head_com_frame, false);

        frame = model_.frames[jointFrameId];
        joint_id = model_.addJoint(joint_id, 
                       revolute_jointZ,  pinocchio::SE3Tpl<double, 0>(upa_in_torso), 
                       "shoulder_rot_Z");
        jointFrameId = model_.addJointFrame(joint_id, (int)prev_frame);
        model_.appendBodyToJoint(joint_id, upa_in, SE3::Identity());
        prev_frame = model_.addBodyFrame("upperarm", joint_id,  SE3::Identity(), (int)jointFrameId);
        Frame upa_com_frame = Frame("upperarm_com", joint_id, prev_frame, pinocchio::SE3Tpl<double, 0>(upa_com_in_upa), FrameType::FIXED_JOINT, null_inertia);
        model_.addFrame(upa_com_frame, false);

        frame = model_.frames[jointFrameId];
        joint_id = model_.addJoint(joint_id, 
                       revolute_jointZ,  pinocchio::SE3Tpl<double, 0>(loa_in_upa), 
                       "elbow_rot_Z");
        jointFrameId = model_.addJointFrame(joint_id, (int)prev_frame);
        model_.appendBodyToJoint(joint_id, loa_in, SE3::Identity());
        prev_frame = model_.addBodyFrame("lowerarm", joint_id,  SE3::Identity(), (int)jointFrameId);
        Frame loa_com_frame = Frame("lowerarm_com", joint_id, prev_frame, pinocchio::SE3Tpl<double, 0>(loa_com_in_loa), FrameType::FIXED_JOINT, null_inertia);
        model_.addFrame(loa_com_frame, false);

        Frame hand_frame = Frame("hand", joint_id, prev_frame, pinocchio::SE3Tpl<double, 0>(hand_in_loa), FrameType::FIXED_JOINT, hand_in);
        prev_frame = model_.addFrame(hand_frame, true);
        Frame hand_com_frame = Frame("hand_com", joint_id, prev_frame, pinocchio::SE3Tpl<double, 0>(hand_com_in_hand), FrameType::FIXED_JOINT, null_inertia);
        model_.addFrame(hand_com_frame, false);

        // Create data required by the algorithms
        data_ = Data(model_);

        Vector true_g(Vector::Zero(3,1));
        true_g(1) = -9.81;

        model_.gravity.linear() = true_g;

        // Sample a random configuration
        Vector q = Vector::Zero(model_.nq);
        
        // std::cout << "Model.nq : " << model_.nq << std::endl;
        // std::cout << "q: " << q.transpose() << std::endl;

        // // Perform the forward kinematics over the kinematic tree
        // forwardKinematics(model_,data_,q);
        // updateFramePlacements(model_, data_);
        // // Print out the placement of each joint of the kinematic tree
        // for(JointIndex joint_id = 0; joint_id < (JointIndex)model_.njoints; ++joint_id)
        //     std::cout << std::setw(24) << std::left
        //             << model_.names[joint_id] << ": "
        //             << std::fixed << std::setprecision(5)
        //             << data_.oMi[joint_id].translation().transpose()
        //             <<"\n"
        //             // << data_.oMi[joint_id].rotation()
        //             << std::endl;
}
