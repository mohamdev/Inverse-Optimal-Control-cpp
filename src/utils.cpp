#include <utils.h>

Matrix vector_To_Eigen(std::vector<Vector> & vector_mat)
{
    /**
     * This function is used to convert a std::vector of Eigen column vectors
     * into an Eigen matrix
     * @param vector_mat The vector of Vectors to convert
     * @see Vector
    */
    Matrix eigenMat = Matrix::Zero(vector_mat[0].rows(), vector_mat.size());
    for (unsigned int i = 0; i < vector_mat.size() ; i++)
    {
        for(unsigned int j = 0; j<vector_mat[i].rows(); j++)
        {
            eigenMat(j,i) = vector_mat[i](j);
        }
    }
    return eigenMat;
}

std::vector<std::vector<double>> eigen_To_Vector(Matrix const& eigen_mat)
{
    /**
     * This function converts an eigen Matrix to a std::vector<vector<double>> 
     * @param eigen_mat Eigen matrix which will be converted to vector of vectors
    */
    std::vector<std::vector<double>> returned_vect(eigen_mat.rows(), std::vector<double>(eigen_mat.cols()));
    int j = 0;
    int i = 0;
    for (i=0; i<eigen_mat.rows(); i++)
    {
        for(j=0; j<eigen_mat.cols(); j++)
        {
            returned_vect[i][j] = eigen_mat(i,j);
        }
    }
    return returned_vect;
}

void write_meas_data_vicon(std::string const& path, std::vector<Vector> const& meas_data)
{
    std::ofstream file(path);

    for(int i = 0; i < meas_data.size(); i++)
    {
        file << meas_data[i](0) << "," << meas_data[i](1) << "," << meas_data[i](2) << 
        "," << meas_data[i](3) << "," << meas_data[i](4) << "," << meas_data[i](5) <<
        "," << meas_data[i](6) << "," << meas_data[i](7) << "," << meas_data[i](8) <<
        "," << meas_data[i](9) << "," << meas_data[i](10) << "," << meas_data[i](11) <<
        "," << meas_data[i](12) << "," << meas_data[i](13) << "," << meas_data[i](14) <<
        "," << meas_data[i](15) << "," << meas_data[i](16) << "," << meas_data[i](17) << 
        "," << meas_data[i](18) << "," << meas_data[i](19) << "," << meas_data[i](20) <<
        "," << meas_data[i](21) << "," << meas_data[i](22) << "," << meas_data[i](23) <<
        "," << meas_data[i](24) << "," << meas_data[i](25) << "," << meas_data[i](26) <<
        "," << meas_data[i](27) << "," << meas_data[i](28) << "," << meas_data[i](29) <<
        "," << meas_data[i](30) << "," << meas_data[i](31) << "," << meas_data[i](32) <<
        "," << meas_data[i](33) << "," << meas_data[i](34) << "," << meas_data[i](35) <<
        "," << meas_data[i](36) << "," << meas_data[i](37) << "," << meas_data[i](38) <<std::endl;
    }
    file.close();
}

void write_meas_data_4vimus(std::string const& path, std::vector<Vector> const& meas_data)
{
    std::ofstream file(path);

    for(int i = 0; i < meas_data.size(); i++)
    {
        file << meas_data[i](0) << "," << meas_data[i](1) << "," << meas_data[i](2) << 
        "," << meas_data[i](3) << "," << meas_data[i](4) << "," << meas_data[i](5) <<
        "," << meas_data[i](6) << "," << meas_data[i](7) << "," << meas_data[i](8) <<
        "," << meas_data[i](9) << "," << meas_data[i](10) << "," << meas_data[i](11) <<
        "," << meas_data[i](12) << "," << meas_data[i](13) << "," << meas_data[i](14) <<
        "," << meas_data[i](15) << "," << meas_data[i](16) << "," << meas_data[i](17) << 
        "," << meas_data[i](18) << "," << meas_data[i](19) << "," << meas_data[i](20) <<
        "," << meas_data[i](21) << "," << meas_data[i](22) << "," << meas_data[i](23) <<
        "," << meas_data[i](24) << "," << meas_data[i](25) << "," << meas_data[i](26) <<
        "," << meas_data[i](27) << "," << meas_data[i](28) << "," << meas_data[i](29) <<
        "," << meas_data[i](30) << "," << meas_data[i](31) << "," << meas_data[i](32) <<
        "," << meas_data[i](33) << "," << meas_data[i](34) << "," << meas_data[i](35) <<
        "," << meas_data[i](36) << "," << meas_data[i](37) << "," << meas_data[i](38) << 
        "," << meas_data[i](39) << "," << meas_data[i](40) << "," << meas_data[i](41) << 
        "," << meas_data[i](42) << "," << meas_data[i](43) << "," << meas_data[i](44) << 
        "," << meas_data[i](45) << "," << meas_data[i](46) << "," << meas_data[i](47) << 
        "," << meas_data[i](48) << "," << meas_data[i](49) << "," << meas_data[i](50) << 
        "," << meas_data[i](51) << std::endl;
    }
    file.close();
}

void write_anat_orientation_data(std::string const& path,
                                std::vector<Vector> const& quats_trunk,
                                std::vector<Vector> const& quats_upperarm,
                                std::vector<Vector> const& quats_lowerarm,
                                std::vector<Vector> const& quats_hand)
{
    std::ofstream file(path);

    for(int i = 0; i < quats_trunk.size(); i++)
    {
        file << quats_trunk[i](0) << "," << quats_trunk[i](1) << "," << quats_trunk[i](2) << 
        "," << quats_trunk[i](3) << "," <<  quats_upperarm[i](0) << "," << quats_upperarm[i](1) << "," << quats_upperarm[i](2) << 
        "," << quats_upperarm[i](3) << "," << quats_lowerarm[i](0) << "," << quats_lowerarm[i](1) << "," << quats_lowerarm[i](2) << 
        "," << quats_lowerarm[i](3) << "," << quats_hand[i](0) << "," << quats_hand[i](1) << "," << quats_hand[i](2) << 
        "," << quats_hand[i](3) << std::endl;
    }

    file.close();
}

void write_anat_orientation_data_2vimus(std::string const& path,
                                std::vector<Vector> const& quats_trunk,
                                std::vector<Vector> const& quats_hand)
{
    std::ofstream file(path);

    for(int i = 0; i < quats_trunk.size(); i++)
    {
        file << quats_trunk[i](1) << "," << quats_trunk[i](2) << "," << quats_trunk[i](3) << 
        "," << quats_trunk[i](0) << "," << quats_hand[i](1) << "," << quats_hand[i](2) << "," << quats_hand[i](3) << 
        "," << quats_hand[i](0) << std::endl;
    }

    file.close();
}

void write_anat_pos_orientation_data_2vimus(std::string const& path,
                                std::vector<Vector> const& quats_trunk,
                                std::vector<Vector> const& quats_hand,
                                std::vector<Vector> const& pos_trunk,
                                std::vector<Vector> const& pos_hand)
{
    std::ofstream file(path);

    for(int i = 0; i < quats_trunk.size(); i++)
    {
        file << pos_trunk[i](0) << "," << pos_trunk[i](1) << "," << pos_trunk[i](2) << "," << quats_trunk[i](1) << "," << quats_trunk[i](2) << "," 
        << quats_trunk[i](3) << "," << quats_trunk[i](0) << "," << pos_hand[i](0) << "," << pos_hand[i](1) << "," << pos_hand[i](2) << "," 
        << quats_hand[i](1) << "," << quats_hand[i](2) << "," << quats_hand[i](3) << "," << quats_hand[i](0) << std::endl;
    }

    file.close();
}

void write_2vimus_data(std::string const& path,
                        std::vector<Vector> const& quats_trunk,
                        std::vector<Vector> const& quats_hand,
                        std::vector<Vector> const& pos_trunk,
                        std::vector<Vector> const& pos_hand,
                        std::vector<Vector> const& gyr_trunk,
                        std::vector<Vector> const& gyr_hand,
                        std::vector<Vector> const& acc_trunk,
                        std::vector<Vector> const& acc_hand)
{
    std::ofstream file(path);

    for(int i = 0; i < quats_trunk.size(); i++)
    {
        file << pos_trunk[i](0) << "," << pos_trunk[i](1) << "," << pos_trunk[i](2) 
        << "," << quats_trunk[i](1) << "," << quats_trunk[i](2) << "," << quats_trunk[i](3) << "," << quats_trunk[i](0) 
        << "," << gyr_trunk[i](0) << "," << gyr_trunk[i](1) << "," << gyr_trunk[i](2) 
        << "," << acc_trunk[i](0) << "," << acc_trunk[i](1) << "," << acc_trunk[i](2) 
        << "," << pos_hand[i](0) << "," << pos_hand[i](1) << "," << pos_hand[i](2) 
        << ","  << quats_hand[i](1) << "," << quats_hand[i](2) << "," << quats_hand[i](3) << "," << quats_hand[i](0) 
        << "," << gyr_hand[i](0) << "," << gyr_hand[i](1) << "," << gyr_hand[i](2) 
        << "," << acc_hand[i](0) << "," << acc_hand[i](1) << "," << acc_hand[i](2) << std::endl;
    }

    file.close();
}

Vector rotm_2_euler(Matrix const&R)
{

    double sy = sqrt(R(0,0) * R(0,0) +  R(1,0) * R(1,0) );

    bool singular = sy < 1e-6; // If

    float x, y, z;
    if (!singular)
    {
        x = std::atan2((double)R(2,1) , (double)R(2,2));
        y = std::atan2(-R(2,0), sy);
        z = std::atan2(R(1,0), R(0,0));
    }
    else
    {
        x = std::atan2(-R(1,2), R(1,1));
        y = std::atan2(-R(2,0), sy);
        z = 0;
    }
    Vector v(Vector::Zero(3,1));
    v(0) = x;
    v(1) = y;
    v(2) = z;
    return v;
}

void write_ioc_solution(std::string const& path, Vector const& weights)
{
    std::ofstream file(path);

    file << weights(0) << "," << weights(1) << "," << weights(2) << 
    "," << weights(3) << "," << weights(4) << "," << weights(5) <<  
    "," << weights(6) << "," << weights(7) << "," << weights(8) << std::endl;

    file.close();
}

void write_q_data(std::string const& path, std::vector<Vector> const& q_data)
{
    std::ofstream file(path);

    for(int i = 0; i < q_data.size(); i++)
    {
        file << q_data[i](0) << "," << q_data[i](1) << "," << q_data[i](2) << 
        "," << q_data[i](3) << "," << q_data[i](4) << "," << q_data[i](5) <<
        "," << q_data[i](6) << "," << q_data[i](7) << "," << q_data[i](8) <<
        "," << q_data[i](9) << "," << q_data[i](10) << "," << q_data[i](11) <<
        "," << q_data[i](12) << std::endl;
    }
    file.close();
}

void write_q_data_6dof(std::string const& path, std::vector<Vector> const& q_data)
{
    std::ofstream file(path);

    for(int i = 0; i < q_data.size(); i++)
    {
        file << q_data[i](0) << "," << q_data[i](1) << "," << q_data[i](2) << 
        "," << q_data[i](3) << "," << q_data[i](4) << "," << q_data[i](5) << std::endl;
    }
    file.close();
}

void write_q_data_6dof(std::string const& path, std::vector<CGVector> const& q_data)
{
    std::ofstream file(path);

    for(int i = 0; i < q_data.size(); i++)
    {
        file << q_data[i](0) << "," << q_data[i](1) << "," << q_data[i](2) << 
        "," << q_data[i](3) << "," << q_data[i](4) << "," << q_data[i](5) << std::endl;
    }
    file.close();
}

std::vector<Vector> read_q_data(std::string const& path)
{
    using namespace std;

    std::vector<Vector> q_traj;
    fstream fin;

    string line;
    // Open an existing file
    fin.open(path, ios::in);

    while(std::getline(fin, line))
    {
        // Read the Data from the file
        // as String Vector
        std::vector<string> row;
        string word;

        /**Read anat points data**/
        std::vector<Vector> anat_points;

        row.clear();

        // used for breaking words
        stringstream s(line);

        // read every column data of a row and
        // store it in a string variable, 'word'
        while (getline(s, word, ',')) 
        {
            // add all the column data
            // of a row to a vector
            row.push_back(word);
            // std::cout << "Word : " << word << std::endl;
        }

        Vector q_curr(Vector::Zero(13,1));
        for(int i = 0; i < 13; i+=1)
        {
            q_curr(i) = std::stod(row[i]);
        }
        q_traj.push_back(q_curr);
    }

    fin.close();
    return q_traj;
}

std::vector<Vector> read_q7_data(std::string const& path)
{
    using namespace std;

    std::vector<Vector> q_traj;
    fstream fin;

    string line;
    // Open an existing file
    fin.open(path, ios::in);

    while(std::getline(fin, line))
    {
        // Read the Data from the file
        // as String Vector
        std::vector<string> row;
        string word;

        /**Read anat points data**/
        std::vector<Vector> anat_points;

        row.clear();

        // used for breaking words
        stringstream s(line);

        // read every column data of a row and
        // store it in a string variable, 'word'
        while (getline(s, word, ',')) 
        {
            // add all the column data
            // of a row to a vector
            row.push_back(word);
            // std::cout << "Word : " << word << std::endl;
        }

        Vector q_curr(Vector::Zero(7,1));
        for(int i = 0; i < 7; i+=1)
        {
            q_curr(i) = std::stod(row[i]);
        }
        q_traj.push_back(q_curr);
    }

    fin.close();
    return q_traj;
}

std::vector<Vector> read_q6_data(std::string const& path)
{
    using namespace std;

    std::vector<Vector> q_traj;
    fstream fin;

    string line;
    // Open an existing file
    fin.open(path, ios::in);

    int k = 0;
    while(std::getline(fin, line))
    {
        if (k>0)
        {
            // Read the Data from the file
            // as String Vector
            std::vector<string> row;
            string word;

            /**Read anat points data**/
            std::vector<Vector> anat_points;

            row.clear();

            // used for breaking words
            stringstream s(line);

            // read every column data of a row and
            // store it in a string variable, 'word'
            while (getline(s, word, ',')) 
            {
                // add all the column data
                // of a row to a vector
                row.push_back(word);
                // std::cout << "Word : " << word << std::endl;
            }

            Vector q_curr(Vector::Zero(6,1));
            for(int i = 0; i < 6; i+=1)
            {
                q_curr(i) = std::stod(row[i]);
            }
            q_traj.push_back(q_curr);
        }
        k++;

    }

    fin.close();
    return q_traj;
}


void write_calib_data(std::string const& path_vimu_local,
                      std::string const& path_anat_local,
                      std::map<std::string, Matrix> calib_data_map)
{
    Matrix vimu_trunk = calib_data_map["VIMU_TRUNK"];
    Matrix vimu_upperarm = calib_data_map["VIMU_UPPERARM"];
    Matrix vimu_lowerarm = calib_data_map["VIMU_LOWERARM"];
    Matrix vimu_hand = calib_data_map["VIMU_HAND"];
    Vector GH = calib_data_map["GH"];
    Vector ELBOW = calib_data_map["ELBOW"];
    Vector WRIST = calib_data_map["WRIST"];
    

    std::ofstream file_vimu(path_vimu_local);
    file_vimu << vimu_trunk(0,0) << "," << vimu_trunk(0,1) << "," << vimu_trunk(0,2) 
    << "," << vimu_trunk(1,0) << "," << vimu_trunk(1,1) << "," << vimu_trunk(1,2)
    << "," << vimu_trunk(2,0) << "," << vimu_trunk(2,1) << "," << vimu_trunk(2,2)
    << "," << vimu_trunk(0,3) << "," << vimu_trunk(1,3) << "," << vimu_trunk(2,3) << std::endl; 
    file_vimu << vimu_upperarm(0,0) << "," << vimu_upperarm(0,1) << "," << vimu_upperarm(0,2) 
    << "," << vimu_upperarm(1,0) << "," << vimu_upperarm(1,1) << "," << vimu_upperarm(1,2)
    << "," << vimu_upperarm(2,0) << "," << vimu_upperarm(2,1) << "," << vimu_upperarm(2,2)
    << "," << vimu_upperarm(0,3) << "," << vimu_upperarm(1,3) << "," << vimu_upperarm(2,3) << std::endl; 
    file_vimu << vimu_lowerarm(0,0) << "," << vimu_lowerarm(0,1) << "," << vimu_lowerarm(0,2) 
    << "," << vimu_lowerarm(1,0) << "," << vimu_lowerarm(1,1) << "," << vimu_lowerarm(1,2)
    << "," << vimu_lowerarm(2,0) << "," << vimu_lowerarm(2,1) << "," << vimu_lowerarm(2,2)
    << "," << vimu_lowerarm(0,3) << "," << vimu_lowerarm(1,3) << "," << vimu_lowerarm(2,3) << std::endl; 
    file_vimu << vimu_hand(0,0) << "," << vimu_hand(0,1) << "," << vimu_hand(0,2) 
    << "," << vimu_hand(1,0) << "," << vimu_hand(1,1) << "," << vimu_hand(1,2)
    << "," << vimu_hand(2,0) << "," << vimu_hand(2,1) << "," << vimu_hand(2,2)
    << "," << vimu_hand(0,3) << "," << vimu_hand(1,3) << "," << vimu_hand(2,3) << std::endl;

    file_vimu.close();

    std::ofstream file_anat(path_anat_local);
    file_anat << GH(0) << "," << GH(1) << "," << GH(2) << std::endl;
    file_anat << ELBOW(0) << "," << ELBOW(1) << "," << ELBOW(2) << std::endl;
    file_anat << WRIST(0) << "," << WRIST(1) << "," << WRIST(2) << std::endl;
    file_anat.close();
}


void write_calib_data_2vimus(std::string const& path_vimu_local,
                      std::string const& path_anat_local,
                      std::map<std::string, Matrix> calib_data_map)
{
    Matrix vimu_trunk = calib_data_map["VIMU_TRUNK"];
    Matrix vimu_hand = calib_data_map["VIMU_HAND"];
    Vector GH = calib_data_map["GH"];
    Vector ELBOW = calib_data_map["ELBOW"];
    Vector WRIST = calib_data_map["WRIST"];
    

    std::ofstream file_vimu(path_vimu_local);
    file_vimu << vimu_trunk(0,0) << "," << vimu_trunk(0,1) << "," << vimu_trunk(0,2) 
    << "," << vimu_trunk(1,0) << "," << vimu_trunk(1,1) << "," << vimu_trunk(1,2)
    << "," << vimu_trunk(2,0) << "," << vimu_trunk(2,1) << "," << vimu_trunk(2,2)
    << "," << vimu_trunk(0,3) << "," << vimu_trunk(1,3) << "," << vimu_trunk(2,3) << std::endl; 
    file_vimu << vimu_hand(0,0) << "," << vimu_hand(0,1) << "," << vimu_hand(0,2) 
    << "," << vimu_hand(1,0) << "," << vimu_hand(1,1) << "," << vimu_hand(1,2)
    << "," << vimu_hand(2,0) << "," << vimu_hand(2,1) << "," << vimu_hand(2,2)
    << "," << vimu_hand(0,3) << "," << vimu_hand(1,3) << "," << vimu_hand(2,3) << std::endl;

    file_vimu.close();

    std::ofstream file_anat(path_anat_local);
    file_anat << GH(0) << "," << GH(1) << "," << GH(2) << std::endl;
    file_anat << ELBOW(0) << "," << ELBOW(1) << "," << ELBOW(2) << std::endl;
    file_anat << WRIST(0) << "," << WRIST(1) << "," << WRIST(2) << std::endl;
    file_anat.close();
}

Vector compute_xyz_rmse(std::vector<Vector> est_xyz,
                        std::vector<Vector> ref_xyz)
{
    /**
     * @brief Takes as input two collections of 3D vectors, and returns a 4D 
     * vector containing rmse on x, y, z and mean
     */

    Matrix error_mat(Matrix::Zero(3,est_xyz.size()));

    for(int i = 0; i < error_mat.cols(); i++)
    {
        error_mat.col(i) = est_xyz[i] - ref_xyz[i];
    }

    for(int i = 0; i < error_mat.rows(); i++)
    {
        // std::cout << "before : " << error_mat.col(i).transpose() << std::endl;
        for(int j = 0; j < error_mat.cols(); j++)
        {
            error_mat(i,j) = error_mat(i,j)*error_mat(i,j);
        }
        // std::cout << "after : " << error_mat.col(i).transpose() << std::endl;
    }

    double x_rms = std::sqrt(error_mat.row(0).mean());
    double y_rms = std::sqrt(error_mat.row(1).mean());
    double z_rms = std::sqrt(error_mat.row(2).mean());
    double mean = (x_rms + y_rms + z_rms)/3.0;

    std::cout << "x_rms : " << x_rms << " m" << std::endl;
    std::cout << "y_rms : " << y_rms << " m" << std::endl;
    std::cout << "z_rms : " << z_rms << " m" << std::endl;
    std::cout << "mean : " << mean << " m" << std::endl;

    Vector rms_vec(Vector::Zero(4,1));
    rms_vec(0) = x_rms; rms_vec(1) = y_rms; rms_vec(2) = z_rms; rms_vec(3) = mean;

    return rms_vec;
}


Vector compute_xyz_covariance(std::vector<Vector> est_xyz,
                              std::vector<Vector> ref_xyz)
{
    /**
     * @brief Takes as input two collections of 3D vectors, and returns a 4D 
     * vector containing rmse on x, y, z and mean
     */
    Matrix X_mat(Matrix::Zero(3,est_xyz.size()));
    Matrix Y_mat(Matrix::Zero(3,ref_xyz.size()));


    for(int i = 0; i < X_mat.rows(); i++)
    {
        for(int j = 0; j < Y_mat.cols(); j++)
        {
            X_mat(i,j) = est_xyz[j][i];
            Y_mat(i,j) = ref_xyz[j][i];
        }
    }

    double x_0_mean = X_mat.row(0).mean();
    double x_1_mean = X_mat.row(1).mean();
    double x_2_mean = X_mat.row(2).mean();
    double y_0_mean = Y_mat.row(0).mean();
    double y_1_mean = Y_mat.row(1).mean();
    double y_2_mean = Y_mat.row(2).mean();

    for(int j = 0; j < Y_mat.cols(); j++)
    {
        X_mat(0,j) = X_mat(0,j) - x_0_mean;
        X_mat(1,j) = X_mat(1,j) - x_1_mean;
        X_mat(2,j) = X_mat(2,j) - x_2_mean;
        Y_mat(0,j) = Y_mat(0,j) - y_0_mean;
        Y_mat(1,j) = Y_mat(1,j) - y_1_mean;
        Y_mat(2,j) = Y_mat(2,j) - y_2_mean;
    }

    Matrix cov_mat(Matrix::Zero(3,est_xyz.size()));

    for(int i = 0; i < cov_mat.rows(); i++)
    {
        // std::cout << "before : " << error_mat.col(i).transpose() << std::endl;
        for(int j = 0; j < cov_mat.cols(); j++)
        {
            cov_mat(i,j) = X_mat(i,j)*Y_mat(i,j);
        }
        // std::cout << "after : " << error_mat.col(i).transpose() << std::endl;
    }

    double x_rms = cov_mat.row(0).mean();
    double y_rms = cov_mat.row(1).mean();
    double z_rms = cov_mat.row(2).mean();
    double mean = (x_rms + y_rms + z_rms)/3.0;

    std::cout << "x cov : " << x_rms << " m2" << std::endl;
    std::cout << "y cov : " << y_rms << " m2" << std::endl;
    std::cout << "z cov : " << z_rms << " m2" << std::endl;
    std::cout << "mean cov : " << mean << " m2" << std::endl;

    Vector rms_vec(Vector::Zero(4,1));
    rms_vec(0) = x_rms; rms_vec(1) = y_rms; rms_vec(2) = z_rms; rms_vec(3) = mean;

    return rms_vec;
}

Vector compute_xyz_PC(std::vector<Vector> est_xyz,
                      std::vector<Vector> ref_xyz)
{
    /**
     * @brief Takes as input two collections of 3D vectors, and returns a 4D 
     * vector containing rmse on x, y, z and mean
     */
    Matrix X_mat(Matrix::Zero(3,est_xyz.size()));
    Matrix Y_mat(Matrix::Zero(3,ref_xyz.size()));


    for(int i = 0; i < X_mat.rows(); i++)
    {
        for(int j = 0; j < Y_mat.cols(); j++)
        {
            X_mat(i,j) = est_xyz[j][i];
            Y_mat(i,j) = ref_xyz[j][i];
        }
    }

    double PC_x_rms = calculatePC(X_mat.row(0).transpose(), Y_mat.row(0).transpose());
    double PC_y_rms = calculatePC(X_mat.row(1).transpose(), Y_mat.row(1).transpose());
    double PC_z_rms = calculatePC(X_mat.row(2).transpose(), Y_mat.row(2).transpose());
    double PC_mean = (PC_x_rms + PC_y_rms + PC_z_rms)/3.0;

    std::cout << "x_rms : " << PC_x_rms << std::endl;
    std::cout << "y_rms : " << PC_y_rms << std::endl;
    std::cout << "z_rms : " << PC_z_rms << std::endl;
    std::cout << "mean : " << PC_mean << std::endl;

    Vector rms_vec(Vector::Zero(4,1));
    rms_vec(0) = PC_x_rms; rms_vec(1) = PC_y_rms; rms_vec(2) = PC_z_rms; rms_vec(3) = PC_mean;

    return rms_vec;
}

void write_arm_3D_pos_data(std::string const& path,
                            std::vector<Vector> const& pos_shoulder,
                            std::vector<Vector> const& pos_elbow,
                            std::vector<Vector> const& pos_wrist)
{

    std::ofstream file(path);

    for(int i = 0; i < pos_shoulder.size(); i++)
    {
        file << pos_shoulder[i](0) << "," << pos_shoulder[i](1) << "," << pos_shoulder[i](2) << 
        "," <<  pos_elbow[i](0) << "," << pos_elbow[i](1) << "," << pos_elbow[i](2) << 
        "," << pos_wrist[i](0) << "," << pos_wrist[i](1) << "," << pos_wrist[i](2) << std::endl;
    }

    file.close();
}

double calculatePC(Vector v1, Vector v2){
        double s12 = 0.0,s1 =0.0,s2 =0.0,s1s =0.0,s2s =0.0;
        double pc = 0.0;
        int n = v1.rows();
        for(int i=0;i < n;i++){
            s12 += v1(i)*v2(i);
            s1 += v1(i); s1s += v1(i)*v1(i);
            s2 += v2(i); s2s += v2(i)*v2(i);
        }
        pc = (n*s12 - s1*s2)/std::sqrt((n*s1s-s1*s1)*(n*s2s-s2*s2));
        return pc;
}

void read_control_points(std::vector<double> & t_cp,
                         std::vector<Vector> & q_cp,
                         std::vector<Vector> & dq_cp,
                         std::string const& cp_path)
{
    std::fstream fin;

    std::string line;
    // Open an existing file
    fin.open(cp_path, std::ios::in);

    int k =0;
    while(std::getline(fin, line))
    {
        if(k>0)
        {
            // Read the Data from the file
            // as String Vector
            std::vector<std::string> row;
            std::string word;

            /**Read anat points data**/
            std::vector<Vector> anat_points;

            row.clear();

            // used for breaking words
            std::stringstream s(line);

            // read every column data of a row and
            // store it in a string variable, 'word'
            while (getline(s, word, ',')) 
            {
                // add all the column data
                // of a row to a vector
                row.push_back(word);
                // std::cout << "Word : " << word << std::endl;
            }
            double t = std::stod(row[0]);
            Vector q_curr(Vector::Zero(6,1)),
                dq_curr(Vector::Zero(6,1));
            for(int i = 1; i < 7; i+=1)
            {
                q_curr(i-1) = std::stod(row[i]);
                dq_curr(i-1) = std::stod(row[i+6]);
            }
            t_cp.push_back(t);
            q_cp.push_back(q_curr);
            dq_cp.push_back(dq_curr);
        }
        k++;
    }

    fin.close();
}

void read_box_params(double & mass,
                     Vector & pos,
                     std::string const& box_params_path)
{
    std::fstream fin;

    std::string line;
    // Open an existing file
    fin.open(box_params_path, std::ios::in);

    int k =0;
    while(std::getline(fin, line))
    {
        if(k>0)
        {
            // Read the Data from the file
            // as String Vector
            std::vector<std::string> row;
            std::string word;

            /**Read anat points data**/
            std::vector<Vector> anat_points;

            row.clear();

            // used for breaking words
            std::stringstream s(line);

            // read every column data of a row and
            // store it in a string variable, 'word'
            while (getline(s, word, ',')) 
            {
                // add all the column data
                // of a row to a vector
                row.push_back(word);
                // std::cout << "Word : " << word << std::endl;
            }
            mass = std::stod(row[2]);
            pos(0) = std::stod(row[3]);
            pos(1) = std::stod(row[4]);
        }
        k++;
    }

    fin.close();
}

void read_box_params(double & mass,
                     Vector & pos,
                     double & box_width,
                     std::string const& box_params_path)
{
    std::fstream fin;

    std::string line;
    // Open an existing file
    fin.open(box_params_path, std::ios::in);

    int k =0;
    while(std::getline(fin, line))
    {
        if(k>0)
        {
            // Read the Data from the file
            // as String Vector
            std::vector<std::string> row;
            std::string word;

            /**Read anat points data**/
            std::vector<Vector> anat_points;

            row.clear();

            // used for breaking words
            std::stringstream s(line);

            // read every column data of a row and
            // store it in a string variable, 'word'
            while (getline(s, word, ',')) 
            {
                // add all the column data
                // of a row to a vector
                row.push_back(word);
                // std::cout << "Word : " << word << std::endl;
            }
            box_width = std::stod(row[0]);
            mass = std::stod(row[2]);
            pos(0) = std::stod(row[3]);
            pos(1) = std::stod(row[4]);
        }
        k++;
    }

    fin.close();
}

void read_table_params(double & tw,
                       double & th,
                       Vector & pos,
                       std::string const& table_params_path)
{
    std::fstream fin;

    std::string line;
    // Open an existing file
    fin.open(table_params_path, std::ios::in);

    int k =0;
    while(std::getline(fin, line))
    {
        if(k>0)
        {
            // Read the Data from the file
            // as String Vector
            std::vector<std::string> row;
            std::string word;

            /**Read anat points data**/
            std::vector<Vector> anat_points;

            row.clear();

            // used for breaking words
            std::stringstream s(line);

            // read every column data of a row and
            // store it in a string variable, 'word'
            while (getline(s, word, ',')) 
            {
                // add all the column data
                // of a row to a vector
                row.push_back(word);
                // std::cout << "Word : " << word << std::endl;
            }
            tw = std::stod(row[5]);
            th = std::stod(row[6]);
            pos(0) = std::stod(row[7]);
            pos(1) = std::stod(row[8]);
        }
        k++;
    }
    fin.close();
}

void read_final_hand_pos(Vector & pos,
                         std::string const& box_params_path)
{
    std::fstream fin;

    std::string line;
    // Open an existing file
    fin.open(box_params_path, std::ios::in);

    int k =0;
    while(std::getline(fin, line))
    {
        if(k>0)
        {
            // Read the Data from the file
            // as String Vector
            std::vector<std::string> row;
            std::string word;

            /**Read anat points data**/
            std::vector<Vector> anat_points;

            row.clear();

            // used for breaking words
            std::stringstream s(line);

            // read every column data of a row and
            // store it in a string variable, 'word'
            while (getline(s, word, ',')) 
            {
                // add all the column data
                // of a row to a vector
                row.push_back(word);
                // std::cout << "Word : " << word << std::endl;
            }
            pos(0) = std::stod(row[11]);
            pos(1) = std::stod(row[12]);
        }
        k++;
    }

    fin.close();
}
