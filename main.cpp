/**
* @ingroup idyntree_tutorials
*
* \defgroup icub_genericChainController Controller for a
* Generic Kinematic Chain
*
* A tutorial on how to obtain the parametric model for
*   iCub arm FT sensor measurement.
*
* \author Silvio Traversaro
* \author Raffaello Camoriano
*
* CopyPolicy: Released under the terms of GPL 2.0 or later
*/

///////////////////////////////////////////////////////////
//  Headers
///////////////////////////////////////////////////////////

// YARP headers
#include <yarp/os/ResourceFinder.h>
#include <yarp/os/LogStream.h>
#include <yarp/os/Property.h>

// iDynTree headers
#include <kdl_codyco/regressors/dynamicRegressorGenerator.hpp>
#include "kdl_format_io/urdf_import.hpp"
#include "kdl_format_io/urdf_export.hpp"
#include <kdl_codyco/regressors/DynamicDatasetFile.hpp>
#include <kdl_codyco/regressors/dynamicRegressorGenerator.hpp>
#include <kdl_codyco/regressor_loops.hpp>
#include <kdl_codyco/regressor_utils.hpp>
#include <kdl_codyco/six_axis_ft_sensor.hpp>

// iCub headers
#include <iCub/ctrl/filters.h>
#include <iCub/ctrl/adaptWinPolyEstimator.h>

#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <fstream>

// Others
#include "multitaskRecursiveLinearEstimator.h"
#include "multitaskSVDLinearEstimator.h"

///////////////////////////////////////////////////////////
// Namespaces
///////////////////////////////////////////////////////////

using namespace std;
using namespace KDL;
using namespace KDL::CoDyCo;
using namespace KDL::CoDyCo::Regressors;


///////////////////////////////////////////////////////////
// Defines
///////////////////////////////////////////////////////////

// DoFs
#define MODEL_DOFS 32
#define DATASET_DOFS 41

// Torque indexes
#define LEFT_ELBOW_TORQUE_INDEX 2
#define RIGHT_ELBOW_TORQUE_INDEX 8
#define LEFT_SHOULDER_PITCH_TORQUE_INDEX 6

// FT indexes
#define LEFT_ARM_FT_INDEX 0
#define RIGHT_ARM_FT_INDEX 1

// Gains
#define LEFT_ELBOW_GAIN 3.5290
#define LEFT_SHOULDER_PITCH_GAIN 2.4191

///////////////////////////////////////////////////////////
/** Tolerance for considering two values equal */
///////////////////////////////////////////////////////////
const double ZERO_TOL = 1e-5;

///////////////////////////////////////////////////////////
/// Useful functions
///////////////////////////////////////////////////////////

/**
 * The set of joints that are found in the dataset (i.e. head+torso+right_arm+left_arm) is different from
 * the joints found in the model, so we have to convert them
 */
bool convertJointsFromDatasetToModel(const KDL::JntArray & q_dataset, KDL::JntArray & q_model, bool verbose = false)
{
    if( q_dataset.rows() != DATASET_DOFS ) {
        if( verbose ) cout << "Q from dataset has size " << q_dataset.rows() << std::endl;
    }

    assert(q_dataset.rows() == DATASET_DOFS);
    assert(q_model.rows() == MODEL_DOFS);
    //Joints in dataset are : Torso (3), head(6), left arm (16), right_arm (16)
    //Joints in the model are : Torso (3), head(3), left arm (7), right arm (7), left_leg (6), right_leg(6)
    
    SetToZero(q_model);
    
    //Torso should be inverted
    q_model(0) = q_dataset(2);
    q_model(1) = q_dataset(1);
    q_model(2) = q_dataset(0);

    //Copyng head
    for(int i = 0; i < 3; i++ ) {
        q_model(3+i) = q_dataset(3+i);
    }
    
    //Copyng left arm
    for(int i = 0; i < 7; i++ ) {
        q_model(3+3+i) = q_dataset(3+6+i);
    }
    
    //Copyng right arm 
    for(int i = 0; i < 7; i++ ) {
        q_model(3+3+7+i) = q_dataset(3+6+16+i);
    }
    
    return true;
}

bool toYarp(const KDL::Wrench & ft, yarp::sig::Vector & ft_yrp)
{
    if( ft_yrp.size() != 6 ) { ft_yrp.resize(6); }
    for(int i=0; i < 6; i++ ) {
        ft_yrp[i] = ft[i];
    }
}

bool toKDL(const yarp::sig::Vector & ft_yrp, KDL::Wrench & ft)
{
    for(int i=0; i < 6; i++ ) {
        ft[i] = ft_yrp[i];
    }
}

bool toYarp(const Eigen::VectorXd & vec_eigen, yarp::sig::Vector & vec_yrp)
{
    if( vec_yrp.size() != vec_eigen.size() ) { vec_yrp.resize(vec_eigen.size()); }
    if( memcpy(vec_yrp.data(),vec_eigen.data(),sizeof(double)*vec_eigen.size()) != NULL ) {
        return true;
    } else {
        return false;
    }
}

bool toEigen(const yarp::sig::Vector & vec_yrp, Eigen::VectorXd & vec_eigen)
{
 if( vec_yrp.size() != vec_eigen.size() ) { vec_eigen.resize(vec_yrp.size()); }
    if( memcpy(vec_eigen.data(),vec_yrp.data(),sizeof(double)*vec_eigen.size()) != NULL ) {
        return true;
    } else {
        return false;
    }
}

SensorsTree generateSensorsTree(const KDL::CoDyCo::UndirectedTree & undirected_tree,
                                const std::vector<std::string> & ft_names,
                                const std::vector<bool> & is_measure_direction_child_to_parent)
{
    SensorsTree sensors_tree;
    for(int i=0; i < ft_names.size(); i++ )
    {
        //Creating a new ft sensor to be added in the ft sensors structure
        KDL::CoDyCo::SixAxisForceTorqueSensor new_sens;


        if( undirected_tree.getJunction(ft_names[i]) != undirected_tree.getInvalidJunctionIterator() )
        {
            //Set the sensor name (for the time being equal to the junction name)
            new_sens.setName(ft_names[i]);
            //Set the junction name
            new_sens.setParent(ft_names[i]);
            int junction_index = undirected_tree.getJunction(ft_names[i])->getJunctionIndex();
            new_sens.setParentIndex(junction_index);
            KDL::CoDyCo::JunctionMap::const_iterator junct_it = undirected_tree.getJunction(ft_names[i]);

            int parent_index = junct_it->getParentLink()->getLinkIndex();
            int child_index = junct_it->getChildLink()->getLinkIndex();
            std::string parent_link_name = junct_it->getParentLink()->getName();
            std::string child_link_name = junct_it->getChildLink()->getName();

            if( is_measure_direction_child_to_parent[i] )
            {
                new_sens.setAppliedWrenchLink(parent_index);
            }
            else
            {
                new_sens.setAppliedWrenchLink(child_index);
            }

            // Currently we support only the case where the ft sensor frame is equal
            // to the child link frame
            new_sens.setSecondLinkSensorTransform(child_index,KDL::Frame::Identity());

            // Then, the parent_link_H_sensor transform is simply parent_link_H_child_link transform
            KDL::Frame parent_link_H_sensor = junct_it->pose(0.0,false);
            new_sens.setFirstLinkSensorTransform(parent_index,parent_link_H_sensor);

        }
        else
        {
            std::cerr << "[ERR] DynTree::generateSensorsTree: problem generating sensor for ft "
                      << ft_names[i] << std::endl;
            assert(false);
        }

        int ret = sensors_tree.addSensor(new_sens);

        assert(ret == i);
    }

    return sensors_tree;
}

typedef double T;

/*****************************************************************/
int main(int argc, char *argv[])
{
    yarp::os::Property opt;
    
    opt.fromCommand(argc,argv);
   
    if( !opt.check("dataset") ) {
        cout << "icubArmRegressor: " << endl;
        cout << "usage: ./icubArmRegressor --dataset dataset.csv --results results.csv --verbose --parameters estimatedParameters.csv" << endl;
        cout << "additional options: --n_samples the number of random regressors to generate for the numerical identifable subspace computation" << endl;
        cout << "additional options: --cut_freq the cut frequency to use for filtering torques and ft measures (in Hz, default 3 Hz)" << endl;
        cout << "additional options: --sampleTime the sample time of the data (in s, default: 0.01 s)" << endl;
        cout << "additional options: --lambda regularization parameter" << endl;
        cout << "additional options: --skip number of initial samples to skip" << endl;
        cout << "additional options: --verbose print debug information" << endl;
        cout << "additional options: --parameters param.csv output a file of trajectories of the estimated parameters" << endl;
        cout << "additional options: --results results.csv output a file of estimated outputs" << endl;
        return 0;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// Parameter parsing
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int n_samples;
    if(!opt.check("n_samples")) {
        n_samples = 1000;
    } 
    else 
    {
        n_samples = opt.find("n_samples").asInt();
    }
    
    double sampleTime;
    if(!opt.check("sampleTime")) {
        sampleTime = 0.01;
    } 
    else 
    {
        sampleTime = opt.find("sampleTime").asDouble();
    }
    
    double cut_freq;
    if(!opt.check("cut_freq")) {
        cut_freq = 3;
    } 
    else 
    {
        cut_freq = opt.find("cut_freq").asDouble();
    }
    
    double lambda;
    if(!opt.check("lambda")) 
    {
        lambda = 1.0; // Default lambda value   
    }
    else 
    {
        // Get lambda from parsed options
        lambda = opt.find("lambda").asDouble();
    }
    
    cout << "WARNING: Lambda is not optimized on the training data."<< endl;
    cout << "Selected lambda: " << lambda << endl;
    
    int skip;
    if(!opt.check("skip")) 
    {
        skip = 0;      // Default # of skipped samples
    } 
    else 
    {
        skip = opt.find("skip").asInt();
    }
    
    bool verbose;
    if(!opt.check("verbose")) 
    {
        verbose = false;
    } 
    else 
    {
        verbose = true;
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    /// Opening the dataset in csv format. 
    /// The format of the dataset is described in https://github.com/traversaro/dir-dataset-format/blob/master/dynamic_dataset_format.md 
    //////////////////////////////////////////////////////////////////////////////////

    std::string dataset_file_name(opt.find("dataset").asString());
    DynamicDatasetFile test_dataset;
    test_dataset.loadFromFile(dataset_file_name);

    ///////////////////////////////////////////////////////////////////////////////////
    //// Open a result file to dump output
    ///////////////////////////////////////////////////////////////////////////////////
    
    std::string results_file_name(opt.find("results").asString());
    std::ofstream results_file;
    results_file.open(results_file_name.c_str());
    
    //////////////////////////////////////////////////////////////
    ////////// Get file name for estimated parameters dumping
    //////////////////////////////////////////////////////////
    
    bool output_param = false;
    std::string params_file_name;
    
    if(  opt.check("parameters") ) {
        params_file_name = opt.find("parameters").asString().c_str();
        output_param = true;
    }
    
    std::ofstream param_file;
    
    if( output_param ) {
        param_file.open(params_file_name.c_str());
    }
    
    yarp::os::ResourceFinder rf;
    rf.setVerbose(true);
    rf.setDefault("name","icubArmRegressor");
    rf.setDefault("config","config.ini");
    rf.configure(argc,argv);

    if (rf.check("help"))
    {
        yInfo() << "Options:\n\n";
        yInfo() << "\t--urdf    file: name of the URDF file containing"
                   " the model of the iCub (default: model.urdf).  \n";
        return 0;
    }
        
    int c;

    string urdf_file_name;
    if(  opt.check("urdf") ) {
        urdf_file_name = opt.find("urdf").asString().c_str();
        cout << endl<< endl<< "urdf_file_name: " << urdf_file_name << endl<< endl;
//         cin>>c;
    }
//     cin>>c;
//     std::string urdf_filename = rf.findFile("urdf");
//     std::string urdf_filename = rf.findFile("model.urdf");

    yInfo() << "icubArmRegressor: Trying to open " << urdf_file_name
                                                  << " as iCub model";
                                                  
                                                  
    //////////////////////////////////////////////////
    /// Some data structures used to generate regressors
    //////////////////////////////////////////////////
    
    // Get KDL tree from the specified URDF file
    KDL::Tree icub_kdl_tree;
        
    if (!kdl_format_io::treeFromUrdfFile(urdf_file_name,icub_kdl_tree))
    {
        std::cerr << "Could not generate robot model and extract kdl tree" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "icub_kdl_tree information :" << endl
              << "Number of segments: " << icub_kdl_tree.getNrOfSegments() <<std::endl<<std::endl;
    
    // Get Root link of the overall robot
    std::string root_link_name = "root_link";

    
    /////// Define the names of the ft sensors
    //Can be l_arm_ft_sensor or r_arm_ft_sensor 
    std::string ft_sensor_name = "l_arm_ft_sensor";
    std::vector<std::string> ft_names;
    ft_names.push_back(ft_sensor_name);

    //Can be RIGHT_ELBOW_TORQUE_INDEX or LEFT_ELBOW_TORQUE_INDEX
    const int elbow_torque_global_id = LEFT_ELBOW_TORQUE_INDEX;

    //Can be RIGHT_ARM_FT_INDEX or LEFT_ARM_FT_INDEX
    const int ft_dataset_sensor_index = LEFT_ARM_FT_INDEX;
    
    //Can be right or left (gain between raw outputs)
    const double elbow_torque_gain = LEFT_ELBOW_GAIN; 

    //In the regressor generator, the ft is identified by an index
    //We add only a FT sensor for regressor, so it should be 0
    const int ft_regressor_generator_sensor_index = 0;

    /////// Define the names of the fake links
    ///// Some links that are always considered to have 0 mass and not considered in parametric estimation
    vector<string> fake_names;
    fake_names.push_back("torso");
    fake_names.push_back("imu_frame");
    fake_names.push_back("l_gripper");
    fake_names.push_back("r_gripper");    
    fake_names.push_back("l_sole");
    fake_names.push_back("r_sole");
    fake_names.push_back("l_wrist_1");
    fake_names.push_back("r_wrist_1");


    // Get KDL CoDyCo TreeSerialization
//     KDL::CoDyCo::TreeSerialization icub_serialization = icub_kdl_undirected_tree.getSerialization();
    
    std::vector<std::string> custom_serialization;
    
    custom_serialization.push_back("torso_pitch");
    custom_serialization.push_back("torso_roll");
    custom_serialization.push_back("torso_yaw");
    custom_serialization.push_back("neck_pitch");
    custom_serialization.push_back("neck_roll");
    custom_serialization.push_back("neck_yaw");
    custom_serialization.push_back("l_shoulder_pitch");
    custom_serialization.push_back("l_shoulder_roll");
    custom_serialization.push_back("l_shoulder_yaw");
    custom_serialization.push_back("l_elbow");
    custom_serialization.push_back("l_wrist_prosup");
    custom_serialization.push_back("l_wrist_pitch");
    custom_serialization.push_back("l_wrist_yaw");
    custom_serialization.push_back("r_shoulder_pitch");
    custom_serialization.push_back("r_shoulder_roll");
    custom_serialization.push_back("r_shoulder_yaw");
    custom_serialization.push_back("r_elbow");
    custom_serialization.push_back("r_wrist_prosup");
    custom_serialization.push_back("r_wrist_pitch");
    custom_serialization.push_back("r_wrist_yaw");
    
    // Following DOFs in random order, unused in this experiment
    custom_serialization.push_back("l_hip_pitch");
    custom_serialization.push_back("l_hip_roll");
    custom_serialization.push_back("l_hip_yaw");
    custom_serialization.push_back("r_hip_pitch");
    custom_serialization.push_back("r_hip_roll");
    custom_serialization.push_back("r_hip_yaw");
    custom_serialization.push_back("l_knee");
    custom_serialization.push_back("r_knee");
    custom_serialization.push_back("r_ankle_pitch");
    custom_serialization.push_back("r_ankle_roll");
    custom_serialization.push_back("l_ankle_pitch");
    custom_serialization.push_back("l_ankle_roll");
    
    // New tree serialization based on custom_serialization list
    KDL::CoDyCo::TreeSerialization icub_serialization(icub_kdl_tree, custom_serialization );
    
    cout << "icub_serialization number of links: " << icub_serialization.getNrOfLinks() << endl;
    cout << "icub_serialization: " << icub_serialization.toString() << endl;
    
    // Get KDL CoDyCo UndirectedTree
    KDL::CoDyCo::UndirectedTree icub_kdl_undirected_tree(icub_kdl_tree,icub_serialization);

    // Set dynamic base name
    // Can be l_upper_arm or r_upper_arm
    string subtree_dynamic_base = "l_upper_arm";
    
    // Compute dynamic_traversal
    Traversal dynamic_traversal;
    icub_kdl_undirected_tree.compute_traversal(dynamic_traversal,subtree_dynamic_base);

    // Set measure direction
    std::vector<bool> is_measure_direction_child_to_parent(ft_names.size(),false);
    SensorsTree sensors_tree = generateSensorsTree(icub_kdl_undirected_tree,ft_names,is_measure_direction_child_to_parent);
     
    // Set consider_ft_offset
    bool consider_ft_offset = true;
    
    /**
      * The dynamics regressor generator is a class for calculating arbitrary regressor
      * related to robot dynamics, for identification of dynamics parameters, such
      * as inertial parameters (masses, centers of mass, inertia tensor elements) or
      * other related parameters (for example force/torque sensor offset).
      */

    KDL::CoDyCo::Regressors::DynamicRegressorGenerator ft_regressor_generator(icub_kdl_undirected_tree,
                                                                              sensors_tree,
                                                                              root_link_name,
                                                                              consider_ft_offset,
                                                                              fake_names,
                                                                              verbose);
    // Change the base link to left or right arm
//     ft_regressor_generator.changeDynamicBase(subtree_dynamic_base);
//     ft_regressor_generator.changeKinematicBase(root_link_name);

    /////// Define the subtree we are interested into
    
    std::vector<std::string> l_arm_subtree, r_arm_subtree;
    l_arm_subtree.push_back("l_upper_arm");
    r_arm_subtree.push_back("r_upper_arm");
    
    //Can be r_arm_subtreee or l_arm_subtree
    std::vector<std::string> arm_subtree = l_arm_subtree;
    
    //Add one regressor to the generator:
    //    (1) 6 for ft sensors (with offset)
    int ret = ft_regressor_generator.addSubtreeRegressorRows(arm_subtree);
    assert(ret == 0);
    
    // Make pointers to F/T LPF and velocity and acceleration estimators
    iCub::ctrl::FirstOrderLowPassFilter *ft_filt;
    iCub::ctrl::AWLinEstimator *velocity_estimator;
    iCub::ctrl::AWPolyEstimator *acceleration_estimator;
    
    //Defining some variables used in the estimation loop
    KDL::JntArray q, dq, ddq;
    yarp::sig::Vector q_yarp, dq_yarp, ddq_yarp;
    
    q.resize(ft_regressor_generator.getNrOfDOFs());
    SetToZero(q);    
    dq = q;
    ddq = q;
    
    q_yarp.resize(q.rows(),0.0);
    dq_yarp = q_yarp;
    ddq_yarp = q_yarp;
    
    KDL::Wrench ft_sample;
    yarp::sig::Vector ft_sample_raw(6,0.0), ft_sample_filtered(6,0.0);
    
    //Defining gravity vector
//     const double g = 9.806;
    const double g = -9.806;
    KDL::Twist gravity(KDL::Vector(0.0,0.0,g) , KDL::Vector(0.0,0.0,0.0));
    
    cout << "ft_regressor_generator.getNrOfOutputs() [regressor rows] = " << ft_regressor_generator.getNrOfOutputs() << endl;
    cout << "ft_regressor_generator.getNrOfParameters() [regressor cols] = " << ft_regressor_generator.getNrOfParameters() << endl;
    
    //Defining the matrix that I will use for getting the regressor
    Eigen::MatrixXd ft_regressor(ft_regressor_generator.getNrOfOutputs(),ft_regressor_generator.getNrOfParameters());
    cout <<"rows of ft_regressor :" <<ft_regressor.rows()<<"; cols of ft_regressor :" << ft_regressor.cols()<<endl;
    
    // Create known terms vector
    Eigen::VectorXd ft_kt(ft_regressor_generator.getNrOfOutputs());
    cout <<"rows of ft_kt :" << ft_kt.rows()<<" cols of ft_kt :" << ft_kt.cols()<<endl;
    
    //Computing the identifiable subspace basis matrix
    Eigen::MatrixXd base_parameters_subspace; 
    int ret_value = 0;
    ret_value = ft_regressor_generator.computeNumericalIdentifiableSubspace(base_parameters_subspace);
    assert( ret_value == 0 );
    
    cout << endl <<"------------------------------------" << endl;
    cout << "Identifiable subspace computation" << endl;
    cout << "------------------------------------" << endl;
    cout << "n_samples: " << n_samples <<endl;
    cout << "ret_value: " << ret_value << endl;
    cout << "base_parameters_subspace rows: " << base_parameters_subspace.rows() << endl;
    cout << "base_parameters_subspace cols: " << base_parameters_subspace.cols() << endl;
    cout << "------------------------------------" << endl<< endl;
    
    //Defining the output standard deviation hyperparameter
    Eigen::VectorXd output_stddev(6);
    output_stddev << 0.0707119 ,  0.07460326,  0.11061799,  0.00253377,  0.00295331, 0.00281101;
   
    //Defining the estimator objects
    unsigned int nrOfBase_parameters = base_parameters_subspace.cols();
    multiTaskRecursiveLinearEstimator estimator_dynamic(nrOfBase_parameters,ft_regressor_generator.getNrOfOutputs(),lambda);
    estimator_dynamic.setOutputErrorStandardDeviation(output_stddev);
    
    // Get CAD parameters from model
    Eigen::VectorXd cad_parameters(ft_regressor_generator.getNrOfParameters());
    cad_parameters.setZero();
    KDL::CoDyCo::inertialParametersVectorLoopFakeLinks(icub_kdl_undirected_tree,cad_parameters,fake_names);
    
    // For comparing estimates between CAD and estimated parameters, we have to estimate
    // the offset for the CAD
    
    // Define objects for offset estimation
    multiTaskSVDLinearEstimator estimator_static_offset(6,6);
    int sample_nr;  // Current sample index
    int nmbr_of_samples_for_offset_calibration = 30;
    int nmbr_of_samples_for_offset_calibration_obtained = 0;
    Eigen::VectorXd offset = Eigen::VectorXd(6);
    Eigen::MatrixXd regressor_offset;
    Eigen::VectorXd offset_kt;
    
    cout << "test_dataset rows: " << test_dataset.getNrOfSamples() <<endl;
    
    
    
    
    
    /////////////////////////////
    //    Consistency Tests    //
    /////////////////////////////
    
    // Test 1

    cout << "------------------------------------------------------------------" <<endl;
    cout << " Home position test: q = [ 0 0 0 0 0 0 0 ]" << endl;
    cout << "------------------------------------------------------------------" <<endl;
    
    KDL::JntArray zero_q = q;
    SetToZero(zero_q);
        
    //Compute the regressors
    int rt = ft_regressor_generator.setRobotState(zero_q,zero_q,zero_q,gravity);
    cout << "rt: " << rt << endl;
    KDL::Wrench ft_sample_zero;
    rt = ft_regressor_generator.setFTSensorMeasurement(ft_regressor_generator_sensor_index,ft_sample_zero);
    cout << "rt: " << rt << endl;
    rt = ft_regressor_generator.computeRegressor(ft_regressor,ft_kt);
    cout << "rt: " << rt << endl;
    Eigen::Matrix<double,6,1> ft_cad_prediction =  ft_regressor * cad_parameters;

    cout << "ft_sample_zero: \nForce = " << endl << ft_sample_zero.force(0) << ft_sample_zero.force(1) << ft_sample_zero.force(2) << endl << "Torque = " << endl << ft_sample_zero.torque(0) << ft_sample_zero.torque(1) << ft_sample_zero.torque(2) << endl;
    cout << "ft_kt: "<< endl << ft_kt << endl;
    cout << "ft_cad_prediction: " << endl<< ft_cad_prediction << endl << endl;    
    cout << "Interesting regressor columns:" << endl;
    for(int k  = 0; k < 6; ++k)
    {
        for(int l  = 270; l < ft_regressor.cols(); l = l+10)
        {
            cout<< setw(20) << ft_regressor(k,l) /*"\t\t\t\t"*/;
        }
        cout << endl;
    }
    
    cout << "------------------------------------------------------------------" <<endl;
    
    

    
    // Test 2

    cout << "------------------------------------------------------------------" <<endl;
    cout << " Test 2: q = [ PI/2 0 0 0 0 0 0 ]" << endl;
    cout << "------------------------------------------------------------------" <<endl;
    
    KDL::JntArray q_t2 = zero_q;
    q_t2(0+6) = PI/2;
    
    //Compute the regressors
    rt= ft_regressor_generator.setRobotState(q_t2,zero_q,zero_q,gravity);
    cout << "rt: " << rt << endl;
    rt= ft_regressor_generator.setFTSensorMeasurement(ft_regressor_generator_sensor_index,ft_sample_zero);
    cout << "rt: " << rt << endl;
    rt= ft_regressor_generator.computeRegressor(ft_regressor,ft_kt);
    cout << "rt: " << rt << endl;
    ft_cad_prediction =  ft_regressor * cad_parameters;

    cout << "ft_sample_zero: \nForce = " << endl << ft_sample_zero.force(0) << ft_sample_zero.force(1) << ft_sample_zero.force(2) << endl << "Torque = " << endl << ft_sample_zero.torque(0) << ft_sample_zero.torque(1) << ft_sample_zero.torque(2) << endl;
    cout << "ft_kt: "<< endl << ft_kt << endl;
    cout << "ft_cad_prediction: " << endl<< ft_cad_prediction << endl;
    
    cout << "Interesting regressor columns:" << endl;
    for(int k  = 0; k < 6; ++k)
    {
        for(int l  = 270; l < ft_regressor.cols(); l = l+10)
        {
            cout<< setw(20) << ft_regressor(k,l) /*"\t\t\t\t"*/;
        }
        cout << endl;
    }
    
    cout << "------------------------------------------------------------------" <<endl;

    
    
    
    
    // Test 3
    
    cout << "------------------------------------------------------------------" <<endl;
    cout << " Test 3: q = [ 0 PI/2 0 0 0 0 0 ]" << endl;
    cout << "------------------------------------------------------------------" <<endl;

    KDL::JntArray q_t3 = zero_q;
    q_t3(1+6) = PI/2;
    
    //Compute the regressors
    ft_regressor_generator.setRobotState(q_t3,zero_q,zero_q,gravity);
    ft_regressor_generator.setFTSensorMeasurement(ft_regressor_generator_sensor_index,ft_sample_zero);
    ft_regressor_generator.computeRegressor(ft_regressor,ft_kt);
    ft_cad_prediction =  ft_regressor * cad_parameters;

    cout << "ft_sample_zero: \nForce = " << endl << ft_sample_zero.force(0) << ft_sample_zero.force(1) << ft_sample_zero.force(2) << endl << "Torque = " << endl << ft_sample_zero.torque(0) << ft_sample_zero.torque(1) << ft_sample_zero.torque(2) << endl;
    cout << "ft_kt: "<< endl << ft_kt << endl;
    cout << "ft_cad_prediction: " << endl<< ft_cad_prediction << endl;
    
    cout << "Interesting regressor columns:" << endl;
    for(int k  = 0; k < 6; ++k)
    {
        for(int l  = 270; l < ft_regressor.cols(); l = l+10)
        {
            cout<< setw(20) << ft_regressor(k,l) /*"\t\t\t\t"*/;
        }
        cout << endl;
    }
    
    cout << "------------------------------------------------------------------" <<endl;


        
    
    // Test 4
    cout << "------------------------------------------------------------------" <<endl;
    cout << " Test 4: q = [ 0 0 PI/2 0 0 0 0 ]" << endl;
    cout << "------------------------------------------------------------------" <<endl;
    
    KDL::JntArray q_t4 = zero_q;
    q_t4(2+6) = PI/2;
    
    //Compute the regressors
    ft_regressor_generator.setRobotState(q_t4,zero_q,zero_q,gravity);
    ft_regressor_generator.setFTSensorMeasurement(ft_regressor_generator_sensor_index,ft_sample_zero);
    ft_regressor_generator.computeRegressor(ft_regressor,ft_kt);
    ft_cad_prediction =  ft_regressor * cad_parameters;

    cout << "ft_sample_zero: \nForce = " << endl << ft_sample_zero.force(0) << ft_sample_zero.force(1) << ft_sample_zero.force(2) << endl << "Torque = " << endl << ft_sample_zero.torque(0) << ft_sample_zero.torque(1) << ft_sample_zero.torque(2) << endl;
    cout << "ft_kt: "<< endl << ft_kt << endl;
    cout << "ft_cad_prediction: " << endl<< ft_cad_prediction << endl;
    
    cout << "Interesting regressor columns:" << endl;
    for(int k  = 0; k < 6; ++k)
    {
        for(int l  = 270; l < ft_regressor.cols(); l = l+10)
        {
            cout<< setw(20) << ft_regressor(k,l) /*"\t\t\t\t"*/;
        }
        cout << endl;
    }
    
    cout << "------------------------------------------------------------------" <<endl;

    // Test 5
    
    cout << "------------------------------------------------------------------" <<endl;
    cout << " Test 5: q = [ 0 0 PI/2 0 0 0 0 ]" << endl;
    cout << "------------------------------------------------------------------" <<endl;
    KDL::JntArray q_t5 = zero_q;
    q_t5(2+6) = PI/2;
    Vector f_5(1,2,3), t_5(4,5,6);
    KDL::Wrench ft_sample_t5(f_5,t_5);

    //Compute the regressors
    ft_regressor_generator.setRobotState(q_t4,zero_q,zero_q,gravity);
    ft_regressor_generator.setFTSensorMeasurement(ft_regressor_generator_sensor_index,ft_sample_t5);
    ft_regressor_generator.computeRegressor(ft_regressor,ft_kt);
    ft_cad_prediction =  ft_regressor * cad_parameters;

    cout << "ft_sample_t5: \nForce = " << endl << ft_sample_t5.force(0) << ft_sample_t5.force(1) << ft_sample_t5.force(2) << endl << "Torque = " << endl << ft_sample_t5.torque(0) << ft_sample_t5.torque(1) << ft_sample_t5.torque(2) << endl;
    cout << "ft_kt: "<< endl << ft_kt << endl;
    cout << "ft_cad_prediction: " << endl<< ft_cad_prediction << endl;
        
    cout << "Interesting regressor columns:" << endl;
    for(int k  = 0; k < 6; ++k)
    {
        for(int l  = 270; l < ft_regressor.cols(); l = l+10)
        {
            cout<< setw(20) << ft_regressor(k,l) /*"\t\t\t\t"*/;
        }
        cout << endl;
    }
    
    cout << "------------------------------------------------------------------" <<endl;

    // Test 6
    cout << "------------------------------------------------------------------" <<endl;
    cout << " Test 6: q = [ 0 0 0 0 0 0 0 ]" << endl;
    cout << "------------------------------------------------------------------" <<endl;
    
    KDL::JntArray q_t6 = zero_q;

    //Compute the regressors
    ft_regressor_generator.setRobotState(q_t6,zero_q,zero_q,gravity);
    ft_regressor_generator.setFTSensorMeasurement(ft_regressor_generator_sensor_index,ft_sample_t5);
    ft_regressor_generator.computeRegressor(ft_regressor,ft_kt);
    ft_cad_prediction =  ft_regressor * cad_parameters;

    cout << "ft_sample_t5: \nForce = " << endl << ft_sample_t5.force(0) << ft_sample_t5.force(1) << ft_sample_t5.force(2) << endl << "Torque = " << endl << ft_sample_t5.torque(0) << ft_sample_t5.torque(1) << ft_sample_t5.torque(2) << endl;
    cout << "ft_kt: "<< endl << ft_kt << endl;
    cout << "ft_cad_prediction: " << endl<< ft_cad_prediction << endl;
    
    cout << "Interesting regressor columns:" << endl;
    for(int k  = 0; k < 6; ++k)
    {
        for(int l  = 270; l < ft_regressor.cols(); l = l+10)
        {
            cout<< setw(20) << ft_regressor(k,l) /*"\t\t\t\t"*/;
        }
        cout << endl;
    }
    
    cout << "------------------------------------------------------------------" <<endl;

    
//     int c;
//     cin >> c;
    
    
    
    
    
    /////////////////////
    // Estimation loop //
    /////////////////////
    
    for(sample_nr = skip; sample_nr < test_dataset.getNrOfSamples(); ++sample_nr ) 
    {
        if( verbose && sample_nr % 1000 == 0 )
            cout << "Looping on sample: " << sample_nr << endl;
        
        // Joints in dataset are : Torso (3), head(6), left arm (16), right_arm (16)
        // Depending on the dataset we use different parts, here we use only the arm, the torso and the head (which however should be still)
        DynamicSample sample;
        bool ret = test_dataset.getSample(sample_nr,sample);
        if( !ret ) return -1;
        assert(sample.getNrOfDOFs() == DATASET_DOFS);
        assert( sample.getJointPosition().rows() == DATASET_DOFS);
       
        convertJointsFromDatasetToModel(sample.getJointPosition(),q);
        
        //Get F/T reading
        ft_sample = sample.getWrenchMeasure(ft_dataset_sensor_index);
        
        // Print F/T reading
        if( verbose && sample_nr % 1000 == 0 ) {
            cout << "------------------------------------------------------------------" <<endl;
            cout << "FT sample from dataset" << endl;
            cout << "Force norm: "<<  ft_sample.force.Norm() << endl;
            cout << "Torque norm: "<<  ft_sample.torque.Norm() << endl;
            cout << "F: " << *ft_sample.force.data << "   " << *ft_sample.force.data+1 << "   " << *ft_sample.force.data+2 << endl << "T: " << *ft_sample.torque.data<< "   " << *ft_sample.torque.data+1 << "   " << *ft_sample.torque.data+2  << endl;
            cout << "------------------------------------------------------------------" <<endl;
        }
        
        // Filter F/T reading
        toYarp(ft_sample,ft_sample_raw);    // ft_sample_raw is in YARP format and has to be filtered
        
        //If the sample is the first one (skip) allocate the memory, with all the parameters
        if( sample_nr == skip ) 
        {
            // Instantiate LP filter and vel and acc estimators
            ft_filt = new iCub::ctrl::FirstOrderLowPassFilter(cut_freq,sampleTime,ft_sample_raw);
            velocity_estimator = new iCub::ctrl::AWLinEstimator(16,0.02);
            acceleration_estimator = new iCub::ctrl::AWQuadEstimator(50,0.02);
        } 
        else 
        {
            // Filter raw sample
            ft_sample_filtered = ft_filt->filt(ft_sample_raw);
            toKDL(ft_sample_filtered,ft_sample);    // update ft_sample in KDL format
        }
        
        // Print filtered F/T samples
        if( verbose && sample_nr % 1000 == 0 ) {
            cout << "------------------------------------------------------------------" <<endl;
            cout << "FT sample filtered" << endl;
            cout<< "F: " << *ft_sample.force.data << "   " << *ft_sample.force.data+1 << "   " << *ft_sample.force.data+2 << endl << "T: " << *ft_sample.torque.data<< "   " << *ft_sample.torque.data+1 << "   " << *ft_sample.torque.data+2  << endl;
            cout << "Torque norm: " << ft_sample.torque.Norm() << endl;
            cout << "Force norm: "  <<  ft_sample.force.Norm() << endl;
            cout << "------------------------------------------------------------------" <<endl;
        }
        
        
        
        // Velocities and accelerations estimation
        iCub::ctrl::AWPolyElement el;
        toYarp(q.data,el.data);
        el.time = sample.getTimestamp();
        
        // Estimate dq and ddq
        toEigen(velocity_estimator->estimate(el) , dq.data);
        toEigen(acceleration_estimator->estimate(el) , ddq.data);
        
        
        
        //Compute the regressors
        int rv =ft_regressor_generator.setRobotState(q,dq,ddq,gravity);
        rv = ft_regressor_generator.setFTSensorMeasurement(ft_regressor_generator_sensor_index,ft_sample);
        rv = ft_regressor_generator.computeRegressor(ft_regressor,ft_kt);
        
        // Print computed regressor info
        if( verbose && sample_nr % 1000 == 0) {
            cout << "------------------------------------------------------------------" <<endl;
            cout << "Interesting regressor columns:" << endl;
            for(int k  = 0; k < 6; ++k)
            {
                for(int l  = 270; l < ft_regressor.cols(); l = l+10)
                {
                    cout<< setw(20) << ft_regressor(k,l) /*"\t\t\t\t"*/;
                }
                cout << endl;
            }
//             cout << "Entire ft_regressor:" << endl<< ft_regressor << endl;
//             cout << "Entire ft_regressor * base_parameters_subspace:" << endl<< ft_regressor * base_parameters_subspace << endl;
            cout << "ft_regressor norm:" << endl<< ft_regressor.norm() << endl;
            cout << "ft_kt: "<< ft_kt << endl;
            cout << "Force norm: " << ft_kt.head(3).norm() << endl;
            cout << "------------------------------------------------------------------" <<endl;
        }
        
        
        
        //Compute offset for CAD parameters
        if( nmbr_of_samples_for_offset_calibration_obtained < nmbr_of_samples_for_offset_calibration ) {
            regressor_offset = ft_regressor.rightCols<6>();     
            offset_kt = ft_kt - ft_regressor * cad_parameters;
            
            estimator_static_offset.feedSample(regressor_offset,offset_kt);
            ++nmbr_of_samples_for_offset_calibration_obtained;
            
            if( nmbr_of_samples_for_offset_calibration_obtained == nmbr_of_samples_for_offset_calibration ) {
                estimator_static_offset.updateParameterEstimate();
                offset = estimator_static_offset.getParameterEstimate();
                offset.setZero();
                assert(offset.size() == 6);
                cad_parameters.tail<6>() = offset;
            }
        }
        
        
        
        // Predict output given the current estimate
        Eigen::Matrix<double,6,1> ft_estimated_prediction =  ft_regressor * base_parameters_subspace * estimator_dynamic.getParameterEstimate();
        Eigen::Matrix<double,6,1> ft_cad_prediction =  ft_regressor * cad_parameters;
        
        //Update the parameteres estimate
        estimator_dynamic.feedSampleAndUpdate(ft_regressor * base_parameters_subspace,ft_kt);
               
        
        
        // Add the F/T prediction to the results file
        for(int i=0; i < 6; ++i ) 
        {
            results_file << ft_estimated_prediction[i] << ",";
        }
        for(int i=0; i < 6; ++i ) 
        {
            results_file << ft_cad_prediction[i]; 
            if( i != 6 ) 
            {
                results_file << ",";
            }
        }
        results_file << endl;
      
        
        
        // Add the estimated parameters at this iteration to the output file
        if( output_param ) {
            Eigen::VectorXd params1 = estimator_dynamic.getParameterEstimate();
            for(int param_nr=0; param_nr < params1.size(); ++param_nr ) 
            {
                if( param_nr < params1.size() ) 
                {
                    param_file << params1[param_nr];
                } 
                if( param_nr != params1.size()-1 ) 
                {
                    param_file << ",";
                }
            }
            param_file << endl;
        }
    }
   

   // Close dump files
    results_file.close();
    
    if( output_param ) {
        param_file.close();
    }
    
    cout << "icubArmRegressor: saved results to " << results_file_name << endl;
    if( output_param ) {
        cout << "icubArmRegressor: saved parameters to " << params_file_name << endl;
    }
    cout << "icubArmRegressor: got " << sample_nr << " samples."  << endl;
    
    return 0;
}

