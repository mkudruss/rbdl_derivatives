/*                                                                             *
 * TorqueMuscle 
 * Copyright (c) 2016 Matthew Millard <matthew.millard@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */


//==============================================================================
// INCLUDES
//==============================================================================



#include "../Millard2016TorqueMuscle.h"
#include "../csvtools.h"
#include "../../geometry/tests/numericalTestFunctions.h"
#include <UnitTest++.h>
#include <rbdl/rbdl_math.h>
#include <ctime>
#include <string>
#include <ostream>
#include <sstream>
#include <stdio.h>
#include <exception>
#include <cassert>
#include <vector>


using namespace RigidBodyDynamics::Addons::Muscle;
using namespace RigidBodyDynamics::Addons::Geometry;
using namespace std;
/*
   Constructor tests:
   1. Coefficients are copied over correctly.
   2. Curves are made correctly

   calcTorqueMuscleInfo test
   stiffness calculation
   power calculation

*/

TEST(ConstructorRegularCallCheck)
{


    //Check that the 3 constructors when called properly
    //do not abort
    Millard2016TorqueMuscle test0 = Millard2016TorqueMuscle();

    RigidBodyDynamics::Math::VectorNd c1c2c3c4c5c6(6);
    c1c2c3c4c5c6 << 0.161, 0.958,0.932,1.578,3.19,0.242;
    RigidBodyDynamics::Math::VectorNd b1k1b2k2(4);
    b1k1b2k2 << -1.21, -6.351, 0.476, 5.91;

    Millard2016TorqueMuscle test1 =
            Millard2016TorqueMuscle(c1c2c3c4c5c6,
                                     b1k1b2k2,
                                     1.75,
                                     75.0,
                                     0.0,
                                     1.0,
                                     1.0,
                                     "test_detailedConstructor");

    SubjectInformation subjectInfo;
      subjectInfo.gender          = GenderSet::Male;
      subjectInfo.ageGroup        = AgeGroupSet::Young18To25;
      subjectInfo.heightInMeters  =  1.732;
      subjectInfo.massInKg        = 69.0;


    Millard2016TorqueMuscle test2 =
    Millard2016TorqueMuscle(
          DataSet::Anderson2007,
          subjectInfo,
          Anderson2007::HipExtension,
          0.0,
          1.0,
          1.0,
          "test_easyConstructor");

    CHECK(abs( test2.getPassiveTorqueScale()-1.0) < TOL);

}



TEST(calcJointTorqueCorrectnessTests){


    double jointAngleOffset     = 0;    
    double signOfJointAngle     = 1;
    double signOfJointTorque    = 1;
    double err = 0.0;

    std::string name("test");

    SubjectInformation subjectInfo;
      subjectInfo.gender          = GenderSet::Male;
      subjectInfo.ageGroup        = AgeGroupSet::Young18To25;
      subjectInfo.heightInMeters  =  1.732;
      subjectInfo.massInKg        = 69.0;

    Millard2016TorqueMuscle tq =
            Millard2016TorqueMuscle(DataSet::Anderson2007,
                                     subjectInfo,
                                     Anderson2007::HipExtension,
                                     jointAngleOffset,
                                     signOfJointAngle,
                                     signOfJointTorque,
                                     name);
    //Zero out the passive forces so that calcMuscleTorque reports
    //just the active force - this allows us to test its correctness.
    tq.setPassiveTorqueScale(0.0);

    //Test that the get and set functions work for
    //maximum isometric torque
    double tauMaxOld = tq.getMaximumActiveIsometricTorque();
    double tauMax = tauMaxOld*10.0;
    tq.setMaximumActiveIsometricTorque(tauMax);
    CHECK(abs( tq.getMaximumActiveIsometricTorque()-tauMax)
          < TOL );

    RigidBodyDynamics::Math::VectorNd c1c2c3c4c5c6 =
            tq.getParametersC1C2C3C4C5C6();
    double thetaAtTauMax = c1c2c3c4c5c6[2];

    //TorqueMuscleInfo tqInfo;

    double torque = tq.calcJointTorque(thetaAtTauMax,0.0,1.0);
    err = abs(torque -tauMax);
    CHECK(err< TOL);

    double omegaAt75TauMax = c1c2c3c4c5c6[3];
    torque = tq.calcJointTorque(thetaAtTauMax,omegaAt75TauMax,1.0);
    err = abs(torque - 0.75*tauMax);
    CHECK( err < TOL);

    double omegaAt50TauMax = c1c2c3c4c5c6[4];
    torque = tq.calcJointTorque(thetaAtTauMax,omegaAt50TauMax,1.0);
    err = abs(torque -0.50*tauMax);
    CHECK(err < TOL);


}

TEST(calcTorqueMuscleInfoCorrectnessTests){

  double jointAngleOffset     = 0;
  double signOfJointAngle     = 1;
  double signOfJointTorque    = 1;
  double signOfJointVelocity  = signOfJointTorque;

  std::string name("test");

  SubjectInformation subjectInfo;
    subjectInfo.gender          = GenderSet::Female;
    subjectInfo.ageGroup        = AgeGroupSet::SeniorOver65;
    subjectInfo.heightInMeters  =  1.732;
    subjectInfo.massInKg        = 69.0;

  Millard2016TorqueMuscle tq =
          Millard2016TorqueMuscle(DataSet::Anderson2007,
                                   subjectInfo,
                                   Anderson2007::HipExtension,
                                   jointAngleOffset,
                                   signOfJointAngle,
                                   signOfJointTorque,
                                   name);
  double jointAngle       = 0.;
  double jointVelocity    = 0.;
  double activation       = 1.0;

  double tauMax = tq.getMaximumActiveIsometricTorque();
  RigidBodyDynamics::Math::VectorNd c1c2c3c4c5c6 =
          tq.getParametersC1C2C3C4C5C6();

  RigidBodyDynamics::Math::VectorNd b1k1b2k2 =
          tq.getParametersB1K1B2K2();

  int idx = -1;
  if(b1k1b2k2[0] > 0){
      idx = 0;
  }else if(b1k1b2k2[2] > 0){
      idx = 2;
  }
  double b = b1k1b2k2[idx];
  double k = b1k1b2k2[idx+1];
  double thetaAtPassiveTauMax = log(tauMax/b)/k;

  double thetaAtTauMax    = c1c2c3c4c5c6[2];
  double omegaAt75TauMax  = c1c2c3c4c5c6[3];
  double omegaAt50TauMax  = c1c2c3c4c5c6[4];

  double jointAngleAtTauMax =
          signOfJointAngle*thetaAtTauMax+jointAngleOffset;



  TorqueMuscleInfo tmi;
  tq.calcTorqueMuscleInfo(jointAngleAtTauMax,
                              0.,
                              activation,
                              tmi);

  //Keypoint check: active force components + fiber kinematics

  CHECK(abs(tmi.activation-1)                 < EPSILON);
  CHECK(abs(tmi.jointAngle-jointAngleAtTauMax)< TOL);
  CHECK(abs(tmi.jointAngularVelocity-0.)      < TOL);

  CHECK(abs(tmi.activation-1)                 < EPSILON);
  CHECK(abs(tmi.fiberAngle-thetaAtTauMax)     < TOL);
  CHECK(abs(tmi.fiberAngularVelocity-0.)      < TOL);

  CHECK(abs(tmi.fiberActiveTorque - tauMax)   < TOL);
  CHECK(abs(tmi.fiberActiveTorqueAngleMultiplier-1.0) < TOL);
  CHECK(abs(tmi.fiberActiveTorqueAngularVelocityMultiplier-1.0)<TOL);

  //Total force check
  double torque = tq.calcJointTorque(jointAngleAtTauMax,0,activation);
  double err    = abs(torque 
                      - signOfJointTorque*(
                          tmi.fiberActiveTorque+tmi.fiberPassiveTorque));
  CHECK(abs(torque 
            - signOfJointTorque*(
                tmi.fiberActiveTorque+tmi.fiberPassiveTorque)) < TOL);

  //Total active force scales with activation
  tq.calcTorqueMuscleInfo( jointAngleAtTauMax,
                           0.,
                           activation*0.5,
                           tmi);

  CHECK(abs(tmi.fiberActiveTorque - tauMax*0.5)   < TOL);

  //Keypoint check - at omega at 50% tau max
  tq.calcTorqueMuscleInfo(jointAngleAtTauMax,
                          signOfJointVelocity*omegaAt50TauMax,
                          activation,
                          tmi);

  CHECK(abs(tmi.jointAngularVelocity
            -signOfJointVelocity*omegaAt50TauMax)      < TOL);

  CHECK(abs(tmi.fiberAngularVelocity-omegaAt50TauMax)  < TOL);
  CHECK(abs(tmi.fiberActiveTorque - tauMax*0.5)   < TOL);

  //Keypoint check - power
  CHECK(abs(tmi.jointPower-tmi.fiberPower) < TOL);
  CHECK(abs(tmi.jointPower-tmi.fiberTorque*omegaAt50TauMax) < TOL);

  //Keypoint check - where passive curve hits tau max.
  double jointAngleAtPassiveTauMax =
          signOfJointAngle*thetaAtPassiveTauMax+jointAngleOffset;

  tq.calcTorqueMuscleInfo(jointAngleAtPassiveTauMax,
                          signOfJointVelocity*omegaAt50TauMax,
                          activation,
                          tmi);

  CHECK(abs(tmi.fiberPassiveTorqueAngleMultiplier-1) < TOL);

  //Numerically check the active and passive fiber stiffnesses
  double h = sqrt(EPSILON);
  tq.calcTorqueMuscleInfo(jointAngleAtTauMax*0.75,
                          0.,
                          activation,
                          tmi);


  TorqueMuscleInfo tmiL;
  TorqueMuscleInfo tmiR;

  tq.calcTorqueMuscleInfo(jointAngleAtTauMax*0.75-h,
                          0.,
                          activation,
                          tmiL);

  tq.calcTorqueMuscleInfo(jointAngleAtTauMax*0.75+h,
                          0.,
                          activation,
                          tmiR);

  double jointK = (tmiR.jointTorque-tmiL.jointTorque)/(2*h);
  err = tmi.jointStiffness - jointK;

  CHECK(abs(tmi.jointStiffness-jointK) < 1e-5);

  double fiberK = signOfJointAngle*(tmiR.fiberTorque-tmiL.fiberTorque)/(2*h);
  err = tmi.fiberStiffness - fiberK;
  CHECK(abs(tmi.fiberStiffness - fiberK) < 1e-5);

  tq.setPassiveTorqueScale(1.5);

  CHECK(abs(tq.getPassiveTorqueScale()-1.5)<TOL);

  tq.setPassiveTorqueScale(1.0);

  TorqueMuscleInfo tmi0;
  TorqueMuscleInfo tmi1;
  TorqueMuscleInfo tmi2;

  tq.calcTorqueMuscleInfo(jointAngleAtPassiveTauMax,
                          signOfJointVelocity*omegaAt50TauMax,
                          activation,
                          tmi0);

  tq.setPassiveTorqueScale(2.0);
  tq.calcTorqueMuscleInfo(jointAngleAtPassiveTauMax,
                          signOfJointVelocity*omegaAt50TauMax,
                          activation,
                          tmi1);

  CHECK(abs(tmi0.fiberPassiveTorqueAngleMultiplier -
            0.5*tmi1.fiberPassiveTorqueAngleMultiplier) < TOL);

  CHECK(abs(tmi0.fiberPassiveTorque -
            0.5*tmi1.fiberPassiveTorque) < TOL);

  double jtq = tq.calcJointTorque(jointAngleAtPassiveTauMax,
                                  signOfJointVelocity*omegaAt50TauMax,
                                  activation);
  err = jtq-tmi1.jointTorque;
  CHECK(abs(jtq-tmi1.jointTorque) < TOL );


  tq.calcTorqueMuscleInfo(jointAngleAtPassiveTauMax,
                          signOfJointVelocity*omegaAt50TauMax,
                          activation-SQRTEPSILON,
                          tmi0);

  tq.calcTorqueMuscleInfo(jointAngleAtPassiveTauMax,
                          signOfJointVelocity*omegaAt50TauMax,
                          activation,
                          tmi1);

  tq.calcTorqueMuscleInfo(jointAngleAtPassiveTauMax,
                          signOfJointVelocity*omegaAt50TauMax,
                          activation+SQRTEPSILON,
                          tmi2);

  double DtqDa = tmi1.DjointTorqueDactivation;
  double DtqDa_NUM = (tmi2.jointTorque-tmi0.jointTorque)/(2*SQRTEPSILON);
  err = abs(DtqDa-DtqDa_NUM);
  CHECK(abs(DtqDa-DtqDa_NUM) < abs(DtqDa)*1e-5 );


  tq.calcTorqueMuscleInfo(jointAngleAtPassiveTauMax-SQRTEPSILON,
                          signOfJointVelocity*omegaAt50TauMax,
                          activation,
                          tmi0);

  tq.calcTorqueMuscleInfo(jointAngleAtPassiveTauMax+SQRTEPSILON,
                          signOfJointVelocity*omegaAt50TauMax,
                          activation,
                          tmi2);

  double DtqDq = tmi1.DjointTorqueDjointAngle;
  double DtqDq_NUM = (tmi2.jointTorque-tmi0.jointTorque)/(2*SQRTEPSILON);
  err = abs(DtqDq-DtqDq_NUM);
  CHECK(abs(DtqDq-DtqDq_NUM) < abs(DtqDq)*1e-5 );

  tq.calcTorqueMuscleInfo(jointAngleAtPassiveTauMax,
                          signOfJointVelocity*omegaAt50TauMax-SQRTEPSILON,
                          activation,
                          tmi0);

  tq.calcTorqueMuscleInfo(jointAngleAtPassiveTauMax,
                          signOfJointVelocity*omegaAt50TauMax+SQRTEPSILON,
                          activation,
                          tmi2);

  double DtqDqdot = tmi1.DjointTorqueDjointAngularVelocity;
  double DtqDqdot_NUM = (tmi2.jointTorque-tmi0.jointTorque)/(2*SQRTEPSILON);
  err = abs(DtqDqdot-DtqDqdot_NUM);
  CHECK(abs(DtqDqdot-DtqDqdot_NUM) < abs(DtqDqdot)*1e-5 );

}

TEST(exampleUsage){


  bool printCurves = false;
  bool printAllCurves = false;

  //int dataSet = DataSetAnderson2007;

  //int gender  = 0; //male
  //int age     = 0; //young

  double jointAngleOffset     = 0;
  double signOfJointAngle     = 1;
  double signOfJointTorque    = 1;

  std::string name("test");

  SubjectInformation subjectInfo;
    subjectInfo.gender          = GenderSet::Male;
    subjectInfo.ageGroup        = AgeGroupSet::Young18To25;
    subjectInfo.heightInMeters  =  1.732;
    subjectInfo.massInKg        = 69.0;

  std::vector< Millard2016TorqueMuscle > muscleVector;

  bool exception = false;

  double angleTorqueSigns[][2] = {{-1, 1},
                                  {-1,-1},
                                  { 1,-1},
                                  { 1, 1},
                                  {-1, 1},
                                  {-1,-1}};

  Millard2016TorqueMuscle tqMuscle;
  std::stringstream tqName;
  int tqIdx;

  for(int i=0; i < Anderson2007::LastJointTorque; ++i){

      tqName.str("");
      tqName << DataSet.names[0]
             <<Anderson2007::JointTorqueNames[i];

      tqMuscle = Millard2016TorqueMuscle(
                    DataSet::Anderson2007,
                    subjectInfo,
                    Anderson2007::JointTorque(i),
                    0.0,
                    1.0,
                    1.0,
                    tqName.str() );


      if(printCurves)
          tqMuscle.printJointTorqueProfileToFile("",tqMuscle.getName() ,100);
  }

  for(int i=0; i < Gymnast::LastJointTorque; ++i){

      tqName.str("");
      tqName << DataSet.names[1]
             << Gymnast::JointTorqueNames[i];

      tqMuscle = Millard2016TorqueMuscle(
                    DataSet::Gymnast,
                    subjectInfo,
                    Gymnast::JointTorque(i),
                    0.0,
                    1.0,
                    1.0,
                    tqName.str() );


      if(printCurves)
          tqMuscle.printJointTorqueProfileToFile("",tqMuscle.getName() ,100);
  }

  tqIdx = -1;


  if(printAllCurves){
      std::stringstream muscleName;

      Millard2016TorqueMuscle muscle;
      int genderIdx,ageIdx,tqIdx;

      for(int age =0; age < Anderson2007::LastAgeGroup; ++age){
          for(int gender=0; gender < Anderson2007::LastGender; ++gender){
            for( int tqDir = 0; tqDir < Anderson2007::LastJointTorque; ++tqDir){
              //for(int joint=0; joint < 3; ++joint){
              //   for(int dir = 0; dir < 2; ++dir){
                      muscleName.str(std::string());

                      genderIdx = Anderson2007::Gender(gender);
                      ageIdx    = Anderson2007::AgeGroup(age);
                      tqIdx     = Anderson2007::JointTorque(tqDir);


                      muscleName  << "Anderson2007_"
                                  << AgeGroupSet::names[ageIdx]
                                  << "_"
                                  << GenderSet::names[genderIdx]
                                  << "_"
                                  << JointTorqueSet::names[tqIdx];

                      subjectInfo.ageGroup = AgeGroupSet::item(age);
                      subjectInfo.gender   = GenderSet::item(gender);

                       muscle = Millard2016TorqueMuscle(
                                 DataSet::Anderson2007,
                                 subjectInfo,
                                 Anderson2007::JointTorque(tqDir),
                                 0,
                                 1.0,
                                 1.0,
                                 muscleName.str());

                      const SmoothSegmentedFunction &tp
                          = muscle.getPassiveTorqueAngleCurve();
                      const SmoothSegmentedFunction &ta
                          = muscle.getActiveTorqueAngleCurve();
                      const SmoothSegmentedFunction &tv
                          = muscle.getTorqueAngularVelocityCurve();

                      RigidBodyDynamics::Math::VectorNd tpDomain
                          = tp.getCurveDomain();
                      RigidBodyDynamics::Math::VectorNd tvDomain
                          = tv.getCurveDomain();
                      RigidBodyDynamics::Math::VectorNd taDomain
                          = ta.getCurveDomain();



                      tp.printCurveToCSVFile(
                        "",
                        tp.getName(),
                        tpDomain[0]-0.1*(tpDomain[1]-tpDomain[0]),
                        tpDomain[1]+0.1*(tpDomain[1]-tpDomain[0]));
                      tv.printCurveToCSVFile(
                        "",
                        tv.getName(),
                        tvDomain[0]-0.1*(tvDomain[1]-tvDomain[0]),
                        tvDomain[1]+0.1*(tvDomain[1]-tvDomain[0]));
                      ta.printCurveToCSVFile(
                        "",
                        ta.getName(),
                        taDomain[0]-0.1*(taDomain[1]-taDomain[0]),
                        taDomain[1]+0.1*(taDomain[1]-taDomain[0]));

                  //}
              //}
            }
          }
      }


  }

   //catch(...){
   //     exceptionThrown = true;
   // }
   CHECK(true);




}

int main (int argc, char *argv[])
{
    return UnitTest::RunAllTests ();
}

