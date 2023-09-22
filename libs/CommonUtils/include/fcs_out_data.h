//#############AUTOGENERATED##########
#include "fcsModel.h"


typedef struct fcsModel::ExtY_fcsModel_T FcsOutput;


struct DataFields{

	// Time it took for the EKF to run
	float dt_s;
	// GPS Time if available
	float time_of_week_ms;
	// Raw Sensor Data (Mag is calibrated and normalized)
	float imu_data[9];
	// GPS NED position and velocity if available;
	float ned_pos_m[3];
	float ned_vel_mps[3];
	// Current EKF State
	float ekf_current_state[16];
	// GPS Valid Flag
	float gps_valid_flag[6];

	// RC In
	float rc_in[7];

	// float actuator commands
	float actuator_cmds[4];

	//Autogenerated fields

	float fcsOut0;

	float fcsOut1;

	float fcsOut2;

	float fcsOut3;

	float fcsOut4;

	float fcsOut5;

	float fcsOut6;

	float fcsOut7;

	float fcsOut8;

	float fcsOut9;

	float fcsOut10;

	float fcsOut11;

	float fcsOut12;

	float fcsOut13;

	float fcsOut14;

	float fcsOut15;

	float fcsOut16;

	float fcsOut17;

	float fcsOut18;

	float fcsOut19;

	float fcsOut20;

	float fcsOut21;

	float fcsOut22;

	float fcsOut23;

	float fcsOut24;

};


void AssignFcsData(const FcsOutput &fcs_output, DataFields data_to_save[], size_t data_index);
