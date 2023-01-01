#ifndef EKF16DOFQUAT_H
#define EKF16DOFQUAT_H


#include<iostream>
#include<string>
#include <functional>

#include <ekf_base_class.h>

using namespace std;

/* In State array indices correspond to the following
0 		-> q0
1:3 	-> [q1, q2, q3]
4:6 	-> [bwx, bwy, bwz] Gyro biases
7:9 	-> [X, Y, Z] NED Position
10:12 	-> [Vn, Ve, Vd] NED velocity
13:15 	-> [bax, bay, baz] Accel biases
*/
#define NUM_STATES 			16
/* In Meas array indices correspond to the following
0:2 	-> [mx, my, mz] unit mag vector in body frame
3:5 	-> [X, Y, Z] GPS Measured NED Position
6:8 	-> [Vn, Ve, Vd] GPS Measured NED velocity
*/
#define NUM_MEAS 			9

/* First three indices correspond to body angular rates and last 3 sensors corresponds to
accelerometers measurement without the gravity resolved in body frame*/
#define NUM_STATES_SENSOR 	6 

template <typename T>
class Ekf16DofQuat:public EkfBase<T>{
	public:
		// constructors
		Ekf16DofQuat(T sample_time_s, MatrixInv<T> initial_state, MatrixInv<T> process_noise_q, MatrixInv<T> meas_noise_r, MatrixInv<T> initial_covariance_p);

		// Run method specific to this EKF formulation
		void Run(MatrixInv<T> state_sensor_val, MatrixInv<T> meas_sensor_val, bool meas_indices[]);

		// Method to get Eulr angle from quaternion
		MatrixInv<T> GetEulerAngle();

		// destructor
		~Ekf16DofQuat();	

	protected:
		// member functions
		void PropagateState(MatrixInv<T> state_sensor_val);
		void ComputeStateJacobian(MatrixInv<T> state_sensor_val);
		void ComputeStateNoiseJacobian(MatrixInv<T> previous_state);
		void GetMeas(MatrixInv<T> meas_sensor_val);
		void ComputeMeasJacobian(MatrixInv<T> meas_sensor_val);
		void ComputeMeasNoiseJacobian(MatrixInv<T> meas_sensor_val);
		void ComputeMeasFromState(MatrixInv<T> time_propagated_state);

		MatrixInv<T> g_;			
};
#endif