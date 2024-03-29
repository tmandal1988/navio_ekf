#ifndef EKF15DOF_H
#define EKF15DOF_H


#include<iostream>
#include<string>
#include <functional>

#include <ekf_base_class.h>

using namespace std;

#define NUM_STATES 			15
#define NUM_MEAS 			7

/* First three indices correspond to body angular rates and last 3 sensors corresponds to
accelerometers measurement without the gravity resolved in body frame*/
#define NUM_STATES_SENSOR 	6 

template <typename T>
class Ekf15Dof:public EkfBase<T>{
	public:
		// constructors
		Ekf15Dof(T sample_time_s, MatrixInv<T> initial_state, MatrixInv<T> process_noise_q, MatrixInv<T> meas_noise_r, MatrixInv<T> initial_covariance_p,
			bool compute_q_each_iter = false, float process_noise_eps = 1e-8);

		// destructor
		~Ekf15Dof();	

	protected:
		// member functions
		void PropagateState(const MatrixInv<T> &state_sensor_val);
		void ComputeStateJacobian(const MatrixInv<T> &state_sensor_val);
		void ComputeStateNoiseJacobian(const MatrixInv<T> &previous_state);
		void GetMeas(const MatrixInv<T> &meas_sensor_val);
		void ComputeMeasJacobian(const MatrixInv<T> &meas_sensor_val);
		void ComputeMeasNoiseJacobian(const MatrixInv<T> &meas_sensor_val);
		void ComputeMeasFromState(size_t r_idx);

		void ComputeControlToStateMap();
		void ComputeTrignometricValues();
		T WrapToPi(T ang_rad);

		// trignometric values
		T s_phi;
		T c_phi;

		T s_theta;
		T c_theta;
		T t_theta;
		T sc_theta;

		T s_psi;
		T c_psi;

		// Intermediate matrices
		MatrixInv<T> body_rates_2_euler_rates_;
		MatrixInv<T> c_b2ned_;


		MatrixInv<T> g;

			
};
#endif