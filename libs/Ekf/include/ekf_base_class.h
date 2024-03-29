#ifndef EKFBASE_H
#define EKFBASE_H
#endif

#include<iostream>
#include<string>
#include<functional>
#include<cmath>

#include<Matrix/matrix_inv_class.h>
#include<constants.h>


using namespace std;

template <typename T>
class EkfBase{
	
	public:
		// constructors
		EkfBase(size_t num_states, size_t num_meas, size_t num_states_sensor, T sample_time_s, 
			MatrixInv<T> initial_state, MatrixInv<T> process_noise_q, MatrixInv<T> meas_noise_r, MatrixInv<T> initial_covariance_p, 
			bool compute_q_each_iter, float process_noise_eps);

		// destructor
		~EkfBase();

		void Run(const MatrixInv<T> &state_sensor_val, const MatrixInv<T> &meas_sensor_val, const bool meas_indices []);
		MatrixInv<T> GetCurrentState(){ return current_state_; }
		MatrixInv<T> GetStateJacobian(){ return state_jacobian_; }
		MatrixInv<T> GetSensorMeasurement(){ return computed_meas_; }
		MatrixInv<T> GetMeasComputedFromState(){ return meas_from_propogated_state_; }
		MatrixInv<T> GetComputedMeas(){ return computed_meas_; }
		MatrixInv<T> GetMeasJacobian(){ return meas_jacobian_; }
		MatrixInv<T> GetCovariance(){ return covariance_p_; }
	protected:

		// EKF dimension and sample time
		const size_t num_states_;
		const size_t num_meas_;
		const size_t num_states_sensor_;
		const T sample_time_s_;


		// member functions
		virtual void PropagateState(const MatrixInv<T> &state_sensor_val) = 0;
		virtual void ComputeStateJacobian(const MatrixInv<T> &state_sensor_val) = 0;
		virtual void ComputeStateNoiseJacobian(const MatrixInv<T> &previous_state) = 0;
		virtual void GetMeas(const MatrixInv<T> &meas_sensor_val) = 0;
		virtual void ComputeMeasJacobian(const MatrixInv<T> &meas_sensor_val) = 0;
		virtual void ComputeMeasNoiseJacobian(const MatrixInv<T> &meas_sensor_val) = 0;
		virtual void ComputeMeasFromState(size_t r_idx) = 0;
		virtual void ComputeControlToStateMap() = 0;

		void ComputeKalmanGainSequential(size_t r_idx);

		bool compute_q_each_iter_;

		// EKF variables
		MatrixInv<T> initial_state_;
		MatrixInv<T> current_state_;
		MatrixInv<T> computed_meas_;

		const MatrixInv<T> process_noise_q_;
		const MatrixInv<T> meas_noise_r_;
		T process_noise_eps_;
		MatrixInv<T> covariance_p_;

		// Variables to run EKF
		MatrixInv<T> time_propagated_state_;
		MatrixInv<T> meas_from_propogated_state_;
		MatrixInv<T> state_jacobian_;
		MatrixInv<T> state_noise_jacobian_;
		MatrixInv<T> meas_jacobian_;
		MatrixInv<T> meas_noise_jacobian_;
		MatrixInv<T> inv_part_of_kalman_gain_;
		MatrixInv<T> kalman_gain_;
		MatrixInv<T> kalman_gain_seq_;
		MatrixInv<T> map_controls_to_state_;
		MatrixInv<T> process_noise_eps_matrix_;

		MatrixInv<T> kalman_eye_;
};