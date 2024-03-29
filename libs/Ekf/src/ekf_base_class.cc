#include "ekf_base_class.h"

template <typename T>
EkfBase<T>::EkfBase(size_t num_states, size_t num_meas, size_t num_states_sensor, T sample_time_s,
					MatrixInv<T> initial_state, MatrixInv<T> process_noise_q, MatrixInv<T> meas_noise_r, MatrixInv<T> initial_covariance_p,
					bool compute_q_each_iter, float process_noise_eps):
			 num_states_(num_states),
			 num_meas_(num_meas),
			 num_states_sensor_(num_states_sensor),
			 sample_time_s_(sample_time_s),
			 process_noise_q_(process_noise_q),
			 meas_noise_r_(meas_noise_r),
			 compute_q_each_iter_(compute_q_each_iter),
			 process_noise_eps_(process_noise_eps){	
		
		initial_state_ = initial_state;
		current_state_ = initial_state;
		computed_meas_ = MatrixInv<T> (num_meas_, 1);

		covariance_p_ = initial_covariance_p;

		state_jacobian_ = MatrixInv<T> (num_states_, num_states_);
		state_noise_jacobian_ = MatrixInv<T> (num_states_, num_states_, "eye");

		time_propagated_state_ = initial_state;
		meas_from_propogated_state_ = computed_meas_;
		meas_jacobian_ = MatrixInv<T> (num_meas_, num_states_);
		meas_noise_jacobian_ = MatrixInv<T> (num_meas_, num_meas_, "eye");

		kalman_eye_ = MatrixInv<T> (num_states_, num_states_, "eye");
		kalman_gain_seq_ = MatrixInv<T> (num_states_, 1);

		map_controls_to_state_ = MatrixInv<T> (num_states_, num_states_sensor_);
		process_noise_eps_matrix_ = MatrixInv<T> (num_states_, num_states_, "eye");
		process_noise_eps_matrix_ = process_noise_eps_matrix_*process_noise_eps_;
}

template <typename T>
EkfBase<T>::~EkfBase(){

}

template <typename T>
void EkfBase<T>::Run(const MatrixInv<T> &state_sensor_val, const MatrixInv<T> &meas_sensor_val, const bool meas_indices []){
	PropagateState(state_sensor_val);
	ComputeStateJacobian(state_sensor_val);

	GetMeas(meas_sensor_val);
	ComputeMeasJacobian(meas_sensor_val);

	//P = F*P*F' + L*Q*L';
	// if (compute_q_each_iter_){
	// 	ComputeControlToStateMap();
	// 	covariance_p_ = state_jacobian_*covariance_p_*state_jacobian_.Transpose() + map_controls_to_state_*process_noise_q_*map_controls_to_state_.Transpose() + process_noise_eps_matrix_;
	// }else{
	// 	covariance_p_ = state_jacobian_*covariance_p_*state_jacobian_.Transpose() + process_noise_q_;
	// }
	covariance_p_ = state_jacobian_*covariance_p_*state_jacobian_.Transpose() + process_noise_q_;

	current_state_ = time_propagated_state_;
	// sequentially update state with measurement
	for(size_t idx_r = 0; idx_r < num_meas_; idx_r++){
		if (meas_indices[idx_r]){
			ComputeMeasFromState(idx_r);
			ComputeKalmanGainSequential(idx_r);
			MatrixInv<T> meas_jacobian_row = meas_jacobian_.GetRow(idx_r);
			current_state_ = current_state_ + kalman_gain_seq_*( computed_meas_(idx_r) - meas_from_propogated_state_(idx_r) );
			covariance_p_ = covariance_p_ - kalman_gain_seq_*meas_jacobian_row*covariance_p_;
		}
	}
}

template <typename T>
inline void EkfBase<T>::ComputeKalmanGainSequential(size_t r_idx){
	kalman_gain_seq_ = covariance_p_*meas_jacobian_.Transpose().GetCol(r_idx);
	MatrixInv<T> k_denominator = meas_jacobian_.GetRow(r_idx)*kalman_gain_seq_ + meas_noise_r_(r_idx, r_idx);
	kalman_gain_seq_ = kalman_gain_seq_/k_denominator(0);
}

// Explicit template instantiation
template class EkfBase<float>;
template class EkfBase<double>;