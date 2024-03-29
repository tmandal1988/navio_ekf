#include "ekf_15dof_class.h"

template <typename T>
Ekf15Dof<T>::Ekf15Dof(T sample_time_s, MatrixInv<T> initial_state, MatrixInv<T> process_noise_q, MatrixInv<T> meas_noise_r, MatrixInv<T> initial_covariance_p,
					  bool compute_q_each_iter, float process_noise_eps): EkfBase<T>(NUM_STATES, NUM_MEAS, NUM_STATES_SENSOR, 
					  sample_time_s, initial_state, process_noise_q, meas_noise_r, initial_covariance_p, compute_q_each_iter, process_noise_eps){	

	this->meas_jacobian_(0, 2) = 1;
	this->meas_jacobian_(1, 6) = 1;
	this->meas_jacobian_(2, 7) = 1;
	this->meas_jacobian_(3, 8) = 1;
	this->meas_jacobian_(4, 9) = 1;
	this->meas_jacobian_(5, 10) = 1;
	this->meas_jacobian_(6, 11) = 1;


	T s_phi = 0;;
	T c_phi = 0;

	T s_theta = 0;
	T c_theta = 0;
	T t_theta = 0;
	T sc_theta = 0;

	T s_psi = 0;
	T c_psi = 0;

	g = MatrixInv<T>(3, 1);
	g(2) = 9.81;
}

template <typename T>
Ekf15Dof<T>::~Ekf15Dof(){

}

template <typename T>
inline void Ekf15Dof<T>::ComputeTrignometricValues(){
	// precalculate trignometric values
	s_phi = sin(this->current_state_(0));
	c_phi = cos(this->current_state_(0));

	s_theta = sin(this->current_state_(1));
	c_theta = cos(this->current_state_(1));
	t_theta = s_theta/c_theta; //tan(this->current_state_(1));
	sc_theta = 1/c_theta;

	s_psi = sin(this->current_state_(2));
	c_psi = cos(this->current_state_(2));

}
template <typename T>
inline T Ekf15Dof<T>::WrapToPi(T ang_rad){
	return ang_rad - 2*PI*floor( (ang_rad + PI) / PIx2 ); 
}

template <typename T>
void Ekf15Dof<T>::PropagateState(const MatrixInv<T> &state_sensor_val){
	// wrap the yaw angle	
	this->current_state_(2) = WrapToPi( this->current_state_(2) );

	ComputeTrignometricValues();
	T dt = this->sample_time_s_;

	this->time_propagated_state_(0) = this->current_state_(0) - dt*(this->current_state_(3) - state_sensor_val(0) + c_phi*t_theta*(this->current_state_(5) - state_sensor_val(2)) + s_phi*t_theta*(this->current_state_(4) - state_sensor_val(1)));

	this->time_propagated_state_(1) = this->current_state_(1) - dt*(c_phi*(this->current_state_(4) - state_sensor_val(1)) - s_phi*(this->current_state_(5) - state_sensor_val(2)));

	this->time_propagated_state_(2) = this->current_state_(2) - dt*((s_phi*(this->current_state_(4) - state_sensor_val(1)))/c_theta + (c_phi*(this->current_state_(5) - state_sensor_val(2)))/c_theta);

	this->time_propagated_state_(3) = this->current_state_(3);

	this->time_propagated_state_(4) = this->current_state_(4);

	this->time_propagated_state_(5) = this->current_state_(5);

	this->time_propagated_state_(6) = this->current_state_(6) + this->current_state_(9)*dt;

	this->time_propagated_state_(7) = this->current_state_(7) + this->current_state_(10)*dt;

	this->time_propagated_state_(8) = this->current_state_(8) + this->current_state_(11)*dt;

	this->time_propagated_state_(9) = this->current_state_(9) + dt*((s_phi*s_psi + c_phi*c_psi*s_theta)*(state_sensor_val(5) - this->current_state_(14)) - (c_phi*s_psi - c_psi*s_phi*s_theta)*(state_sensor_val(4) - this->current_state_(13)) + c_psi*c_theta*(state_sensor_val(3) - this->current_state_(12)));

	this->time_propagated_state_(10) = this->current_state_(10) + dt*((c_phi*c_psi + s_phi*s_psi*s_theta)*(state_sensor_val(4) - this->current_state_(13)) - (c_psi*s_phi - c_phi*s_psi*s_theta)*(state_sensor_val(5) - this->current_state_(14)) + c_theta*s_psi*(state_sensor_val(3) - this->current_state_(12)));

	this->time_propagated_state_(11) = this->current_state_(11) + dt*(g(2) - s_theta*(state_sensor_val(3) - this->current_state_(12)) + c_phi*c_theta*(state_sensor_val(5) - this->current_state_(14)) + c_theta*s_phi*(state_sensor_val(4) - this->current_state_(13)));

	this->time_propagated_state_(12) = this->current_state_(12);

	this->time_propagated_state_(13) = this->current_state_(13);

	this->time_propagated_state_(14) = this->current_state_(14);


}

template <typename T>
void Ekf15Dof<T>::ComputeControlToStateMap(){
	for(size_t idx_r = 0; idx_r < 3; idx_r++){
		for(size_t idx_c = 0; idx_c < 3; idx_c++){
			/* body_rates_2_euler_rates_ is calculated in PropagateState function
			*/
			this->map_controls_to_state_(idx_r, idx_c) = this->sample_time_s_*body_rates_2_euler_rates_(idx_r, idx_c);
		}
	}

	for(size_t idx_r = 0; idx_r < 3; idx_r++){
		for(size_t idx_c = 0; idx_c < 3; idx_c++){
			/* c_b2ned_ is calculated in PropagateState function
			*/
			this->map_controls_to_state_(idx_r + 9, idx_c + 3) = this->sample_time_s_*c_b2ned_(idx_r, idx_c);
		}
	}


}

template <typename T>
void Ekf15Dof<T>::ComputeStateJacobian(const MatrixInv<T> &state_sensor_val){
	T dt = this->sample_time_s_;

	// 1 row
	this->state_jacobian_(0 , 0) = 1 - dt*(c_phi*t_theta*(this->current_state_(4) - state_sensor_val(1)) - s_phi*t_theta*(this->current_state_(5) - state_sensor_val(2)));
	this->state_jacobian_(0 , 1) = -dt*(c_phi*(this->current_state_(5) - state_sensor_val(2))*sc_theta*sc_theta + s_phi*(this->current_state_(4) - state_sensor_val(1))*sc_theta*sc_theta);
	this->state_jacobian_(0 , 3) = -dt;
	this->state_jacobian_(0 , 4) = -dt*s_phi*t_theta;
	this->state_jacobian_(0 , 5) = -dt*c_phi*t_theta;

	// 2 row
	this->state_jacobian_(1 , 0) = dt*(c_phi*(this->current_state_(5) - state_sensor_val(2)) + s_phi*(this->current_state_(4) - state_sensor_val(1)));
	this->state_jacobian_(1 , 1) = 1;
	this->state_jacobian_(1 , 4) = -dt*c_phi;
	this->state_jacobian_(1 , 5) = dt*s_phi;

	// 3 row
	this->state_jacobian_(2 , 0) = dt*((s_phi*(this->current_state_(5) - state_sensor_val(2)))*sc_theta - (c_phi*(this->current_state_(4) - state_sensor_val(1)))*sc_theta);
	this->state_jacobian_(2 , 1) = -dt*((c_phi*t_theta*(this->current_state_(5) - state_sensor_val(2)))*sc_theta + (s_phi*t_theta*(this->current_state_(4) - state_sensor_val(1)))*sc_theta);
	this->state_jacobian_(2 , 2) = 1;
	this->state_jacobian_(2 , 4) = -(dt*s_phi)*sc_theta;
	this->state_jacobian_(2 , 5) = -(dt*c_phi)*sc_theta;

	// 4 row
	this->state_jacobian_(3 , 3) = 1;

	// 5 row
	this->state_jacobian_(4 , 4) = 1;

	// 6 row
	this->state_jacobian_(5 , 5) = 1;

	// 7 row
	this->state_jacobian_(6 , 6) = 1;
	this->state_jacobian_(6 , 9) = dt;

	// 8 row
	this->state_jacobian_(7 , 7) = 1;
	this->state_jacobian_(7 , 10) = dt;

	// 9 row
	this->state_jacobian_(8 , 8) = 1;
	this->state_jacobian_(8 , 11) = dt;

	// 10 row
	this->state_jacobian_(9 , 0) = dt*((s_phi*s_psi + c_phi*c_psi*s_theta)*(state_sensor_val(4) - this->current_state_(13)) + (c_phi*s_psi - c_psi*s_phi*s_theta)*(state_sensor_val(5) - this->current_state_(14)));
	this->state_jacobian_(9 , 1) = dt*(c_phi*c_psi*c_theta*(state_sensor_val(5) - this->current_state_(14)) - c_psi*s_theta*(state_sensor_val(3) - this->current_state_(12)) + c_psi*c_theta*s_phi*(state_sensor_val(4) - this->current_state_(13)));
	this->state_jacobian_(9 , 2) = -dt*((c_phi*c_psi + s_phi*s_psi*s_theta)*(state_sensor_val(4) - this->current_state_(13)) - (c_psi*s_phi - c_phi*s_psi*s_theta)*(state_sensor_val(5) - this->current_state_(14)) + c_theta*s_psi*(state_sensor_val(3) - this->current_state_(12)));
	this->state_jacobian_(9 , 9) = 1;
	this->state_jacobian_(9 , 12) = -dt*c_psi*c_theta;
	this->state_jacobian_(9 , 13) = dt*(c_phi*s_psi - c_psi*s_phi*s_theta);
	this->state_jacobian_(9 , 14) = -dt*(s_phi*s_psi + c_phi*c_psi*s_theta);

	// 11 row
	this->state_jacobian_(10 , 0) = -dt*((c_psi*s_phi - c_phi*s_psi*s_theta)*(state_sensor_val(4) - this->current_state_(13)) + (c_phi*c_psi + s_phi*s_psi*s_theta)*(state_sensor_val(5) - this->current_state_(14)));
	this->state_jacobian_(10 , 1) = dt*(c_phi*c_theta*s_psi*(state_sensor_val(5) - this->current_state_(14)) - s_psi*s_theta*(state_sensor_val(3) - this->current_state_(12)) + c_theta*s_phi*s_psi*(state_sensor_val(4) - this->current_state_(13)));
	this->state_jacobian_(10 , 2) = dt*((s_phi*s_psi + c_phi*c_psi*s_theta)*(state_sensor_val(5) - this->current_state_(14)) - (c_phi*s_psi - c_psi*s_phi*s_theta)*(state_sensor_val(4) - this->current_state_(13)) + c_psi*c_theta*(state_sensor_val(3) - this->current_state_(12)));
	this->state_jacobian_(10 , 10) = 1;
	this->state_jacobian_(10 , 12) = -dt*c_theta*s_psi;
	this->state_jacobian_(10 , 13) = -dt*(c_phi*c_psi + s_phi*s_psi*s_theta);
	this->state_jacobian_(10 , 14) = dt*(c_psi*s_phi - c_phi*s_psi*s_theta);

	// 12 row
	this->state_jacobian_(11 , 0) = dt*(c_phi*c_theta*(state_sensor_val(4) - this->current_state_(13)) - c_theta*s_phi*(state_sensor_val(5) - this->current_state_(14)));
	this->state_jacobian_(11 , 1) = -dt*(c_theta*(state_sensor_val(3) - this->current_state_(12)) + c_phi*s_theta*(state_sensor_val(5) - this->current_state_(14)) + s_phi*s_theta*(state_sensor_val(4) - this->current_state_(13)));
	this->state_jacobian_(11 , 11) = 1;
	this->state_jacobian_(11 , 12) = dt*s_theta;
	this->state_jacobian_(11 , 13) = -dt*c_theta*s_phi;
	this->state_jacobian_(11 , 14) = -dt*c_phi*c_theta;

	// 13 row
	this->state_jacobian_(12 , 12) = 1;

	// 14 row
	this->state_jacobian_(13 , 13) = 1;

	// 15 row
	this->state_jacobian_(14 , 14) = 1;
}

template <typename T>
void Ekf15Dof<T>::ComputeStateNoiseJacobian(const MatrixInv<T> &previous_state){
	// Not used right now
	//1st row
	this->state_noise_jacobian_(0, 15) = this->sample_time_s_;
	this->state_noise_jacobian_(0, 16) = this->sample_time_s_*s_phi*t_theta;
	this->state_noise_jacobian_(0, 17) = this->sample_time_s_*c_phi*t_theta;

	//2nd row
	this->state_noise_jacobian_(1, 16) = this->sample_time_s_*c_phi;
	this->state_noise_jacobian_(1, 17) = -this->sample_time_s_*s_phi;

	//3rd row
	this->state_noise_jacobian_(2, 16) = this->sample_time_s_*s_phi*sc_theta;
	this->state_noise_jacobian_(2, 17) = -this->sample_time_s_*c_phi*sc_theta;

	//4th row
	this->state_noise_jacobian_(3, 3) = 1;
	//5th row
	this->state_noise_jacobian_(4, 4) = 1;
	//6th row
	this->state_noise_jacobian_(5, 5) = 1;

	//7th row
	this->state_noise_jacobian_(6, 18) = pow(this->sample_time_s_, 2)*c_theta*c_psi;
	this->state_noise_jacobian_(6, 19) = pow(this->sample_time_s_, 2)*c_theta*s_psi;
	this->state_noise_jacobian_(6, 20) = -pow(this->sample_time_s_, 2)*s_theta;

	//8th row
	this->state_noise_jacobian_(7, 18) = pow(this->sample_time_s_, 2)*(s_phi*s_theta*c_psi - c_phi*s_psi);
	this->state_noise_jacobian_(7, 19) = pow(this->sample_time_s_, 2)*(s_phi*s_theta*s_psi + c_phi*c_psi);
	this->state_noise_jacobian_(7, 20) = pow(this->sample_time_s_, 2)*c_theta*s_phi;

	//9th row
	this->state_noise_jacobian_(8, 18) = pow(this->sample_time_s_, 2)*(c_phi*s_theta*c_psi + s_phi*s_psi);
	this->state_noise_jacobian_(8, 19) = pow(this->sample_time_s_, 2)*(c_phi*s_theta*s_psi - s_phi*c_psi);
	this->state_noise_jacobian_(8, 20) = pow(this->sample_time_s_, 2)*c_theta*c_phi;

	//10th row
	this->state_noise_jacobian_(9, 18) = this->sample_time_s_*c_theta*c_psi;
	this->state_noise_jacobian_(9, 19) = this->sample_time_s_*c_theta*s_psi;
	this->state_noise_jacobian_(9, 20) = -this->sample_time_s_*s_theta;

	//11th row
	this->state_noise_jacobian_(10, 18) = this->sample_time_s_*(s_phi*s_theta*c_psi - c_phi*s_psi);
	this->state_noise_jacobian_(10, 19) = this->sample_time_s_*(s_phi*s_theta*s_psi + c_phi*c_psi);
	this->state_noise_jacobian_(10, 20) = this->sample_time_s_*c_theta*s_phi;

	//12th row
	this->state_noise_jacobian_(11, 18) = this->sample_time_s_*(c_phi*s_theta*c_psi + s_phi*s_psi);
	this->state_noise_jacobian_(11, 19) = this->sample_time_s_*(c_phi*s_theta*s_psi - s_phi*c_psi);
	this->state_noise_jacobian_(11, 20) = this->sample_time_s_*c_theta*c_phi;

	//13th row
	this->state_noise_jacobian_(12, 12) = 1;
	//14th row
	this->state_noise_jacobian_(13, 13) = 1;
	//15th row
	this->state_noise_jacobian_(14, 14) = 1;
}

template <typename T>
void Ekf15Dof<T>::GetMeas(const MatrixInv<T> &meas_sensor_val){
	/* First 3 indices should be magnetometer x, y, z readings in body frame,
	next 3 indices should be GPS measured position translated into NED frame and next
	3 indices should be GPS measured NED velocity*/

	//Compute heading from magnetometer
	T mx = meas_sensor_val(0);
	T my = meas_sensor_val(1);
	T mz = meas_sensor_val(2);

	// MatrixInv<T> mag2d_projection = { {c_theta, s_theta*s_phi, s_theta*c_phi},
	// 								  {   0   ,     c_phi    ,   -s_phi     } };

	// MatrixInv<T> mag_vector = { {m1}, {m2}, {m3} };


	//Assume we are in bay area and use a declination 
	 T magnetic_declination = 13.01*DEG2RAD;

	// MatrixInv<T> mag2d = mag2d_projection*mag_vector;
	// this->computed_meas_(0) = atan2(-mag2d(1), mag2d(0)) + magnetic_declination;
	this->computed_meas_(0) = atan2(-my*c_phi + mz*s_phi, mx*c_theta + (my*s_phi + mz*c_phi)*s_theta) + magnetic_declination;

	while (this->computed_meas_(0) > this->current_state_(2) + PI){
			this->computed_meas_(0) = this->computed_meas_(0) - PIx2;
	}
	while (this->computed_meas_(0) < this->current_state_(2) - PI){
			this->computed_meas_(0) = this->computed_meas_(0) + PIx2;
	}

	this->computed_meas_(1) = meas_sensor_val(3);//Px
	this->computed_meas_(2) = meas_sensor_val(4);//Py
	this->computed_meas_(3) = meas_sensor_val(5);//Pz
	this->computed_meas_(4) = meas_sensor_val(6);//Vx
	this->computed_meas_(5) = meas_sensor_val(7);//Vy
	this->computed_meas_(6) = meas_sensor_val(8);//Vz
}

template <typename T>
void Ekf15Dof<T>::ComputeMeasJacobian(const MatrixInv<T> &meas_sensor_val){
	// Do nothing
}

template <typename T>
void Ekf15Dof<T>::ComputeMeasNoiseJacobian(const MatrixInv<T> &meas_sensor_val){
	// Do Nothing
}

template <typename T>
inline void Ekf15Dof<T>::ComputeMeasFromState(size_t r_idx){
	// In this case it is y = H*x
	this->meas_from_propogated_state_(r_idx) = (this->meas_jacobian_.GetRow(r_idx)*this->current_state_)(0);
}

// Explicit template instantiation
template class Ekf15Dof<float>;
template class Ekf15Dof<double>;