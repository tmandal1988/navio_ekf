#include "ekf_16dof_quat_class.h"

template <typename T>
Ekf16DofQuat<T>::Ekf16DofQuat(T sample_time_s, MatrixInv<T> initial_state, MatrixInv<T> process_noise_q, MatrixInv<T> meas_noise_r, MatrixInv<T> initial_covariance_p):
			EkfBase<T>(NUM_STATES, NUM_MEAS, NUM_STATES_SENSOR, sample_time_s, initial_state, process_noise_q, meas_noise_r, initial_covariance_p){	

	g_ = MatrixInv<T>(3, 1);
	g_(2) = 9.81;
}

template <typename T>
Ekf16DofQuat<T>::~Ekf16DofQuat(){

}

template <typename T>
void Ekf16DofQuat<T>::PropagateState(MatrixInv<T> state_sensor_val){
	// wrap the yaw angle	
	//this->current_state_(2) = remainder(this->current_state_(2) + PI, PIx2) - PI;

	//ComputeTrignometricValues();
	T dt = this->sample_time_s_;

	// Easy to use names
	T q0 = this->current_state_(0);
	T q1 = this->current_state_(1);
	T q2 = this->current_state_(2);
	T q3 = this->current_state_(3);


	MatrixInv<T> propogated_quat(4, 1);
	MatrixInv<T> propogated_velocity(3, 1);

	MatrixInv<T> body_rates_2_quat_rates = { {-q1, -q2, -q3}, {q0, -q3, q2}, {q3, q0, -q1}, {-q2, q1, q0} };

	MatrixInv<T> c_b2ned = { { static_cast<T>( 1 - 2.0*( pow(q2, 2) + pow(q3, 2) ) ), 2*( q1*q2 - q3*q0 ), 2*( q1*q3 + q2*q0 ) },
								   { 2*( q1*q2 + q3*q0 ), static_cast<T>( 1 - 2.0*( pow(q1, 2) + pow(q3, 2) ) ), 2*( q2*q3 - q1*q0 ) },
								   {2*( q1*q3 - q2*q0 ), 2*( q2*q3 + q1*q0 ), static_cast<T>( 1 - 2.0*( pow(q1, 2) + pow(q2, 2) ) ) } };

	MatrixInv<T> angular_rate_vector = { {state_sensor_val(0) - this->current_state_(4)},
										{state_sensor_val(1) - this->current_state_(5)},
										{state_sensor_val(2) - this->current_state_(6)} };
	MatrixInv<T> accel_vector = { {state_sensor_val(3) - this->current_state_(13)},
								  {state_sensor_val(4) - this->current_state_(14)},
								  {state_sensor_val(5) - this->current_state_(15)} };


	propogated_quat = (body_rates_2_quat_rates*angular_rate_vector);
	//propogated_quat = 0.5*dt*(body_rates_2_quat_rates*angular_rate_vector);

	propogated_velocity = (c_b2ned*accel_vector + g_)*dt;

	// propagate quaternions
	this->time_propagated_state_(0) = this->current_state_(0) + propogated_quat(0);
	this->time_propagated_state_(1) = this->current_state_(1) + propogated_quat(1);
	this->time_propagated_state_(2) = this->current_state_(2) + propogated_quat(2);
	this->time_propagated_state_(3) = this->current_state_(3) + propogated_quat(3);

	// if ( abs(this->time_propagated_state_(2)) > 180/57.29578 ){
	// 	this->time_propagated_state_(2) = -( this->time_propagated_state_(2)/abs(this->time_propagated_state_(2)) )*( 360/57.29578 - abs(this->time_propagated_state_(2)) ) ;
	// }
	
	// propagate gyro bias
	this->time_propagated_state_(4) = this->current_state_(4);
	this->time_propagated_state_(5) = this->current_state_(5);
	this->time_propagated_state_(6) = this->current_state_(6);

	// propagate NED X position
	this->time_propagated_state_(7) = this->current_state_(7) + dt*this->current_state_(10);
	// propagate NED Y position
	this->time_propagated_state_(8) = this->current_state_(8) + dt*this->current_state_(11);
	// propagate NED Z position
	this->time_propagated_state_(9) = this->current_state_(9) + dt*this->current_state_(12);

	// propagate NED X velocity
	this->time_propagated_state_(10) = this->current_state_(10) + propogated_velocity(0);
	this->time_propagated_state_(11) = this->current_state_(11) + propogated_velocity(1);
	this->time_propagated_state_(12) = this->current_state_(12) + propogated_velocity(2);

	// // propagate accel bias
	this->time_propagated_state_(13) = this->current_state_(13);
	this->time_propagated_state_(14) = this->current_state_(14);
	this->time_propagated_state_(15) = this->current_state_(15);
}

template <typename T>
void Ekf16DofQuat<T>::ComputeStateJacobian(MatrixInv<T> state_sensor_val){
	T dt = this->sample_time_s_;

	// Easy to use names
	T q0 	= this->current_state_(0);
	T q1 	= this->current_state_(1);
	T q2 	= this->current_state_(2);
	T q3 	= this->current_state_(3);

	T wx 	= state_sensor_val(0);
	T wy 	= state_sensor_val(1);
	T wz 	= state_sensor_val(2);

	T bwx 	= this->current_state_(4);
	T bwy	= this->current_state_(5);
	T bwz	= this->current_state_(6);

	T ax 	= state_sensor_val(3);
	T ay 	= state_sensor_val(4);
	T az	= state_sensor_val(5);

	T bax 	= this->current_state_(13);
	T bay 	= this->current_state_(14);
	T baz 	= this->current_state_(15);

	// 1st row
	this->state_jacobian_(0, 0) = 1;
	this->state_jacobian_(0, 1) = dt*(bwx - wx)*0.5;	
	this->state_jacobian_(0, 2) = dt*(bwy - wy)*0.5;
	this->state_jacobian_(0, 3) = dt*(bwz - wz)*0.5;

	this->state_jacobian_(0, 4) = q1*0.5;
	this->state_jacobian_(0, 5) = q2*0.5;
	this->state_jacobian_(0, 6) = q3*0.5;


	//2nd row
	this->state_jacobian_(1, 0) = dt*(wx - bwx)*0.5;
	this->state_jacobian_(1, 1) = 1;
	this->state_jacobian_(1, 2) = dt*(wz - bwz)*0.5;
	this->state_jacobian_(1, 3) = dt*(bwy - wy)*0.5;

	this->state_jacobian_(1, 4) = -q0*0.5;
	this->state_jacobian_(1, 5) = q3*0.5;
	this->state_jacobian_(1, 6) = -q2*0.5;

	//3rd row
	this->state_jacobian_(2, 0) = dt*(wy - bwy)*0.5;
	this->state_jacobian_(2, 1) = dt*(wz - bwz)*0.5;
	this->state_jacobian_(2, 2) = 1;
	this->state_jacobian_(2, 3) = dt*(wx - bwx)*0.5;

	this->state_jacobian_(2, 4) = -q3*0.5;
	this->state_jacobian_(2, 5) = -q0*0.5;
	this->state_jacobian_(2, 6) = q1*0.5;

	//4th row
	this->state_jacobian_(3, 0) = dt*(wz - bwz)*0.5;
	this->state_jacobian_(3, 1) = dt*(wy - bwy)*0.5;
	this->state_jacobian_(3, 2) = dt*(wx - bwx)*0.5;
	this->state_jacobian_(3, 3) = 1;

	this->state_jacobian_(3, 4) = q2*0.5;
	this->state_jacobian_(3, 5) = -q1*0.5;
	this->state_jacobian_(3, 6) = -q0*0.5;

	//5th row
	this->state_jacobian_(4, 4) = 1;
	//6th row
	this->state_jacobian_(5, 5) = 1;
	//7th row
	this->state_jacobian_(6, 6) = 1;

	//8th row
	this->state_jacobian_(7, 7) = 1;
	this->state_jacobian_(7, 10) = this->sample_time_s_;

	//9th row
	this->state_jacobian_(8, 8) = 1;
	this->state_jacobian_(8, 11) = this->sample_time_s_;

	//10th row
	this->state_jacobian_(9, 9) = 1;
	this->state_jacobian_(9, 12) = this->sample_time_s_;

	//11th row
	this->state_jacobian_(10, 0) = dt*( 2*q3*(bay - ay) - 2*q2*(baz - az) );
	this->state_jacobian_(10, 1) = dt*( -2*q2*(bay - ay) - 2*q3*(baz - az) );
	this->state_jacobian_(10, 2) = dt*( 4*q2*(bax - ax) - 2*q1*(bay - ay) - 2*q0*(baz - az) );
	this->state_jacobian_(10, 3) = dt*( 2*q0*(bay - ay) + 4*q3*(bax - ax) - 2*q1*(baz - az) );

	this->state_jacobian_(10, 10) = 1;

	this->state_jacobian_(10, 13) = dt*(2*pow(q2, 2) + 2*pow(q3, 2) - 1);
	this->state_jacobian_(10, 14) = dt*(2*q0*q3 - 2*q1*q2);
	this->state_jacobian_(10, 15) = dt*(-2*q0*q2 - 2*q1*q3);


	//12th row
	this->state_jacobian_(11, 0) = dt*( 2*q1*(baz - az) - 2*q3*(bax - ax) );
	this->state_jacobian_(11, 1) = dt*( 4*q1*(bay - ay) - 2*q2*(bax - ax) + 2*q0*(baz - az) );
	this->state_jacobian_(11, 2) = dt*( -2*q1*(bax - ax) - 2*q3*(baz - az) );
	this->state_jacobian_(11, 3) = dt*( 4*q3*(bay - ay) - 2*q0*(bax - ax) - 2*q2*(baz - az) );

	this->state_jacobian_(11, 11) = 1;

	this->state_jacobian_(11, 13) = dt*(-2*q0*q3 - 2*q1*q2);
	this->state_jacobian_(11, 14) = dt*(2*pow(q1, 2) + 2*pow(q3, 2) - 1);
	this->state_jacobian_(11, 15) = dt*(2*q0*q1 - 2*q2*q3);

	//13th row
	this->state_jacobian_(12, 0) = dt*( 2*q2*(bax - ax) - 2*q1*(bay - ay) );
	this->state_jacobian_(12, 1) = dt*( 4*q1*(baz - az) - 2*q3*(bax - ax) - 2*q0*(bay - ay) );
	this->state_jacobian_(12, 2) = dt*( 2*q0*(bax - ax) - 2*q3*(bay - ay) + 4*q2*(baz - az) );
	this->state_jacobian_(12, 3) = dt*( -2*q1*(bax - ax) - 2*q2*(bay - ay) );

	this->state_jacobian_(12, 12) = 1;

	this->state_jacobian_(12, 13) = dt*( 2*q0*q2 - 2*q1*q3 );
	this->state_jacobian_(12, 14) = dt*( -2*q0*q1 - 2*q2*q3 );
	this->state_jacobian_(12, 15) = dt*( 2*pow(q1, 2) + 2*pow(q2, 2) - 1 );


	//14th row
	this->state_jacobian_(13, 13) = 1;
	//15th row
	this->state_jacobian_(14, 14) = 1;
	//16th row
	this->state_jacobian_(15, 15) = 1;
}

template <typename T>
void Ekf16DofQuat<T>::ComputeStateNoiseJacobian(MatrixInv<T> previous_state){
	// Do Nothing
}

template <typename T>
void Ekf16DofQuat<T>::GetMeas(MatrixInv<T> meas_sensor_val){
	/* First 3 indices should be normalized magnetometer x, y, z readings in body frame,
	next 3 indices should be GPS measured position translated into NED frame and last
	3 indices should be GPS measured NED velocity*/

	for(size_t idx = 0; idx < 9; idx++)
		this->computed_meas_(idx) = meas_sensor_val(idx);
	
}

template <typename T>
void Ekf16DofQuat<T>::ComputeMeasJacobian(MatrixInv<T> meas_sensor_val){
	// Easy to use names
	T q0 	= this->time_propagated_state_(0);
	T q1 	= this->time_propagated_state_(1);
	T q2 	= this->time_propagated_state_(2);
	T q3 	= this->time_propagated_state_(3);

	T magnetic_declination = 13.01*DEG2RAD; // this can be an input later
	T sin_md = sin(magnetic_declination);
	T cos_md = cos(magnetic_declination);

	// 1st row	
	this->meas_jacobian_(0, 0) = 2*q3*sin_md;
	this->meas_jacobian_(0, 1) = 2*q2*sin_md;
	this->meas_jacobian_(0, 2) = 2*q1*sin_md  - 4*q2*cos_md;
	this->meas_jacobian_(0, 3) = 2*q0*sin_md  - 4*q3*cos_md;

	//2nd row
	this->meas_jacobian_(1, 0) = -2*q3*cos_md;
	this->meas_jacobian_(1, 1) = 2*q2*cos_md - 4*q1*sin_md;
	this->meas_jacobian_(1, 2) = 2*q1*cos_md;
	this->meas_jacobian_(1, 3) = -2*q0*cos_md - 4*q3*sin_md;

	//3rd row
	this->meas_jacobian_(2, 0) = 2*q2*cos_md - 2*q1*sin_md;
	this->meas_jacobian_(2, 1) = 2*q3*cos_md - 2*q0*sin_md;
	this->meas_jacobian_(2, 2) = 2*q0*cos_md + 2*q3*sin_md;
	this->meas_jacobian_(2, 3) = 2*q1*cos_md + 2*q2*sin_md;


	//4th row
	this->meas_jacobian_(3, 7) = 1;
	//5th row
	this->meas_jacobian_(4, 8) = 1;
	//6th row
	this->meas_jacobian_(5, 9) = 1;

	//7th row
	this->meas_jacobian_(6, 10) = 1;
	//8th row
	this->meas_jacobian_(7, 11) = 1;
	//9th row
	this->meas_jacobian_(8, 12) = 1;
}

template <typename T>
void Ekf16DofQuat<T>::ComputeMeasNoiseJacobian(MatrixInv<T> meas_sensor_val){
	// Do Nothing
}

template <typename T>
void Ekf16DofQuat<T>::ComputeMeasFromState(MatrixInv<T> time_propagated_state){

	// Easy to use names
	T q0 	= this->time_propagated_state_(0);
	T q1 	= this->time_propagated_state_(1);
	T q2 	= this->time_propagated_state_(2);
	T q3 	= this->time_propagated_state_(3);

	MatrixInv<T> c_ned2b = { { static_cast<T>( 1 - 2*( pow(q2, 2) + pow(q3,2) ) ), 2*( q1*q2 + q3*q0 ), 2*( q1*q3 - q2*q0 ) },
								   { 2*( q1*q2 - q3*q0 ), static_cast<T>( 1 - 2*( pow(q1, 2) + pow(q3, 2) ) ), 2*( q2*q3 + q1*q0 ) },
								   { 2*( q1*q3 + q2*q0 ), 2*( q2*q3 - q1*q0 ), static_cast<T>( 1 - 2*( pow(q1, 2) + pow(q2,2) ) ) } };




	T mag_declination = 13.01*DEG2RAD; // this can be an input later
	MatrixInv<T> c_mag2ned = { {cos(-mag_declination), sin(-mag_declination), 0},
              				   {-sin(-mag_declination), cos(-mag_declination), 0},
                       		   {0                         ,0                         , 0} };
    MatrixInv<T> unit_mag_vector = {{1}, {0}, {0}};
    MatrixInv<T> mag3D_unitVector_in_body = c_ned2b*c_mag2ned*unit_mag_vector;

    for(size_t idx = 0; idx < 3; idx++){
		this->meas_from_propogated_state_(idx) = mag3D_unitVector_in_body(idx);
		this->meas_from_propogated_state_(idx + 3) = this->time_propagated_state_(7 + idx);
		this->meas_from_propogated_state_(idx + 6) = this->time_propagated_state_(10 + idx);
	}
}

template <typename T>
void Ekf16DofQuat<T>::Run(MatrixInv<T> state_sensor_val, MatrixInv<T> meas_sensor_val, bool meas_indices[]){
	// Call the base class run method
	EkfBase<T>::Run(state_sensor_val, meas_sensor_val, meas_indices);
	//Normalize the quaternion
	T quat_mag = sqrt( pow(this->current_state_(0), 2) + pow(this->current_state_(1), 2) + pow(this->current_state_(1), 2) + pow(this->current_state_(3), 2) );
	for(size_t idx = 0; idx < 4; idx++)
		this->current_state_(idx) = this->current_state_(idx)/quat_mag;
}

template <typename T>
MatrixInv<T> Ekf16DofQuat<T>::GetEulerAngle(){
	MatrixInv<T> euler_ang(3, 1);
	// Easy to use names
	T q0 	= this->current_state_(0);
	T q1 	= this->current_state_(1);
	T q2 	= this->current_state_(2);
	T q3 	= this->current_state_(3);

	MatrixInv<T> c_ned2b 	= { { static_cast<T>( 1 - 2*( pow(q2, 2) + pow(q3,2) ) ), 2*( q1*q2 + q3*q0 ), 2*( q1*q3 - q2*q0 ) },
								{ 2*( q1*q2 - q3*q0 ), static_cast<T>( 1 - 2*( pow(q1, 2) + pow(q3, 2) ) ), 2*( q2*q3 + q1*q0 ) },
								{ 2*( q1*q3 + q2*q0 ), 2*( q2*q3 - q1*q0 ), static_cast<T>( 1 - 2*( pow(q1, 2) + pow(q2,2) ) ) } };

	euler_ang(0)  			= atan2( c_ned2b(1, 2), c_ned2b(2, 2) );
    euler_ang(1) 			= asin( -c_ned2b(0, 2) );
    euler_ang(2)   			= atan2( c_ned2b(0, 1), c_ned2b(0, 0) );

    if ( abs( c_ned2b(0, 2) ) > 1 - 1e-8 ){
        //Pitch=+/-90deg case.  Underdetermined, so assume roll is zero,
        //and solve for pitch and yaw as follows:
        euler_ang(0)   	= 0;
        euler_ang(1)  	= atan2( -c_ned2b(0, 2), c_ned2b(2, 2) );
        euler_ang(2)   	= atan2( -c_ned2b(1, 0), c_ned2b(1, 1) );
    }

    return euler_ang;

}

// Explicit template instantiation
template class Ekf16DofQuat<float>;
template class Ekf16DofQuat<double>;