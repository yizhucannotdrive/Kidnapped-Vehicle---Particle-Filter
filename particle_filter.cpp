/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include "particle_filter.h"
#include <array>

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 100; // number of particles which could be adjusted.
    default_random_engine gen;
    weights.resize(num_particles);
    normal_distribution<double> dist_init_x(x, std[0]);
    normal_distribution<double> dist_init_y(y, std[1]);
    normal_distribution<double> dist_init_theta(theta, std[2]);
    struct Particle Particle_temp;
    for (int i=0; i<num_particles; i++){
        double sample_x = dist_init_x(gen);
        double sample_y = dist_init_y(gen);
        double sample_theta = dist_init_theta(gen);
        Particle_temp.id = i;
        Particle_temp.x = sample_x;
        Particle_temp.y = sample_y;
        Particle_temp.theta = sample_theta;
        Particle_temp.weight = 1.;
        particles.push_back(Particle_temp);
        weights[i]=1.0;
    }
    is_initialized = true;
    return;
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;
    normal_distribution<double> dist_pred_x(0, std_pos[0]);
    normal_distribution<double> dist_pred_y(0, std_pos[1]);
    normal_distribution<double> dist_pred_theta(0, std_pos[2]);
    for (int i=0; i<num_particles; i++){
        if (yaw_rate !=0){
        particles[i].x+= velocity/yaw_rate * (sin( particles[i].theta+ yaw_rate*delta_t)-sin(particles[i].theta))+dist_pred_x(gen);
        particles[i].y+= velocity/yaw_rate * (cos(particles[i].theta)-cos( particles[i].theta+ yaw_rate*delta_t))+dist_pred_y(gen);
        particles[i].theta += yaw_rate * delta_t;
        }
        else{
            particles[i].x+=velocity*cos(particles[i].theta)*delta_t;
            particles[i].y+=velocity*sin(particles[i].theta)*delta_t;
        }
        
    }
    //cout<<"check1"<< endl;
    return;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    int num_obs = observations.size();
    //cout<<"check1.3.1"<< num_obs<<predicted.size()<<endl;
    //cout<<"32th predicted landmark"<<predicted[31].id<<predicted[31].x<<predicted[31].y<<endl;
    //cout<<"39th predicted landmark"<<predicted[38].id<<predicted[38].x<<predicted[38].y<<endl;
       for (int i=0; i< num_obs; i++){
        double dist_temp = dist( predicted[0].x, predicted[0].y, observations[i].x, observations[i].y);
        observations[i].id = predicted[0].id;
        //cout<<"check1.3.1.1"<<"......"<< i<<endl;
                for (int j=0; j< predicted.size();j++){
            if (dist( predicted[j].x, predicted[j].y, observations[i].x, observations[i].y)< dist_temp){
                dist_temp = dist( predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
                observations[i].id = predicted[j].id;
                //cout<<"check1.3.1.2"<<"....."<< j<<endl;;
            }
        }
    }
    //cout<<"transformed_observation[i]_0_id_x_y"<<observations[0].id<< "....."<<observations[0].x<<" ....."<<observations[0].y<<endl;
    //cout<<"transformed_observation[i]_1_id_x_y"<<observations[1].id<< "....."<<observations[1].x<<" ....."<<observations[1].y<<endl;
    //cout<<"transformed_observation[i]_2_id_x_y"<<observations[2].id<< "....."<<observations[2].x<<" ....."<<observations[2].y<<endl;
    //cout<<"transformed_observation[i]_3_id_x_y"<<observations[3].id<< "....."<<observations[3].x<<" ....."<<observations[3].y<<endl;
    //cout<<"transformed_observation[i]_4_id_x_y"<<observations[4].id<< "....."<<observations[4].x<<" ....."<<observations[4].y<<endl;
    return;

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    double sum_weight = 0.0;
    double denom =1/(2*M_PI*std_landmark[0]*std_landmark[1]);
    double var1 = 2*std_landmark[0]*std_landmark[0];
    double var2 = 2*std_landmark[1]*std_landmark[1];
    std::vector <double> newweights;
    newweights.resize(num_particles);
    if (observations.size()==0){
        return;
    }
    for (int i=0; i<num_particles; i++ ){
        std::vector<LandmarkObs> predicted_particle_temp;
        for (int j =0; j< map_landmarks.landmark_list.size(); j++){
            if(dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, particles[i].x, particles[i].y)< sensor_range){
                LandmarkObs obs_temp;
                obs_temp.id = map_landmarks.landmark_list[j].id_i;
                obs_temp.x = map_landmarks.landmark_list[j].x_f;
                obs_temp.y = map_landmarks.landmark_list[j].y_f;
                predicted_particle_temp.push_back(obs_temp);
            }
        }
        
        std::vector<LandmarkObs> transformed_observations;
        //cout<< "check1.1"<<endl;
        for (int j=0; j< observations.size(); j++){
            double obj_x_m_temp;
            double obj_y_m_temp;
            LandmarkObs obs_temp_;
            obj_x_m_temp = particles[i].x+ cos(particles[i].theta)*  observations[j].x- sin(particles[i].theta)*  observations[j].y;
            obj_y_m_temp = particles[i].y+ sin(particles[i].theta)*  observations[j].x+ cos(particles[i].theta)*  observations[j].y;
            obs_temp_.x = obj_x_m_temp;
            obs_temp_.y = obj_y_m_temp;
            transformed_observations.push_back(obs_temp_);
        }
        //cout<< "check1.2"<<endl;
        double weight_i_new = 1.0;
        if(predicted_particle_temp.size()!=0){
            dataAssociation( predicted_particle_temp, transformed_observations);
            //cout<< transformed_observations.size()<< endl;
            //cout<< "check1.3"<<endl;
            for (int j=0; j< transformed_observations.size(); j++){
                int id_of_jth_obs = transformed_observations[j].id;
                //cout<< id_of_jth_obs<< endl;
                double x_diff = transformed_observations[j].x - map_landmarks.landmark_list[id_of_jth_obs-1].x_f;
                double y_diff = transformed_observations[j].y - map_landmarks.landmark_list[id_of_jth_obs-1].y_f;
                //cout<< x_diff<< endl;
                double likelihood = denom *exp(-(x_diff*x_diff/var1+y_diff*y_diff/var2));
                weight_i_new *= likelihood;
            }

        }
        else{
           
            weight_i_new=0.0;
        }
        //cout<< " current particle position"<< particles[i].x<< particles[i].y<<endl;
        //cout<< " num_ observed landmarks"<< predicted_particle_temp.size()<<endl;
        //cout<<"transformed_observation[i]_0_id_x_y"<<transformed_observations[0].id<< "....."<<transformed_observations[0].x<<" ....."<<transformed_observations[0].y<<endl;
        newweights[i] = weight_i_new;
        //if (weight_i_new==0){
         //   cout<< "which_particle"<< i<< endl;
        //    cout<<"predicted_particle_temp.size"<<predicted_particle_temp.size()<<endl;
        //    cout<<"transformed_observations.size"<<transformed_observations.size()<<endl;
        //}
        sum_weight+= weight_i_new;
    }
    
    if(sum_weight<0.0000000001){
        return;
    }
    else
    //cout<< "check1.4"<<endl;
    {
        for (int i=0; i<num_particles; i++ ){
        weights[i] = newweights[i]/sum_weight;
        particles[i].weight = weights[i];
        
    }
    return;
    }
}
    //cout<<"check2"<< endl;


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::vector<Particle> particles_resample;
    particles_resample.resize(num_particles);
    default_random_engine gen;
    discrete_distribution< int > resample_id(weights.begin(),weights.end());
    struct Particle Particle_temp;
    for (int i=0; i<num_particles; i++){
        int index = resample_id(gen);
        particles_resample.at(i) = particles.at(index);
    }
    particles = particles_resample;    //cout<<"check3"<< endl;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
