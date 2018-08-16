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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    default_random_engine gen;
    num_particles = 100; //42*42; // map size squared

    // create normal distributions for x, y and theta.
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    // sample values for each particle from normal distributions
    for (int i = 0; i < num_particles; i++) {
        double sample_x = dist_x(gen);
        double sample_y = dist_y(gen);
        double sample_theta = dist_theta(gen);
        double init_weight = 1.0;
        //cout << "P (" << i << ") = " << sample_x << ", " << sample_y << ", " << sample_theta << "\n";
        Particle p = {i, sample_x, sample_y, sample_theta, init_weight};
        particles.push_back(p);

        weights.push_back(init_weight);
    }

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;

    for (int i = 0; i < num_particles; i++) {
        double x_0 = particles[i].x;
        double y_0 = particles[i].y;
        double theta_0 = particles[i].theta;
        double new_x, new_y, new_theta;

        if (yaw_rate == 0) {
            new_x = x_0 + velocity*delta_t*cos(theta_0);
            new_y = y_0 + velocity*delta_t*sin(theta_0);
            new_theta = theta_0;
        } else {
            new_x = x_0 + ((velocity / yaw_rate) * (sin(theta_0 + (yaw_rate * delta_t)) - sin(theta_0)));
            new_y = y_0 + ((velocity / yaw_rate) * (cos(theta_0) - cos(theta_0 + (yaw_rate * delta_t))));
            new_theta = theta_0 + (yaw_rate * delta_t);
        }

        normal_distribution<double> dist_new_x(new_x, std_pos[0]);
        particles[i].x = dist_new_x(gen);

        normal_distribution<double> dist_new_y(new_y, std_pos[1]);
        particles[i].y = dist_new_y(gen);

        normal_distribution<double> dist_new_theta(new_theta, std_pos[2]);
        particles[i].theta = dist_new_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    for (unsigned int i = 0; i < observations.size(); i++) {
        double min_dist = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y);
        int closest_predicted_landmark_id = predicted[0].id;

        for (unsigned int k = 1; k < predicted.size(); k++) {
            double curr_dist = dist(observations[i].x, observations[i].y, predicted[k].x, predicted[k].y);
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                closest_predicted_landmark_id = predicted[k].id;
            }
        }
        // set landmark ID of closest predicted measurement
        observations[i].id = closest_predicted_landmark_id;
    }
}

Map::single_landmark_s ParticleFilter::get_landmark_by_id(int id, const Map &map_landmarks) {
    for (unsigned int i = 0; i < map_landmarks.landmark_list.size(); i++) {
        if (map_landmarks.landmark_list[i].id_i == id) {
            return map_landmarks.landmark_list[i];
        }
    }
    //cout << "HERE ---------------------------\n";
    return map_landmarks.landmark_list[0];
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	for (int i = 0; i < num_particles; i++) {
        std::vector<LandmarkObs> predicted_meas;
        std::vector<LandmarkObs> transformed_observations;
        double rotation_angle = particles[i].theta;

        // find predicted measurements to landmarks for each particle
        for (unsigned int k = 0; k < map_landmarks.landmark_list.size(); k++) {
            double dist_to_landmark = dist(particles[i].x, particles[i].y,
                                           map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
            //cout << "DIST: ------ " << sensor_range << "\n";
            if (dist_to_landmark <= sensor_range) {
                LandmarkObs pred_landmark = {map_landmarks.landmark_list[k].id_i,
                                             map_landmarks.landmark_list[k].x_f,
                                             map_landmarks.landmark_list[k].y_f};
                predicted_meas.push_back(pred_landmark);
            }
        }

        // convert observations into map's coordinate system
        for (unsigned int o = 0; o < observations.size(); o++) {
            double x_map = particles[i].x + (cos(rotation_angle)*observations[o].x) - (sin(rotation_angle)*observations[o].y);
            double y_map = particles[i].y + (sin(rotation_angle)*observations[o].x) + (cos(rotation_angle)*observations[o].y);
            LandmarkObs lob = {observations[i].id, x_map, y_map};
            //cout << "LOB: " << lob.id << ", " << lob.x << ", " << lob.y << "\n";
            transformed_observations.push_back(lob);
        }

        dataAssociation(predicted_meas, transformed_observations);

        // set landmark associations to particle
        for (unsigned int n = 0; n < transformed_observations.size(); n++) {
            particles[i].associations.push_back(transformed_observations[n].id);
            particles[i].sense_x.push_back(transformed_observations[n].x);
            particles[i].sense_y.push_back(transformed_observations[n].y);
        }

        // calculate each measurement's Multivariate-Gaussian probability density & compute the final particle weight
        double weight = 1.0;
        for (unsigned int n = 0; n < particles[i].associations.size(); n++) {
            double mgpd;
            double sig_x = std_landmark[0];
            double sig_y = std_landmark[1];
            double x_obs = particles[i].sense_x[n];
            double y_obs = particles[i].sense_y[n];
            Map::single_landmark_s landmark = get_landmark_by_id(particles[i].associations[n], map_landmarks);
            double mu_x = landmark.x_f; //map_landmarks.landmark_list[transformed_observations[n].id - 1].x_f;
            double mu_y = landmark.y_f; //map_landmarks.landmark_list[transformed_observations[n].id - 1].y_f;
            double gauss_norm = (1.0 / (2.0 * M_PI * sig_x * sig_y));
            double exponent = (pow((x_obs - mu_x), 2))/(2 * pow(sig_x, 2)) + (pow((y_obs - mu_y), 2))/(2 * pow(sig_y, 2));
            mgpd = gauss_norm * exp(-exponent);
            weight *= mgpd;
        }

        particles[i].weight = weight;
        weights[i] = weight;
        //cout << "HERE " << getAssociations(particles[i]) << "\n";
	}
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    //std::random_device rd;
    //std::mt19937 gen(rd());
    default_random_engine gen;
    std::discrete_distribution<> d(weights.begin(), weights.end());
    std::vector<Particle> r_particles;

    for (int i = 0; i < num_particles; i++) {
        r_particles.push_back(particles[d(gen)]);
    }

    particles = r_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
