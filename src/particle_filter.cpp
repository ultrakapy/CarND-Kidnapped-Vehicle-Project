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
    num_particles = 42*42; // map size squared

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
        double new_x = x_0 + ((velocity / yaw_rate) * (sin(theta_0 + (yaw_rate * delta_t)) - sin(theta_0)));
        double new_y = y_0 + ((velocity / yaw_rate) * (cos(theta_0) - cos(theta_0 + (yaw_rate * delta_t))));
        double new_theta = theta_0 + (yaw_rate * delta_t);

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
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
