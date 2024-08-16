#ifndef __ANGULAR_REPULSION_H__
#define __ANGULAR_REPULSION_H__

#include <map>
#include <string> // for std::string

#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h"
#include "hoomd/Compute.h"
#include "hoomd/HOOMDMPI.h"

#include <iostream> // for std::cout
#include <stdexcept> // for std::out_of_range
#include <algorithm> // for std::min and std::max
#include <vector>    // for std::vector


class AngularRepulsionManager {

public:
    AngularRepulsionManager() {}  // Default constructor

    //! Define functions that are used to calculate the angular repulsion potential 

    double getHeaviside(double r, double radsum) {
        // create the Heaviside step function (2d-r)
        // (this is for the angular repulsion potential)
        //return (2*radsum - r) < 0 ? 0.0 : 1.0;
        return (2 - r) < 0 ? 0.0 : 1.0;
        }

    double getDiracDelta(double r, double radsum) {
        // create the Heaviside step function (2d-r)
        // create the derivative of H(2d-r), which is
        // a Dirac delta function d(H(2d-r)/dr=-δ(r - 2d)
        // (this is for the angular repulsion force, effectively a second contact force)
        double epsilon = 1e-5;
        //return std::abs(r - 2*radsum) < epsilon ? -1.0 / epsilon : 0.0;
        return std::abs(r - 2) < epsilon ? -1.0 / epsilon : 0.0;
        }

    double getLambda(double x, double radsum) {
        // create the radial modulation function (Lambda function)
        // Λ(x) = x^-10 * [1-(0.5x)^10]^2 * heaviside(2-x)
        // (this is for the angular repulsion potential)
        if (x == 0) return 0;  // Handle x == 0 separately to avoid division by zero
        double r = x/radsum;
        //double h_lambda = getHeaviside(x, radsum);
        double h_lambda = getHeaviside(r, radsum);
        double term = std::pow(0.5 * r, 10);
        return std::pow(r, -10) * std::pow(1 - term, 2) * h_lambda;
        }
            
    double getLambdaPrime(double x, double radsum) {
        // and create it's derivative
        // Λ'(x) = -10x^-11 * [1-(0.5x)^10]^2 * heaviside(2d-x)
        //            + r^-10 * ( 2[1-(r/2)^10]*(-5(r/2)^9 ) * heaviside(2d-x)
        //               + r^-10 * [1-(r/2)^10]^2 * neg_dirac_delta
        // (this is for the angular repulsion force)
        if (x == 0) return 0;  // Handle x == 0 separately to avoid division by zero

        double r = x/radsum;
        //double h_lambdaprime = getHeaviside(x, radsum);
        //double dirac_lambdaprime = getDiracDelta(x, radsum);
        double h_lambdaprime = getHeaviside(r, radsum);
        double dirac_lambdaprime = getDiracDelta(r, radsum);

        double halfterm_ten = std::pow(0.5 * r, 10);
        double halfterm_nine = std::pow(0.5 * r, 9);

        double second_term = std::pow(1 - halfterm_ten, 2);
        double second_term_prime = -(5 * halfterm_nine);

        double part_one = -10 * std::pow(r, -11) * second_term * h_lambdaprime;
        double part_two = std::pow(r, -10) * (2 * second_term * second_term_prime) * h_lambdaprime;
        double part_three = std::pow(r, -10) * second_term * dirac_lambdaprime;

        return part_one + part_two + part_three;
        }


    bool getEnergyAndForce(double r_mag, double r_prime_mag, double dot_product, double radsum, double B, double theta_bar, double w, double theta, double angle_potential, double angle_force) {
        double cos_theta = std::cos(theta);
        //if (cos_theta != std::cos(theta_bar)) 
        if (theta != theta_bar) 
            {
            //std::cout << "Angular repulsion will be applied" << std::endl;

            // calculate energy
            double lambda_r = getLambda(r_mag, radsum);
            double lambda_r_prime = getLambda(r_prime_mag, radsum);
            double wsq_inv = std::pow(w, -2);
            //double heaviside_r = heaviside(r_mag);
            //double heaviside_r_prime = heaviside(r_prime_mag);

            //std::cout << "heaviside(r): " << heaviside_r << std::endl;
            //std::cout << "heaviside(r'): " << heaviside_r_prime << std::endl;
            //std::cout << "lambda(r): " << lambda_r << std::endl;
            //std::cout << "lambda(r'): " << lambda_r_prime << std::endl;
            //std::cout << "1/w^2': " << wsq_inv << std::endl;

            double angular_repulsion_exp = std::exp(-std::pow(cos_theta - std::cos(theta_bar), 2) * wsq_inv ); 
            angle_potential = B * lambda_r * lambda_r_prime * angular_repulsion_exp; 
            //std::cout << "Morse Energy: " << pair_eng << std::endl;
            //std::cout << "Anglular Repulsion Potential: " << angle_potential << std::endl;
            //pair_eng = pair_eng + angle_potential;
            //std::cout << "Total Energy: " << pair_eng << std::endl;
            //std::cout << "--" << std::endl;

            // calculate force
            double lambdaprime_r = getLambdaPrime(r_mag, radsum);

            //double neg_dirac_delta_r = neg_dirac_delta(r_mag);
            //double neg_dirac_delta_r_prime = neg_dirac_delta(r_prime_mag);
            //std::cout << "-dirac(r): " << neg_dirac_delta_r << std::endl;
            //std::cout << "-dirac(r'): " << neg_dirac_delta_r_prime << std::endl;

            double ddr_cos_theta = ( ((r_prime_mag * cos_theta) * (r_mag * r_prime_mag)) - (dot_product * r_prime_mag) ) / std::pow(r_mag * r_prime_mag , 2);
            double ddr_exp_contents = -2 * (cos_theta - std::cos(theta_bar)) * wsq_inv * ddr_cos_theta;
            double ddr_full_exp_term = angular_repulsion_exp * ddr_exp_contents;
            double first_force_term = lambdaprime_r * angular_repulsion_exp;
            double second_force_term = lambda_r * ddr_full_exp_term;
            //angle_forcei_divr = B * lambda_r_prime * (first_force_term + second_force_term) / r_mag;
            angle_force = B * lambda_r_prime * (first_force_term + second_force_term);
            //std::cout << "Morse Force: " << force_divr << std::endl;
            //std::cout << "Anglular Repulsion Force: " << angle_force << std::endl;
            //force_divr = force_divr + angle_force;
            //std::cout << "Total Force: " << force_divr << std::endl;
            //std::cout << "--------" << std::endl;
            return true;
            }
        else 
        {
            return false;
        }
    }

};

extern AngularRepulsionManager angular_repulsion;  // Declaration of a global instance of AngleManager

#endif // __ANGULAR_REPULSION_H__

