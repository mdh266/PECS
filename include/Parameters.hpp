#ifndef _PARAMETERS_H__
#define _PARAMETERS_H__

#include<deal.II/base/parameter_handler.h>
#include<iostream>

namespace PhysicalConstants
{
	const double thermal_voltage  	 = 0.02585; // [V]
	const double electron_charge  	 = 1.62e-19;  // [C]
	const double vacuum_permittivity = 8.85e-14; // [A s V^{-1} cm^{-1} 
}

namespace ParameterSpace
{
	using namespace dealii;

	/// \brief Struct which holds the parameters used for simulations.
	
	
	/** These are the parameters which will be used by
	*		DriftDiffusionPoisson class to do simulations. They will be read
	*		read in using the ParameterReader class to read them from 
	*		input_file.prm and then use the ParameterHandler to 
	* 	set the variables and scale them if necessary.
	*/							 
	struct Parameters
	{
			// computational
			unsigned int n_global_refine;				
			unsigned int n_local_refine;				
			unsigned int time_stamps;
			double 	 	 	 h_max;
			double			 h_min;
			double 			 t_end;
			double 			 t_end_2;
			double 			 delta_t;	
			double 			 penalty;

			double scaled_electron_mobility;
			double scaled_electron_recombo_t;
			double scaled_electron_recombo_v;
			double scaled_k_et;			

			double scaled_hole_mobility;
			double scaled_hole_recombo_t;
			double scaled_hole_recombo_v;
			double scaled_k_ht;

			double scaled_intrinsic_density;
			double semiconductor_permittivity;
					
			double scaled_reductant_mobility;
			double scaled_oxidant_mobility;
			double electrolyte_permittivity;


			double scaled_absorption_coeff;
			double scaled_photon_flux;
	
			double scaled_debeye_length;
			double scaled_boundary_layer;			


			double characteristic_length;
			double characteristic_time;
			double characteristic_denisty;

			double scaled_domain_length;
			double scaled_domain_height;
			double scaled_radius_one;
			double scaled_radius_two;	

			bool illum_or_dark;
			bool insulated;
			bool restart_status;
			bool schottky_status;
	
			double scaled_applied_bias;
			double scaled_built_in_bias;
			double scaled_schottky_bias;

			double rescale_current;
			double rescaled_k_et;
			double rescaled_k_ht;

			/** This functions makes almost everythingn 1.0 that is relevant for testing.*/
			void set_params_for_testing(const unsigned int & n_refine)
			{
				n_global_refine 	 = n_refine;
				t_end 			 = 1.0;
				scaled_electron_mobility = 1.0;
				scaled_hole_mobility	 = 1.0;
				scaled_oxidant_mobility  = 1.0;
				scaled_reductant_mobility = 1.0;
				scaled_absorption_coeff  = 0.0;
	
				scaled_domain_height  = 1.0;	
				scaled_domain_length  = 1.0;
				scaled_radius_one     = 0.5;
				scaled_radius_two     = 0.5;

				scaled_debeye_length   = 1.0;
				characteristic_length  = 1.0;
				characteristic_denisty = 1.0;
				characteristic_time    = 1.0;
			}

			/** This function opens <code>input_file.prm<code/>, 
			*		reads in the specified 
			* 	parameter values and scales appopriate values for singular
			* 	perturbation scaling.
			* 
			*   NOTE: This must be called from the constructor: 
			*		SOLARCELL::SolarCellProblem()
			*/
			void parse_and_scale_parameters(ParameterHandler & prm)
			{

				// read in the parameters
				prm.enter_subsection("computational");
				this->n_global_refine = prm.get_integer("global refinements");
				this->n_local_refine  = prm.get_integer("local refinements");
				this->delta_t  = prm.get_double("time step size");
				this->t_end = prm.get_double("end time");
				this->t_end_2 = prm.get_double("end time 2");
				this->time_stamps = prm.get_integer("time stamps");
				this->restart_status = prm.get_bool("restart status");
				prm.leave_subsection();

				prm.enter_subsection("mesh");
				this->scaled_domain_height = prm.get_double("mesh height");
				this->scaled_domain_length = prm.get_double("mesh length");
				this->scaled_radius_one = prm.get_double("radius one");
				this->scaled_radius_two = prm.get_double("radius two");
				this->scaled_boundary_layer = prm.get_double("boundary layer");
				prm.leave_subsection();

				prm.enter_subsection("physical");	
				this->illum_or_dark = prm.get_bool("illumination status");
				this->insulated = prm.get_bool("insulated");
				this->schottky_status = prm.get_bool("schottky status");
				this->scaled_applied_bias = prm.get_double("applied bias");
				this->scaled_built_in_bias = prm.get_double("built in bias");
				this->scaled_schottky_bias = prm.get_double("schottky bias");
				this->characteristic_length = prm.get_double("characteristic length");
				this->characteristic_denisty = prm.get_double("characteristic density");
				this->characteristic_time = prm.get_double("characteristic time");
				this->scaled_intrinsic_density = prm.get_double("intrinsic density");
				this->scaled_photon_flux = prm.get_double("photon flux");
				this->scaled_absorption_coeff = prm.get_double("absorption coefficient"); 
				this->semiconductor_permittivity = prm.get_double("semiconductor permittivity"); 
				this->electrolyte_permittivity = prm.get_double("electrolyte permittivity"); 
				prm.leave_subsection();

					
				prm.enter_subsection("electrons"); 
				this->scaled_electron_mobility = prm.get_double("mobility");
				this->scaled_k_et = prm.get_double("transfer rate");
				this->scaled_electron_recombo_t = prm.get_double("recombination time");
				this->scaled_electron_recombo_v = prm.get_double("recombination velocity");
				prm.leave_subsection();

				prm.enter_subsection("holes"); 
				this->scaled_hole_mobility = prm.get_double("mobility");
				this->scaled_k_ht = prm.get_double("transfer rate");
				this->scaled_hole_recombo_t = prm.get_double("recombination time");
				this->scaled_hole_recombo_v = prm.get_double("recombination velocity");
				prm.leave_subsection();

				prm.enter_subsection("reductants"); 
				this->scaled_reductant_mobility = prm.get_double("mobility");
				prm.leave_subsection();


				prm.enter_subsection("oxidants"); 
				this->scaled_oxidant_mobility = prm.get_double("mobility");
				prm.leave_subsection();


				// scale the parameters
				this->scaled_intrinsic_density /= this->characteristic_denisty;
				this->scaled_electron_recombo_t /= this->characteristic_time;
				this->scaled_hole_recombo_t /= this->characteristic_time;

				this->scaled_electron_recombo_v *= (this->characteristic_time /
								this->characteristic_length);

				this->scaled_hole_recombo_v *= (this->characteristic_time /
								this->characteristic_length);

				this->scaled_photon_flux *= (this->characteristic_time /
							     this->characteristic_denisty);

				this->scaled_absorption_coeff *= this->characteristic_length;

				// DONT INCLUDE THE MATERIAL DIELECTRIC IN THE DEBEYE LENGTH
				this->scaled_debeye_length = 
						(PhysicalConstants::thermal_voltage * 
						PhysicalConstants::vacuum_permittivity) / //* 
//						this->material_permittivity) /
						(PhysicalConstants::electron_charge * 
						this->characteristic_denisty * 
						this->characteristic_length *
						this->characteristic_length);

				double mobility_scale = (this->characteristic_time *
						PhysicalConstants::thermal_voltage) /
						(this->characteristic_length *
						this->characteristic_length);

				this->scaled_electron_mobility *= mobility_scale;
				this->scaled_hole_mobility *= mobility_scale;
				this->scaled_reductant_mobility *= mobility_scale;
				this->scaled_oxidant_mobility *= mobility_scale;

				this->rescaled_k_et  = PhysicalConstants::electron_charge *
							this->scaled_k_et * 
							this->characteristic_denisty *
							this->characteristic_denisty;

				this->rescaled_k_ht  = PhysicalConstants::electron_charge *
							this->scaled_k_ht *
							this->characteristic_denisty *
							this->characteristic_denisty;
	
				this->scaled_k_et *= (this->characteristic_time *
							this->characteristic_denisty /
							this->characteristic_length);
														
				this->scaled_k_ht *= (this->characteristic_time *
							this->characteristic_denisty /
							this->characteristic_length);

				this->scaled_applied_bias /= PhysicalConstants::thermal_voltage;
				this->scaled_built_in_bias /= PhysicalConstants::thermal_voltage;
				this->scaled_schottky_bias /= PhysicalConstants::thermal_voltage;

				this->rescale_current = (PhysicalConstants::electron_charge 
							* this->characteristic_denisty 
							* this->characteristic_length)
							/ this->characteristic_time;


/*
				std::cout << "debeye length = " 
					<< this->scaled_debeye_length 
					<< std::endl;
				std::cout << "semiconductor perm = "
					<< this->semiconductor_permittivity
					<< std::endl;
				std::cout << "electrolyte perm = "
					<< this->electrolyte_permittivity
					<< std::endl;
				std::cout << "scaled electron mobility = " 
					<< this->scaled_electron_mobility 
					<< std::endl;
				std::cout << "scaled hole mobility = " 
					<< this->scaled_hole_mobility 
					<< std::endl;
				std::cout << "scaled reductant mobility = " 
					<< this->scaled_reductant_mobility 
					<< std::endl;
				std::cout << "scaled oxidant mobility = " 
					<< this->scaled_oxidant_mobility 
					<< std::endl;
				std::cout << "k_et = " 
					<< this->scaled_k_et
					<< std::endl;
				std::cout << "k_ht = " 
					<< this->scaled_k_ht
					<< std::endl;
				std::cout << "built in bias = " 
					<< this->scaled_built_in_bias 
					<< std::endl;
				std::cout << "rescale current = "
					<< this->rescale_current
					<< std::endl;
				std::cout << "photon flux = " 
					<< this->scaled_photon_flux
					<< std::endl;
				std::cout << "absorption coeff = "
					<< this->scaled_absorption_coeff 
					<< std::endl;
*/

		} // parse_and_scale_parameters(prm)
		
	};

	/// @author Michael Harmon
}


#endif
