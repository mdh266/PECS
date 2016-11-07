#include "../include/Generation.hpp"

using namespace dealii;

	
template <int dim>
void 
Generation<dim>::
set_dark_params() 
{	
	scaled_photon_flux = 0.0;
	scaled_absorption_coeff	= 0.0; 
	scaled_incident_location = 0.0;
}

	
template <int dim>
void
Generation<dim>::
set_illuminated_params(const ParameterSpace::Parameters & params)
{	
	scaled_photon_flux = params.scaled_photon_flux;
	scaled_absorption_coeff	= params.scaled_absorption_coeff;
	scaled_incident_location = params.scaled_domain_height;
//	scaled_incident_location = params.scaled_radius_one;

}

template <int dim>
double
Generation<dim>::
value(const Point<dim> &p, const unsigned int ) const
{
	// G_0 * aplpha * exp(\alpha ( 1 - y))
//	const double value = 

	return scaled_absorption_coeff *
				 scaled_photon_flux  *
				exp(scaled_absorption_coeff * 
				 (p[1] - scaled_incident_location));

//	std::cout << "g = " << value << std::endl;
//	return value;
}

