#include "../include/BiasValues.hpp"

template<int dim>
void
Built_In_Bias<dim>::
set_value(const double & bias_value)
{
	built_in_bias = bias_value;
}

template <int dim>
double 
Built_In_Bias<dim>::
value(const dealii::Point<dim> &p, 
	  const unsigned int ) const
{
	// potential applied at y = 0, a positive potential is a foward bias,
	// a negative potential is a reverse bias
	if((p[0] == 0.0)) //|| (p[1] == 0.0))
	{
//		std::cout << "Built in" << std::endl;
		return built_in_bias; // volts

	}
	else 
		return 0;
//	double return_value = 0.0;
//	for(unsigned int i = 0; i < dim; i++)
//		return_value += 20.0*p[i];
//	return return_value;
}

template<int dim>
void
Schottky_Bias<dim>::
set_location(const double & bias_location)
{
	Schottky_location = bias_location;
}

template<int dim>
void
Schottky_Bias<dim>::
set_value(const double & bias_value)
{
	Schottky_bias = bias_value;
}

template <int dim>
double 
Schottky_Bias<dim>::
value(const dealii::Point<dim> &p, 
	  const unsigned int ) const
{
	// potential applied at y = 0, a positive potential is a foward bias,
	// a negative potential is a reverse bias
	if(p[1] == Schottky_location)
	{
//		std::cout << "on schottky" << std::endl;
		return Schottky_bias; // volts

	}
	else 
		return 0;
}

template<int dim>
void
Applied_Bias<dim>::
set_value(const double & bias_value)
{
	applied_bias = bias_value;
}

template <int dim>
double 
Applied_Bias<dim>::
value(const dealii::Point<dim> &p, 
	  const unsigned int ) const
{
	// potential applied at y = 0, a positive potential is a foward bias,
	// a negative potential is a reverse bias
	if(p[0] == 0.0) //|| (p[1] == 0))
	{
//		std::cout << "Applied" << std::endl;
		return applied_bias; // volts
	}
	else 
		return 0;
//	double return_value = 0.0;
//	for(unsigned int i = 0; i < dim; i++)
//		return_value += 20.0*p[i];
//	return return_value;
}


// TODO: THIS DOESNT WORK SO WELL FOR ANY NONZERO SINCE
// NEEDS TO BE ON SCALED_DOMAIN_LENGTH
template <int dim>
double 
Bulk_Bias<dim>::
value(const dealii::Point<dim> &p, 
	  const unsigned int ) const
{
	// potential applied at y = 0, a positive potential is a foward bias,
	// a negative potential is a reverse bias
	if(p[0] == 1.0)
		return 0.0;
	else 
		return 0;
//	double return_value = 0.0;
//	for(unsigned int i = 0; i < dim; i++)
//		return_value += 20.0*p[i];
//	return return_value;
}

