#include "../include/PostProcessor.hpp"

using namespace dealii;


template<int dim>
PostProcessor<dim>::
PostProcessor(const ParameterSpace::Parameters 		& sim_params,
	     const bool					& print_carrier,
	     const std::string				& name)
{
	scale_potential  = 0.02585;
	scale_elec_field = 0.2585/sim_params.characteristic_length;
	scale_density	 = sim_params.characteristic_denisty;
	scale_current	 = 1.6e-19* scale_density * sim_params.characteristic_length /
			   sim_params.characteristic_time;

	printing_carrier = print_carrier;	
	if(print_carrier)
	{
		density_name  = name.c_str();
		density_name += " Density";
		current_name  = name.c_str();
		current_name += " Current";
	}
}

template<int dim>
std::vector<std::string> 
PostProcessor<dim>::
get_names() const
{
	std::vector<std::string> solution_names;

	if(printing_carrier)
	{
		for(unsigned int d=0; d<dim; d++)
			solution_names.push_back(current_name.c_str());

		solution_names.push_back(density_name.c_str());	
	}
	else
	{
		for(unsigned int d=0; d<dim; d++)
			solution_names.push_back("Field");

		solution_names.push_back("Potential");
	} 
	return solution_names;
}

template<int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
PostProcessor<dim>::
get_data_component_interpretation() const
{

	std::vector<DataComponentInterpretation::DataComponentInterpretation>
	interpretation;

	for(unsigned int d=0; d<dim; d++)
	{
		interpretation.push_back(DataComponentInterpretation::
					component_is_part_of_vector);
	}
	interpretation.push_back(DataComponentInterpretation::component_is_scalar);

	return interpretation;
}

template<int dim>
UpdateFlags
PostProcessor<dim>::
get_needed_update_flags() const
{
	return update_values;// | update_gradients;
}


template<int dim>
void 
PostProcessor<dim>::
compute_derived_quantities_vector(
	const std::vector<Vector<double>>			&uh,
	const std::vector<std::vector<Tensor<1,dim>>>		&/*duh*/,
	const std::vector<std::vector<Tensor<2,dim>>>		&/* dduh*/,
	const std::vector<Point<dim>>				&/*normals*/,
	const std::vector<Point<dim>>				&/*evaluated_points*/,
	std::vector<Vector<double>>				&computed_quantities) const
{

	
	const unsigned int n_quadrature_points = uh.size();
	Assert(computed_quantities.size() == n_quadrature_points, ExcInternalError());
	Assert(uh[0].size() == (dim+1), ExcInternalError() );
	if(printing_carrier)
	{
		for(unsigned int q=0; q<n_quadrature_points; q++)
		{
			// copy over current with scaling.
			for(unsigned int d=0; d<dim; d++)
				computed_quantities[q](d) = (scale_current * uh[q](d));
	
			// copy over densities without scaling
			computed_quantities[q](dim) = uh[q](dim);
		} // for q
	}
	else
	{
		for(unsigned int q=0; q<n_quadrature_points; q++)
		{
			// copy over electric field..
			// TODO: rescale with permitivities?
			for(unsigned int d=0; d<dim; d++)
				computed_quantities[q](d) = (scale_elec_field * uh[q](d));
	
			// copy over potential
			const double potential = scale_potential * uh[q](dim);
			computed_quantities[q](dim) = potential;
		} // for q
	} // else

}

