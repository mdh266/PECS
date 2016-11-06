#ifndef _POSTPROCESSOR_H__
#define _POSTPROCESSOR_H__

//#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/exceptions.h>
#include "Parameters.hpp"
#include "CarrierPair.hpp"
#include <string>

using namespace dealii;

/** \brief This will add back in units and calculate the current.*/
/** \note This is based off of step-32 in deal.ii examples.*/
template <int dim>
class PostProcessor : public DataPostprocessor<dim>
{
	public:

		/** \brief Constructor for Poisson post processor.*/
		/** The constructor will assign the scaling values to the object.*/
		PostProcessor(const ParameterSpace::Parameters 			& sim_parms,
					  const bool								& print_carrier,
					  const std::string							& name);

		/** \brief post processing done in this function.*/
		/** \note: This is called internally from DatOut::build_patches*/
		virtual
		void 
		compute_derived_quantities_vector(
				const std::vector<Vector<double>>				&uh,
				const std::vector<std::vector<Tensor<1,dim>>>	&duh,
				const std::vector<std::vector<Tensor<2,dim>>>	&dduh,
				const std::vector<Point<dim>>					&normals,
				const std::vector<Point<dim>>					&evaluated_points,
				std::vector<Vector<double>>						&computed_quantities) const;


		/** \brief Returns a vector containing the names of the data componetnt.*/
		virtual std::vector<std::string> get_names() const;
	
		/** \brief returns a vector which tells whether data comp. is scalar or vector.*/
		virtual
		std::vector<DataComponentInterpretation::DataComponentInterpretation>
		get_data_component_interpretation() const;

		/** \brief Returns the values which must be update every print.*/
		virtual UpdateFlags get_needed_update_flags() const;

	private:
		double scale_density;
		double scale_current;
		double scale_potential;
		double scale_elec_field;
		double material_permittivity;
		std::string density_name;
		std::string current_name;
		bool printing_carrier;

};

#endif 
