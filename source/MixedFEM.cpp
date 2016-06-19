#include "../include/MixedFEM.hpp"
#include "PostProcessor.cpp"
#include "Poisson.cpp"


namespace MixedPoisson
{
	using namespace std;

	//  mixed FEM constructor
	template<int dim>
	MixedFEM<dim>::
	MixedFEM()	:
	test_Poisson_solution(),
	test_Poisson_bc(),
	test_Poisson_rhs(),
	test_coupling_Poisson_rhs()	 //: instantiate function objects
	{
	}

	// mixed FEM destructor
	template<int dim>
	MixedFEM<dim>::
	~MixedFEM()
	{
		// nothing
	}

	template<int dim>
	void 
	MixedFEM<dim>::	
	assemble_local_Poisson_matrix(	
				const typename DoFHandler<dim>::active_cell_iterator & cell,
				Assembly::AssemblyScratch<dim>												 & scratch,
				Assembly::Poisson::CopyData<dim>											 & data,
				const double 																					 & semi_permittivity,
				const double 																					 & elec_permittivity,
				const double 																					 & scaled_debye_length)
	{
		const unsigned int	dofs_per_cell 	=	
												scratch.Poisson_fe_values.dofs_per_cell;
		const unsigned int 	n_q_points 			= 
												scratch.Poisson_fe_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		// Get the actual values for vector field and potential from FEValues
		// Use Extractors instead of having to deal with shapefunctions directly
		const FEValuesExtractors::Vector VectorField(0); // Vector as in Vector field
		const FEValuesExtractors::Scalar Potential(dim);

		scratch.Poisson_fe_values.reinit(cell);
		data.local_matrix=0;
		
		if(	(cell->material_id() == semiconductor_id)
				||
				(cell->material_id() == semi_boundary_layer_id))
		{
			// loop over all the quadrature points in this cell
			for(unsigned int q=0; q<n_q_points; q++)
			{
				//  loop over test functions dofs for this cell
				for(unsigned int i=0; i<dofs_per_cell; i++)
				{
					// i-th VectorField basis functions at point q
					const Tensor<1, dim>	psi_i_field = 
												scratch.Poisson_fe_values[VectorField].value(i,q);
					
					// div. of the i-th VectorField basis functions at the point q
					const double	div_psi_i_field =	
												scratch.Poisson_fe_values[VectorField].divergence(i,q);
		
					// i-th potential basis functions at the point q
					const double psi_i_potential = 
												scratch.Poisson_fe_values[Potential].value(i,q);
		
					// loop over all the trial functions dofs for this cell
					for(unsigned int j=0; j<dofs_per_cell; j++)	
					{
						// j-th VectorField basis functions at point q
						const Tensor<1, dim>	psi_j_field = 
												scratch.Poisson_fe_values[VectorField].value(j,q);
					
						// div. of the j-th VectorField basis functions at the point q
						const double	div_psi_j_field =	
												scratch.Poisson_fe_values[VectorField].divergence(j,q);
		
	
						// i-th potential basis functions at the point q
						const double psi_j_potential = 
												scratch.Poisson_fe_values[Potential].value(j,q);

						// build whole matrix at once, not blocks individually.
						// \int (P * k^{-1} D - div P * phi - v * D ) dx
						data.local_matrix(i,j) += (psi_i_field * 
																(1.0/semi_permittivity) * 
																psi_j_field 
																- div_psi_i_field * psi_j_potential
																- psi_i_potential *
																	scaled_debye_length * div_psi_j_field	
																) * scratch.Poisson_fe_values.JxW(q);
 						
					} // for j
				} // for i
			} // for q
		} // end if cell == semiconductor
		else if(	(cell->material_id() == electrolyte_id)
							||
							(cell->material_id() == elec_boundary_layer_id))
		{
			// loop over all the quadrature points in this cell
			for(unsigned int q=0; q<n_q_points; q++)
			{
				//  loop over test functions dofs for this cell
				for(unsigned int i=0; i<dofs_per_cell; i++)
				{
					// i-th VectorField basis functions at point q
					const Tensor<1, dim>	psi_i_field = 
												scratch.Poisson_fe_values[VectorField].value(i,q);
					
					// div. of the i-th VectorField basis functions at the point q
					const double	div_psi_i_field =	
												scratch.Poisson_fe_values[VectorField].divergence(i,q);
		
					// i-th potential basis functions at the point q
					const double psi_i_potential = 
												scratch.Poisson_fe_values[Potential].value(i,q);
		
					// loop over all the trial functions dofs for this cell
					for(unsigned int j=0; j<dofs_per_cell; j++)	
					{
						// j-th VectorField basis functions at point q
						const Tensor<1, dim>	psi_j_field = 
												scratch.Poisson_fe_values[VectorField].value(j,q);
					
						// div. of the j-th VectorField basis functions at the point q
						const double	div_psi_j_field =	
												scratch.Poisson_fe_values[VectorField].divergence(j,q);
		
	
						// i-th potential basis functions at the point q
						const double psi_j_potential = 
												scratch.Poisson_fe_values[Potential].value(j,q);

						// build whole matrix at once, not blocks individually.
						// \int (P * k^{-1} D - div P * phi - v * D ) dx
						data.local_matrix(i,j) += (psi_i_field * 
																(1.0/elec_permittivity) * 
																psi_j_field 
																- div_psi_i_field * psi_j_potential
																- psi_i_potential *
																	scaled_debye_length * div_psi_j_field	
																) * scratch.Poisson_fe_values.JxW(q);
 						
					} // for j
				} // for i
			} // for q			
		} // end material_id() == electrolyte
		else
		{
			// meed soemthing better
			std::cerr << "CELL TYPE NOT SEMICONDUCTOR OR ELECTROLYTE\n";
			Assert(false, ExcNotImplemented() );
		}
		
	} // assemble_local_Poisson_matrix

	template<int dim>
	void
	MixedFEM<dim>::
	assemble_local_test_rhs(
				const typename DoFHandler<dim>::active_cell_iterator 	& cell,
				Assembly::AssemblyScratch<dim>											  & scratch,
				Assembly::Poisson::CopyData<dim>											& data)
	{
		const unsigned int	dofs_per_cell 	=	
												scratch.Poisson_fe_values.dofs_per_cell;
		const unsigned int 	n_q_points 			= 
												scratch.Poisson_fe_values.n_quadrature_points;

		const unsigned int	n_face_q_points =	
												scratch.Poisson_fe_face_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		// Get the actual values for vector field and potential from FEValues
		// Use Extractors instead of having to deal with shapefunctions directly
		const FEValuesExtractors::Vector VectorField(0); // Vector as in Vector field
		const FEValuesExtractors::Scalar Potential(dim);


		// reinitialize the fe_Values on this cell
		scratch.Poisson_fe_values.reinit(cell);
	
		// reset the local_rhs vector to be zero
		data.local_rhs=0;
		// get the test rhs for poisson
		test_Poisson_rhs.value_list(
														scratch.Poisson_fe_values.get_quadrature_points(),
														scratch.Poisson_rhs_values);


		// Loop over all the quadrature points in this cell
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// loop over the test function dofs for this cell
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				// i-th potential basis functions at the point q
				const double psi_i_potential = 
										scratch.Poisson_fe_values[Potential].value(i,q);
		
				// get the local RHS values for this cell
				data.local_rhs(i) += -psi_i_potential 
															* scratch.Poisson_rhs_values[q]			
													 		* scratch.Poisson_fe_values.JxW(q);
			} // for i
		} // for q

		// loop over all the faces of this cell to calculate the vector
		// from the dirichlet boundary conditions if the face is on the 
		// Dirichlet portion of the boundary
		for(unsigned int face_no=0;
					face_no<GeometryInfo<dim>::faces_per_cell;
					face_no++)
		{
			// obtain the face_no-th face of this cell
			typename DoFHandler<dim>::face_iterator 	face = cell->face(face_no);

			// apply Dirichlet boundary conditions
			if((face->at_boundary()) && (face->boundary_id() == Dirichlet))
			{	
				// get the values of the shape functions at this boundary face
				scratch.Poisson_fe_face_values.reinit(cell,face_no);
			
				// get the values of the dirichlet boundary conditions evaluated
				// on the quadrature points of this face

				test_Poisson_bc.value_list(
													scratch.Poisson_fe_face_values.get_quadrature_points(),
													scratch.Poisson_bc_values);
		
				// loop over all the quadrature points of this face
				for(unsigned int q=0; q<n_face_q_points; q++)
				{
					// loop over all the test function dofs of this face
					for(unsigned int i=0; i<dofs_per_cell; i++)
					{

						// - \int_{face} p * n * (phi_{Dichlet}) dx
						data.local_rhs(i)	+=	
													-(scratch.Poisson_fe_face_values[VectorField].value(i,q) *
														scratch.Poisson_fe_face_values.normal_vector(q) *
														scratch.Poisson_bc_values[q] *
														scratch.Poisson_fe_face_values.JxW(q));
					} // for i
				} // for q
			} // end if
		} // end for face_no
	} // assemble_local_test_rhs

	template<int dim>
	void
	MixedFEM<dim>::
	compute_errors(const Triangulation<dim>			& triangulation,
								 DoFHandler<dim>							& Poisson_dof_handler,
								 Vector<double>					& solution,
								 double 											& potential_error,
								 double 											& field_error) const
	{
		const ComponentSelectFunction<dim> potential_mask(dim, dim+1);
		const ComponentSelectFunction<dim> 
												vectorField_mask(std::make_pair(0,dim), dim+1);

		unsigned int degree = Poisson_dof_handler.get_fe().degree;
		unsigned int n_cells = triangulation.n_active_cells();

		QTrapez<1>				q_trapez;
		QIterated<dim> 		quadrature(q_trapez, degree+2);
		Vector<double> 		cellwise_errors(n_cells);


		VectorTools::integrate_difference(Poisson_dof_handler, 
																			solution, 
																			test_Poisson_solution,
																			 cellwise_errors, quadrature, 
																			 VectorTools::L2_norm,
																			 &potential_mask);

		potential_error = cellwise_errors.l2_norm();

		VectorTools::integrate_difference(Poisson_dof_handler, 
																			solution, 
																			test_Poisson_solution,
																			 cellwise_errors, quadrature, 
																			 VectorTools::L2_norm,
																			 &vectorField_mask);

		field_error = cellwise_errors.l2_norm();
		
	
	}	
	
	template<int dim>
	void
	MixedFEM<dim>::
	output_rescaled_results(DoFHandler<dim>	& dof_handler,
								const Vector<double> & solution,
								const ParameterSpace::Parameters & sim_params,
								const unsigned int time_step_number) const
	{
		PostProcessor<dim> 	postprocessor(sim_params,
																			false,
																			"noname");
		DataOut<dim>	data_out;
		data_out.attach_dof_handler(dof_handler);
		data_out.add_data_vector(solution, 
														 postprocessor);
		
		data_out.build_patches();
		std::string file = "Poisson-"
												+ 
												Utilities::int_to_string(time_step_number,3) 
												+ 
												 ".vtu";

		std::ofstream output(file.c_str());
		data_out.write_vtu(output);
		output.close();
	}


	template<int dim>
	void
	MixedFEM<dim>::
	output_unscaled_results(DoFHandler<dim>	& dof_handler,
								const Vector<double> & solution,
								const unsigned int time_step_number) const
	{
		std::vector<std::string> solution_names;

		for(unsigned int d=0; d<dim; d++)
			solution_names.push_back("Field");	
		solution_names.push_back("Potential");	

		std::vector<DataComponentInterpretation::DataComponentInterpretation>
		component_interpretation(dim,
														 DataComponentInterpretation::component_is_part_of_vector);

		component_interpretation.push_back( 
											DataComponentInterpretation::component_is_scalar);
			
		DataOut<dim>	data_out;
		data_out.attach_dof_handler(dof_handler);
		data_out.add_data_vector(solution, 
														solution_names,
														DataOut<dim>::type_dof_data,
														component_interpretation);

		data_out.build_patches();
		std::string file = "Poisson-"
												+ 
												Utilities::int_to_string(time_step_number,3) 
												+ 
												 ".vtu";

		std::ofstream output(file.c_str());
		data_out.write_vtu(output);
		output.close();
	}



} // end namespace
