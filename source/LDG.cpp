#include "../include/LDG.hpp"
//#include "PostProcessor.cpp"

// TO MAKE A STEADY STATE SOLVER: Replace 1.0/delta_t -> 0.0/delta_t

namespace LDG_System
{
	using namespace std;

	// LDG FEM constructor
	template<int dim>
	LDG<dim>::
	LDG()
	:
	test_Poisson_solution(),
	test_Poisson_bc(),
	test_Poisson_rhs(),
	test_LDG_rhs(),
	test_LDG_bc(),
	test_LDG_interface(),
	test_LDG_solution(),	//: instantiate function objects
	test_DD_rhs(),
	test_DD_bc(),
	test_DD_solution(),	//: instantiate function objects
	test_interface_rhs(),
	test_interface_bc(),
	test_interface_function(),
	test_interface_solution(),	//: instantiate function objects
	test_interface_initial()
	{
	}

	// LDG FEM destructor
	template<int dim>
	LDG<dim>::
	~LDG()
	{
	}

	template<int dim>
	void
	LDG<dim>::
	assemble_local_LDG_mass_matrix(
				const typename DoFHandler<dim>::active_cell_iterator & cell,
				Assembly::AssemblyScratch<dim>			& scratch,
				Assembly::DriftDiffusion::CopyData<dim>		& data,
				const double	 				& delta_t)
	{
		const unsigned int dofs_per_cell 	= scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points		= scratch.carrier_fe_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		const FEValuesExtractors::Scalar Density(dim);

		// reinitialize everything for this cell
		scratch.carrier_fe_values.reinit(cell);
		data.local_mass_matrix=0;
	
		// loop over all the quadrature points of this cell
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// loop over all the test function dofs on this cell
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
	
				const double psi_i_density	= 
								scratch.carrier_fe_values[Density].value(i,q);

				// loop over all the trial function dofs for this cell
				for(unsigned int j=0; j<dofs_per_cell; j++)	
				{	
					const double psi_j_density	= 
									scratch.carrier_fe_values[Density].value(j,q);
					// construct the local mass matrix
					// int_{Omega} (1/dt) * v * u dx
					data.local_mass_matrix(i,j) += 
							(1.0/delta_t) * psi_i_density * psi_j_density 
							* scratch.carrier_fe_values.JxW(q);
				}
			}
		}
	}

	template<int dim>
	void 
	LDG<dim>::
	assemble_local_LDG_cell_and_bc_terms(
				const typename DoFHandler<dim>::active_cell_iterator & cell,
				Assembly::AssemblyScratch<dim>			 & scratch,
				Assembly::DriftDiffusion::CopyData<dim>		 & data,
				const double					 & scaled_mobility_1,
				const double 					 & scaled_mobility_2,
				const double 					 & delta_t,
				const double 					 & transient_or_steady,
				const double					 & penalty)
	{
		// NOTE: this has been called by assemble_carrier_system so it is local to a cell
		//  Assembles the body intergral matrix terms as well as the flux matrices for	
		// the boundary conditions

		const unsigned int dofs_per_cell   = scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points	   = scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points =	
						scratch.carrier_fe_face_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);

		// reinitialize everything for this cell
		scratch.carrier_fe_values.reinit(cell);
		data.local_matrix_1=0;
		data.local_matrix_2=0;
		
		double h = cell->diameter();
	
		// loop over all the quadrature points of this cell
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// loop over all the test function dofs on this cell
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				// get the test functions for this cell at quadrature point q
				const Tensor<1, dim>  psi_i_field	= 
									scratch.carrier_fe_values[Current].value(i,q);
				const double 	      div_psi_i_field	= 
									scratch.carrier_fe_values[Current].divergence(i,q);
				const Tensor<1, dim>  grad_psi_i_density = 
									scratch.carrier_fe_values[Density].gradient(i,q);
				const double 	      psi_i_density	= 
									scratch.carrier_fe_values[Density].value(i,q);	
				

				// loop over all the trial function dofs for this cell
				for(unsigned int j=0; j<dofs_per_cell; j++)	
				{	
					// get the trial functions for this cell at the quadrature point 1
					const Tensor<1, dim>	psi_j_field	  = 
									scratch.carrier_fe_values[Current].value(j,q);
					const double 		psi_j_density = 
									scratch.carrier_fe_values[Density].value(j,q);	
				
			// construct the local LDG stiffness matrix i.e. all the solid integrals
			// int_{Omega} ((t_or_s/dt) * v * u + p * \mu^{-1} * q - div(p) * u - grad(v) * q) dx
					data.local_matrix_1(i,j)  +=  
							(
							((transient_or_steady/delta_t) 
							* psi_i_density * psi_j_density )
							+
							( psi_i_field * (1.0/scaled_mobility_1) * psi_j_field)
							- 
							(div_psi_i_field * psi_j_density)
							- 
							(grad_psi_i_density * psi_j_field)
							) 
							* scratch.carrier_fe_values.JxW(q);
	
					data.local_matrix_2(i,j)  +=  
							(
							((transient_or_steady/delta_t) 
							* psi_i_density * psi_j_density )
							+
							( psi_i_field * (1.0/scaled_mobility_2) * psi_j_field)
							- 
							(div_psi_i_field * psi_j_density)
							- 
							(grad_psi_i_density * psi_j_field)
							) 
							* scratch.carrier_fe_values.JxW(q);
					
					

				} // for j
			} // for i
		}	// for q
		
		// loop over all the faces of this cell to see which one are 
		// on the boundary and calculate the local flux matrices corresponding to 
		// 1.) Dirichlet boundary conditions
		// 2.) Neumann boundary conditions
		for(unsigned int face_no=0; face_no< GeometryInfo<dim>::faces_per_cell; face_no++)
		{
			// get the face_no-th face of this cell
			typename DoFHandler<dim>::face_iterator 	face = cell->face(face_no);
				
			// test to see if it is as the boundary
			if(face->at_boundary() )
			{
				// reinitialize fe_face_values for this face and then 
				// if compute corresponding boundary condition matrix
				scratch.carrier_fe_face_values.reinit(cell, face_no);

				// construct the dirichlet matrix
				if(face->boundary_id() == Dirichlet)
				{
					// loop over alll the quadrature points of this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// loop over test function dofs of this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// get the test function
							const double 	psi_i_density	= 
									scratch.carrier_fe_face_values[Density].value(i,q);	
 	
							// loop over all trial function dofs of this face
							for(unsigned int j=0; j<dofs_per_cell; j++)	
							{	
								// get the trial function
								const Tensor<1, dim>  psi_j_field	= 
									scratch.carrier_fe_face_values[Current].value(j,q);
								// get the test function
								const double psi_j_density	= 
									scratch.carrier_fe_face_values[Density].value(j,q);	
 

								// int_{\Gamma_{D}} v^{-} n^{-} q^{-} ds
								data.local_matrix_1(i,j) += psi_i_density *
										(scratch.carrier_fe_face_values.normal_vector(q) *
										psi_j_field 
										+
										(penalty/h) * psi_j_density) * 
										scratch.carrier_fe_face_values.JxW(q);				

								data.local_matrix_2(i,j) += psi_i_density * 
										(scratch.carrier_fe_face_values.normal_vector(q) *
										psi_j_field 
										+
										(penalty/h) *	psi_j_density)* 
										scratch.carrier_fe_face_values.JxW(q);				


							} // end for j
						} // end for i
					} // for q
				} // if Dirichlet
				else if((face->boundary_id() == Interface) ||
						(face->boundary_id() == Neumann)   ||
						(face->boundary_id() == Schottky) )
				{
					// loop over alll the quadrature points of this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// loop over test function dofs of this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
 							// get the test function
							const Tensor<1, dim>  psi_i_field	=
										scratch.carrier_fe_face_values[Current].value(i,q);
 
							// loop over all the trial function dofs of this face
							for(unsigned int j=0; j<dofs_per_cell; j++)	
							{	
								// get the trial function
								const double 	psi_j_density = 
										scratch.carrier_fe_face_values[Density].value(j,q);	

								// int_{\Gamma_{N}} p^{-} n^{-} u^{-} ds
								data.local_matrix_1(i,j) += 
										psi_i_field * 
										scratch.carrier_fe_face_values.normal_vector(q) *
										psi_j_density * 
										scratch.carrier_fe_face_values.JxW(q);				
	
								data.local_matrix_2(i,j) += 
										psi_i_field * 
										scratch.carrier_fe_face_values.normal_vector(q) *
										psi_j_density * 
										scratch.carrier_fe_face_values.JxW(q);				
							} // end j
						} // end i
					} // end q
				} // end Neumann
				else
					Assert(false, ExcNotImplemented() ); // no other boundary terms
			} //if on boundary

		} // for face_no
	} // assemble cell
	
	template<int dim>
	void 
	LDG<dim>::
	assemble_flux_terms(DoFHandler<dim>			& carrier_dof_handler,
			ChargeCarrierSpace::CarrierPair<dim>	& carrier_pair,
			FESystem<dim>				& Poisson_fe,
			FESystem<dim>				& carrier_fe)
										
	{
		//////////////////////////////////////////////////////////////////////////
		// Sequential Assembly of LDG flux matrices over all the interior 
		// faces 
		//////////////////////////////////////////////////////////////////////////
		Assembly::AssemblyScratch<dim>	scratch(Poisson_fe,
							carrier_fe,
							QGauss<dim>(carrier_fe.degree+2),
							QGauss<dim-1>(carrier_fe.degree+2));


		Assembly::DriftDiffusion::CopyData<dim>	data(carrier_fe);

		typename DoFHandler<dim>::active_cell_iterator
		cell = carrier_dof_handler.begin_active(),
		endc = carrier_dof_handler.end();

		// loop over all cells 
		for(; cell != endc; cell++)
		{
			// get the map for the local dofs to global dofs for this cell
			cell->get_dof_indices(data.local_dof_indices);

			// loop over all the faces of this cell to calculate
			// the local flux matrices corresponding to the LDG central fluxes
			for(unsigned int face_no=0; 
					face_no< GeometryInfo<dim>::faces_per_cell; 
					face_no++)
			{
				// get the face_no-th face of this cell
				typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
		
				// make sure that this face is an interior face
				if( !(face->at_boundary()) )
				{
					// now we are on the interior face elements and we want to make
					// sure that the neighbor cell to this cell is a valid cell
					Assert(cell->neighbor(face_no).state() == IteratorState::valid,
										   ExcInternalError());

					// get the neighbor cell that is adjacent to this cell's face
					typename DoFHandler<dim>::cell_iterator neighbor = 
										cell->neighbor(face_no);

					// if this face has children (more refined faces) then 
					// the neighbor cell to this cell is more refined than
					// this cell is so we have to deal with that case
					if(face->has_children())
					{
						// get the face such that 
						// neighbor->face(neighbor_face_no) = cell->face(face_no)
						const unsigned int neighbor_face_no = 
										cell->neighbor_of_neighbor(face_no);
	
						// loop over all the subfaces of this face
						for(unsigned int subface_no=0;
							 	subface_no < face->number_of_children();
								subface_no++)
						{
							// get the refined neighbor cell that matches this 
							// face and subface number
							typename DoFHandler<dim>::cell_iterator	neighbor_child =
									cell->neighbor_child_on_subface(face_no, subface_no);

							// parent cant be more than one refinement level above 
							// the child
							Assert(!neighbor_child->has_children(), ExcInternalError());

							// reinitialize the fe_subface_values to this cell's subface and
							// neighbor_childs fe_face_values to its face
							scratch.carrier_fe_subface_values.reinit(cell, face_no, subface_no);
							scratch.carrier_fe_neighbor_face_values.reinit(neighbor_child, 
												   neighbor_face_no);
						
							// get the map for the local dofs to global dofs for the neighbor
							neighbor_child->get_dof_indices(data.local_neighbor_dof_indices);					
	
							// get the smaller of the h's
							double h = std::min(cell->diameter(), 
												neighbor->diameter());

							assemble_local_flux_terms(scratch,
											data,
											(carrier_pair.penalty/h) );
					
							// now add the local ldg flux matrices to the global one
							// for the carrier_1s and carrier_2s. 
							distribute_local_fluxes_to_global(carrier_pair,
															  data);


							
						} // for subface_no	
					} // if face has children
					else
					{
						// we now know that the neighbor cell of this cell's face 
						// is on the the same refinement level and therefore
						// cell with the lower index does the work
						if((neighbor->level() == cell->level()) &&
							(neighbor->index() > cell->index()) ) 
						{
						  // get the face of the nighbor such that
						  // neighbor->face(neighbor_face_no) = cell->face(face_no)
						 	const unsigned int neighbor_face_no = 
										cell->neighbor_of_neighbor(face_no);

							// reinitialize the fe_face_values on their respective face
							scratch.carrier_fe_face_values.reinit(cell, face_no);
							scratch.carrier_fe_neighbor_face_values.reinit(neighbor, 
												 neighbor_face_no);
						
							// get the map for the local dofs to global dofs for the neighbor
							neighbor->get_dof_indices(data.local_neighbor_dof_indices);					
	
							// assmble the local LDG flux matrices for this face using 
							// the assemblr object
							
							// get the smaller of the h's
							double h = std::min(cell->diameter(), 
												neighbor->diameter());

							assemble_local_flux_terms(scratch, 
													data,
													(carrier_pair.penalty/h));

							// now add the local ldg flux matrices to the global one
							// for the carrier_1s and carrier_2s. 
							distribute_local_fluxes_to_global(carrier_pair,
															  data);
						}	// end if index() >
					} // else cell not have children
				} // end if interior
			}	// end face_no
		} // for cell
	} // asssemble_flux_terms
	

	template<int dim>	
	void 
	LDG<dim>::
	assemble_local_flux_terms(
				Assembly::AssemblyScratch<dim>		 	 & scratch,
				Assembly::DriftDiffusion::CopyData<dim>	 & data,
				const double 				 			 & penalty)
	{
		// this has been called from a cells face and constructs the local ldg flux
		// matrices across that face
		const unsigned int n_face_points  =	
					scratch.carrier_fe_face_values.n_quadrature_points;
	 	const unsigned int dofs_this_cell = 
					scratch.carrier_fe_face_values.dofs_per_cell;
		const unsigned int dofs_neighbor_cell =
					scratch.carrier_fe_neighbor_face_values.dofs_per_cell;	

		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);

		// reset the local LDG flux matrices to zero
		data.vi_ui_matrix = 0;
		data.vi_ue_matrix = 0;
		data.ve_ui_matrix = 0;
		data.ve_ue_matrix = 0;

		Point<dim>	beta;
		for(unsigned int d=0; d<dim; d++)
			beta(d) = 1.0;
		beta /= sqrt(beta.square());

			// loop over all the quadrature points on this face	
		for(unsigned int q=0; q<n_face_points; q++)
		{
			// loop over all the test functiion dofs of this face
			// and get the test function values at this quadrature point
			for(unsigned int i=0; i<dofs_this_cell; i++)
			{	
				const Tensor<1,dim>  psi_i_field_minus	 = 
								scratch.carrier_fe_face_values[Current].value(i,q);
				const double	    psi_i_density_minus = 
								scratch.carrier_fe_face_values[Density].value(i,q);
		
				// loop over all the trial function dofs of this face
				for(unsigned int j=0; j<dofs_this_cell; j++)
				{
					// loop over all the trial functiion dofs of this face
					// and get the trial function values at this quadrature point

					const Tensor<1,dim>	psi_j_field_minus		= 
									scratch.carrier_fe_face_values[Current].value(j,q);
					const double 		psi_j_density_minus		=
									scratch.carrier_fe_face_values[Density].value(j,q);
	
					// int_{face} n^{-} * ( p_{i}^{-} u_{j}^{-} + v^{-} q^{-} ) dx
					// 					  + penalty v^{-}u^{-} dx
					data.vi_ui_matrix(i,j)	+= (
									 0.5 * (
									psi_i_field_minus * 
									scratch.carrier_fe_face_values.normal_vector(q) *
									psi_j_density_minus
									+ 
									psi_i_density_minus *
									scratch.carrier_fe_face_values.normal_vector(q) *
									psi_j_field_minus )
									+ 
									beta *
									psi_i_field_minus *
									psi_j_density_minus
									-
									beta *
									psi_i_density_minus *
								 	psi_j_field_minus 
									+ 
									penalty * 
									psi_i_density_minus *
									psi_j_density_minus
									) * 
									scratch.carrier_fe_face_values.JxW(q);				
				} // for j
			
				for(unsigned int j=0; j<dofs_neighbor_cell; j++)
				{
					const Tensor<1, dim>	psi_j_field_plus	= 
								scratch.carrier_fe_neighbor_face_values[Current].value(j,q);
					const double 			psi_j_density_plus		=
								scratch.carrier_fe_neighbor_face_values[Density].value(j,q);
							
					// int_{face} n^{-} * ( p_{i}^{-} u_{j}^{+} + v^{-} q^{+} ) dx
					// 					  - penalty v^{-}u^{+} dx
					data.vi_ue_matrix(i,j) += (	
									0.5 * (
									psi_i_field_minus * 
									scratch.carrier_fe_face_values.normal_vector(q) *
									psi_j_density_plus 
									+ 
									psi_i_density_minus *
									scratch.carrier_fe_face_values.normal_vector(q) *
									psi_j_field_plus ) 
									-
		 							beta *
									psi_i_field_minus *
									psi_j_density_plus
									+
									beta *
									psi_i_density_minus *
									psi_j_field_plus
									-
		 	 						penalty * 
									psi_i_density_minus *
									psi_j_density_plus
									) *
								 	scratch.carrier_fe_face_values.JxW(q);				
				} // for j
			} // for i

			for(unsigned int i=0; i<dofs_neighbor_cell; i++)
			{
				const Tensor<1,dim>  psi_i_field_plus = 
							scratch.carrier_fe_neighbor_face_values[Current].value(i,q);
				const double		 psi_i_density_plus = 
							scratch.carrier_fe_neighbor_face_values[Density].value(i,q);

				for(unsigned int j=0; j<dofs_this_cell; j++)
				{
					const Tensor<1, dim>	psi_j_field_minus	= 
									scratch.carrier_fe_face_values[Current].value(j,q);
					const double 			psi_j_density_minus		=
									scratch.carrier_fe_face_values[Density].value(j,q);

					// int_{face} -n^{-} * ( p_{i}^{+} u_{j}^{-} + v^{+} q^{-} )
					// 					  - penalty v^{+}u^{-} dx
				
					data.ve_ui_matrix(i,j) +=	( 
									-0.5 * (
									psi_i_field_plus * 
									scratch.carrier_fe_face_values.normal_vector(q) *
									psi_j_density_minus 
									+
									psi_i_density_plus *
									scratch.carrier_fe_face_values.normal_vector(q) *
									psi_j_field_minus)
									-
									beta *
									psi_i_field_plus *
									psi_j_density_minus
									+
									beta *
									psi_i_density_plus *
									psi_j_field_minus 
									-
									penalty * 
									psi_i_density_plus *
									psi_j_density_minus
									) *
									scratch.carrier_fe_face_values.JxW(q);				
				} // for j
			
				for(unsigned int j=0; j<dofs_neighbor_cell; j++)
				{
					const Tensor<1, dim>	psi_j_field_plus	= 
								scratch.carrier_fe_neighbor_face_values[Current].value(j,q);
					const double 			psi_j_density_plus		=
								scratch.carrier_fe_neighbor_face_values[Density].value(j,q);
						
					// int_{face} -n^{-} * ( p_{i}^{+} u_{j}^{+} + v^{+} q^{+} )
					// 					  + penalty v^{+}u^{+} dx
					data.ve_ue_matrix(i,j) +=	( 
									-0.5 * (
									psi_i_field_plus * 
									scratch.carrier_fe_face_values.normal_vector(q) *
									psi_j_density_plus 
									+
									psi_i_density_plus *
									scratch.carrier_fe_face_values.normal_vector(q) *
									psi_j_field_plus ) 
									+
									beta *
									psi_i_field_plus *
									psi_j_density_plus
									-
									beta *
									psi_i_density_plus *
									psi_j_field_plus 
									+
									penalty * 
									psi_i_density_plus *
									psi_j_density_plus
									) *
									scratch.carrier_fe_face_values.JxW(q);				
				} // for j
			} // for i
		} // for q
	} // end assemble_flux_terms() 

	template<int dim>
	void
	LDG<dim>::
	distribute_local_fluxes_to_global(
					ChargeCarrierSpace::CarrierPair<dim>  & carrier_pair,
					Assembly::DriftDiffusion::CopyData<dim>	& data)
	{

		// NOTE: There are the same for carrier_1s and carrier_2s 
		// only the A matrix will be changed by scaling 
		carrier_pair.constraints.distribute_local_to_global(
											data.vi_ui_matrix,
											data.local_dof_indices,
											carrier_pair.carrier_1.system_matrix);
							
		carrier_pair.constraints.distribute_local_to_global(
											data.vi_ue_matrix,
											data.local_dof_indices,
											data.local_neighbor_dof_indices,
											carrier_pair.carrier_1.system_matrix);
	
		carrier_pair.constraints.distribute_local_to_global(
											data.ve_ui_matrix,
											data.local_neighbor_dof_indices,
											data.local_dof_indices,
											carrier_pair.carrier_1.system_matrix);

		carrier_pair.constraints.distribute_local_to_global(
											data.ve_ue_matrix,
											data.local_neighbor_dof_indices,
											carrier_pair.carrier_1.system_matrix);

		carrier_pair.constraints.distribute_local_to_global(
											data.vi_ui_matrix,
											data.local_dof_indices,
											carrier_pair.carrier_2.system_matrix);
		
		carrier_pair.constraints.distribute_local_to_global(
											data.vi_ue_matrix,
											data.local_dof_indices,
											data.local_neighbor_dof_indices,
											carrier_pair.carrier_2.system_matrix);
				
		carrier_pair.constraints.distribute_local_to_global(
											data.ve_ui_matrix,
											data.local_neighbor_dof_indices,
											data.local_dof_indices,
											carrier_pair.carrier_2.system_matrix);

		carrier_pair.constraints.distribute_local_to_global(
											data.ve_ue_matrix,
											data.local_neighbor_dof_indices,
											carrier_pair.carrier_2.system_matrix);

	}


	template<int dim>
	void
	LDG<dim>::
	assemble_local_test_rhs(
				const typename DoFHandler<dim>::active_cell_iterator & cell,
				Assembly::AssemblyScratch<dim>			 & scratch,
				Assembly::DriftDiffusion::CopyData<dim>		 & data,
				const double					 & penalty)
	{
		
		// this assembles the drift term in the ldg formulation.  it uses the electric field
		// at the current iteration and the density of the carrier at the previous time step
		const unsigned int dofs_per_cell 	=  scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points		=  scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points	=  scratch.carrier_fe_face_values.n_quadrature_points;


		cell->get_dof_indices(data.local_dof_indices);

		// reinitialize the fe_values 
		scratch.carrier_fe_values.reinit(cell);

		// reset the local_rhs to be zero
		data.local_carrier_1_rhs=0;
		data.local_carrier_2_rhs=0;
	
		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);

		// the test version	
		test_Poisson_rhs.value_list(scratch.carrier_fe_values.get_quadrature_points(),
					scratch.generation_values);

		double h = cell->diameter();

		
		// loop over all the quadrature points in this cell and compute body integrals
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// loop over all the test function dofs and get the test functions
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				const double	psi_i_density = scratch.carrier_fe_values[Density].value(i,q);


				// contribution from RHS function + Drift
				// int_{Omega} v * R dx
				data.local_carrier_1_rhs(i) += ( 
							(psi_i_density * scratch.generation_values[q])
							) *
							scratch.carrier_fe_values.JxW(q);

				data.local_carrier_2_rhs(i) += ( 
							(psi_i_density * scratch.generation_values[q])
							) *
							scratch.carrier_fe_values.JxW(q);

			} // for i
		}	// for q

		// loop over all the faces of this cell and compute the contribution from the 
		// boundary conditions
		for(unsigned int face_no=0; 
				face_no< GeometryInfo<dim>::faces_per_cell;
				face_no++)
		{
			// get the face_no-th face of this cell
			typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
			
			// if on boundary apply boundayr conditions	
			if(face->at_boundary() )
			{
				if(face->boundary_id() == Dirichlet)
				{
					// reinitialize the fe_face_values for this cell ONLY if it is as the
					//boundary otherwise its a waste.  then assemble the appropriate
					//boundary conditions
					scratch.carrier_fe_face_values.reinit(cell, face_no);

					// test version		
					test_Poisson_bc.value_list(
								scratch.carrier_fe_face_values.get_quadrature_points(),
								scratch.carrier_1_bc_values,
								dim); // calls the density values of the donor profile
									// not the current ones
					// test version		
					test_Poisson_bc.value_list(
								scratch.carrier_fe_face_values.get_quadrature_points(),
								scratch.carrier_2_bc_values,
								dim); // calls the density values of the donor profile
									// not the current ones


				// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// get the test function
							const Tensor<1, dim>  psi_i_field	= 
										scratch.carrier_fe_face_values[Current].value(i,q);
 
							// get the test function
							const double psi_i_density = 
										scratch.carrier_fe_face_values[Density].value(i,q);

							// int_{\Gamma_{D}} -p^{-} n^{-} u_{D} + penalty v * u_{D} ds
							data.local_carrier_1_rhs(i) += 
										(-1.0 * psi_i_field *
										scratch.carrier_fe_face_values.normal_vector(q) 
										+ 
										(penalty/h) *
										psi_i_density) *
										scratch.carrier_1_bc_values[q] *
										scratch.carrier_fe_face_values.JxW(q);				
					
							data.local_carrier_2_rhs(i) +=  
										(-1.0 * psi_i_field *
										scratch.carrier_fe_face_values.normal_vector(q) *
										+ 
										(penalty/h) *
										psi_i_density) *
										scratch.carrier_2_bc_values[q] *
										scratch.carrier_fe_face_values.JxW(q);				


						} // for i
					}	// for q
				} // end Dirichlet
				else if(face->boundary_id() == Neumann)
				{
					// NOTHIN TO DO IF INSULATING
				}
				else
					Assert(false, ExcNotImplemented() );
			} // end at boundary
		} // end for face_no
	}

	

	template<int dim>
	void
	LDG<dim>::
	assemble_local_test_transient_rhs(
					const typename DoFHandler<dim>::active_cell_iterator & cell,
					Assembly::AssemblyScratch<dim>			 & scratch,
					Assembly::DriftDiffusion::CopyData<dim>	 	 & data,				
					const Vector<double>				 & old_solution,
					const double 					 & time,
					const double					 & penalty)			
	{

		// this assembles the drift term in the ldg formulation.  it uses the 
		// electric field at the current iteration and the density of the 
		// carrier at the previous time step
		const unsigned int dofs_per_cell	 = 
						scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points		 = 
						scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points	 = 
						scratch.carrier_fe_face_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		// reinitialize the fe_values 
		scratch.carrier_fe_values.reinit(cell);
		
		// reset the local_rhs to be zero
		data.local_carrier_1_rhs=0;
		data.local_carrier_2_rhs=0;
	
		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);
		Tensor<1,dim>		field_values;
		field_values[0]	=	1.0;
		field_values[1]	=	0.0;
	
		test_LDG_rhs.set_time(time);
		test_LDG_interface.set_time(time);
		test_LDG_bc.set_time(time);

		// the test version	
		test_LDG_rhs.value_list(scratch.carrier_fe_values.get_quadrature_points(),
					scratch.generation_values);
	
		// get the values of carrier_1 and carrier_2 densities at the pevious time step
		scratch.carrier_fe_values[Density].get_function_values(old_solution,
								scratch.old_carrier_1_density_values);

		double h = cell->diameter();
		// loop over all the quadrature points in this cell and compute body integrals
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// loop over all the test function dofs and get the test functions
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				const double		psi_i_density	= 
								scratch.carrier_fe_values[Density].value(i,q);

				const Tensor<1,dim>	psi_i_field	=
								scratch.carrier_fe_values[Current].value(i,q);
					
				// contribution from RHS function + Drift
				// int_{Omega} v * R +  p * E * u dx
				data.local_carrier_1_rhs(i) += ( 
								(psi_i_density * scratch.generation_values[q])
								-
								psi_i_field *
								field_values *
								scratch.old_carrier_1_density_values[q]
								) *
								scratch.carrier_fe_values.JxW(q);

			} // for i
		}	// for q

		// loop over all the faces of this cell and compute the contribution 
		// from the boundary conditions
		for(unsigned int face_no=0; 
				face_no< GeometryInfo<dim>::faces_per_cell;
				face_no++)
		{
			// get the face_no-th face of this cell
			typename DoFHandler<dim>::face_iterator 	face = cell->face(face_no);
			
			// if on boundary apply boundayr conditions	
			if(face->at_boundary() )
			{
				// reinitialize the fe_face_values for this cell ONLY if it is as the
				//boundary otherwise its a waste.  then assemble the appropriate
				//boundary conditions
				scratch.carrier_fe_face_values.reinit(cell, face_no);
			
				if(face->boundary_id() == Dirichlet)
				{
					// get the density
					test_LDG_bc.value_list(
							scratch.carrier_fe_face_values.get_quadrature_points(),
							scratch.carrier_1_bc_values);
							//			dim);

					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// get the test function
							const Tensor<1, dim>  psi_i_field	= 
									scratch.carrier_fe_face_values[Current].value(i,q);
							const double 		  psi_i_density = 
									scratch.carrier_fe_face_values[Density].value(i,q);
 
							// int_{\Gamma_{D}} -p^{-} n^{-} u_{D} ds
							data.local_carrier_1_rhs(i) += 
									(-1.0 * psi_i_field *
									scratch.carrier_fe_face_values.normal_vector(q) 
									+ 
									(penalty/h) * psi_i_density) * 
									scratch.carrier_1_bc_values[q] *
									scratch.carrier_fe_face_values.JxW(q);	
						} // for i
					}	// for q
				} // end Dirichlet
				else if(face->boundary_id() == Interface)
				{
					test_LDG_interface.value_list(
							scratch.carrier_fe_face_values.get_quadrature_points(),
							scratch.carrier_1_bc_values);

					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// get the test function
							const double psi_i_density =
									scratch.carrier_fe_face_values[Density].value(i,q);
	
							// int_{\Gamm
							data.local_carrier_1_rhs(i) +=
									-1.0 * psi_i_density *
									scratch.carrier_1_bc_values[q] *
									scratch.carrier_fe_face_values.JxW(q);
		
						} // for i
					} // for q
				}
				else if(face->boundary_id() == Neumann)
				{
					// NOTHIN TO DO IF INSULATING
				}
				else
					Assert(false, ExcNotImplemented() );
			} // end at boundary
		} // for face_no
		
	
	} // for assemble_local 

	template<int dim>
	void
	LDG<dim>::
	compute_interface_errors(const Triangulation<dim>	& triangulation,
				DoFHandler<dim>			& carrier_dof_handler,
				Vector<double>			& solution,
				double 				& potential_error,
				double 				& field_error,
				const double			& time)// const
	{
		const ComponentSelectFunction<dim> potential_mask(dim, dim+1);
		const ComponentSelectFunction<dim> vectorField_mask(std::make_pair(0,dim), dim+1);
		
		unsigned int degree = carrier_dof_handler.get_fe().degree;
		unsigned int n_cells = triangulation.n_active_cells();

		QTrapez<1>		q_trapez;
		QIterated<dim> 		quadrature(q_trapez, degree+2);
		Vector<double> 		cellwise_errors(n_cells);

		test_interface_solution.set_time(time);
		VectorTools::integrate_difference(carrier_dof_handler, 
						solution, 
						test_interface_solution,
						cellwise_errors, quadrature, 
						VectorTools::L2_norm,
						&potential_mask);

		potential_error = cellwise_errors.l2_norm();

		test_interface_solution.set_time(time);
		VectorTools::integrate_difference(carrier_dof_handler, 
						solution, 
						test_interface_solution,
						cellwise_errors, quadrature, 
						VectorTools::L2_norm,
						&vectorField_mask);
		
		field_error = cellwise_errors.l2_norm();
	
	}

	template<int dim>
	void
	LDG<dim>::
	compute_coupled_errors(const Triangulation<dim>	& triangulation,
				DoFHandler<dim>		& carrier_dof_handler,
				Vector<double>		& solution,
				double 			& potential_error,
				double 			& field_error,
				const double		& time)// const
	{
		const ComponentSelectFunction<dim> potential_mask(dim, dim+1);
		const ComponentSelectFunction<dim> vectorField_mask(std::make_pair(0,dim), dim+1);
		
		unsigned int degree = carrier_dof_handler.get_fe().degree;
		unsigned int n_cells = triangulation.n_active_cells();

		QTrapez<1>		q_trapez;
		QIterated<dim> 		quadrature(q_trapez, degree+2);
		Vector<double> 		cellwise_errors(n_cells);

		test_DD_solution.set_time(time);

		VectorTools::integrate_difference(carrier_dof_handler, 
						solution, 
						test_DD_solution,
						cellwise_errors, quadrature, 
						VectorTools::L2_norm,
						&potential_mask);
		
		potential_error = cellwise_errors.l2_norm();

		VectorTools::integrate_difference(carrier_dof_handler, 
						solution, 
						test_DD_solution,
						cellwise_errors, quadrature, 
						VectorTools::L2_norm,
						&vectorField_mask);
	
		field_error = cellwise_errors.l2_norm();
	
	}

	template<int dim>
	void
	LDG<dim>::
	compute_errors(const Triangulation<dim>	& triangulation,
		 	DoFHandler<dim>		& carrier_dof_handler,
			Vector<double>		& solution,
			double 			& potential_error,
			double 			& field_error,
			const bool 		& steady_state,
			const double		& time)// const
	{
		const ComponentSelectFunction<dim> potential_mask(dim, dim+1);
		const ComponentSelectFunction<dim> vectorField_mask(std::make_pair(0,dim), dim+1);
		
		unsigned int degree = carrier_dof_handler.get_fe().degree;
		unsigned int n_cells = triangulation.n_active_cells();

		QTrapez<1>		q_trapez;
		QIterated<dim> 		quadrature(q_trapez, degree+2);
		Vector<double> 		cellwise_errors(n_cells);

		if(steady_state)
		{
			VectorTools::integrate_difference(carrier_dof_handler, 
							solution, 
							test_Poisson_solution,
							cellwise_errors, quadrature, 
							VectorTools::L2_norm,
							&potential_mask);
		}
		else
		{
			test_LDG_solution.set_time(time);
			VectorTools::integrate_difference(carrier_dof_handler, 
							solution, 
							test_LDG_solution,
							cellwise_errors, quadrature, 
							VectorTools::L2_norm,
							&potential_mask);
		}


		potential_error = cellwise_errors.l2_norm();

		if(steady_state)
		{
			VectorTools::integrate_difference(carrier_dof_handler, 
							solution, 
							test_Poisson_solution,
							cellwise_errors, quadrature, 
							VectorTools::L2_norm,
							&vectorField_mask);
		}
		else
		{
			test_LDG_solution.set_time(time);
			VectorTools::integrate_difference(carrier_dof_handler, 
							solution, 
							test_LDG_solution,
							cellwise_errors, quadrature, 
							VectorTools::L2_norm,
							&vectorField_mask);
		}	
		field_error = cellwise_errors.l2_norm();
	
	}

	template<int dim>
	void
	LDG<dim>::
	output_unscaled_results(DoFHandler<dim>			     & carrier_dof_handler,
				ChargeCarrierSpace::CarrierPair<dim> & carrier_pair,
				const unsigned int 		   time_step_number) const
	{
		// Cell values
		std::vector<std::string> carrier_solution_names_1;
		std::string carrier_1_current = carrier_pair.carrier_1.name.c_str();
		std::string carrier_1_density= carrier_pair.carrier_1.name.c_str();
		carrier_1_current += "_Current";		
		carrier_1_density += "_Density";

		std::vector<std::string> carrier_solution_names_2;
		std::string carrier_2_current = carrier_pair.carrier_2.name.c_str();
		std::string carrier_2_density= carrier_pair.carrier_2.name.c_str();
		carrier_2_current += "_Current";		
		carrier_2_density += "_Density";

		for(unsigned int d=0; d<dim; d++)
			carrier_solution_names_1.push_back(carrier_1_current.c_str());	
		carrier_solution_names_1.push_back(carrier_1_density.c_str());	

		for(unsigned int d=0; d<dim; d++)
			carrier_solution_names_2.push_back(carrier_2_current.c_str());	
		carrier_solution_names_2.push_back(carrier_2_density.c_str());	


		std::vector<DataComponentInterpretation::DataComponentInterpretation>
		component_interpretation(dim,
					DataComponentInterpretation::component_is_part_of_vector);

		component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
			
		DataOut<dim>	data_out;
		data_out.attach_dof_handler(carrier_dof_handler);

		data_out.add_data_vector(carrier_pair.carrier_1.solution, 
					carrier_solution_names_1,
					DataOut<dim>::type_dof_data,
					component_interpretation);

		data_out.add_data_vector(carrier_pair.carrier_2.solution, 
					carrier_solution_names_2,
					DataOut<dim>::type_dof_data,
					component_interpretation);

		data_out.build_patches();;

		std::string carrier_1_file = carrier_pair.material_name.c_str();
		carrier_1_file += Utilities::int_to_string(time_step_number,3);
		carrier_1_file += ".vtu";

		std::ofstream output(carrier_1_file.c_str());
		data_out.write_vtu(output);
		output.close();

	}

	template<int dim>
	void
	LDG<dim>::
	output_rescaled_results(DoFHandler<dim>		             & carrier_dof_handler,
				ChargeCarrierSpace::CarrierPair<dim> & carrier_pair,
				const ParameterSpace::Parameters     & sim_params,
				const unsigned int 		   time_step_number) const
	{
		DataOut<dim>		data_out;
		PostProcessor<dim>	postprocessor_1(sim_params,
							true,
							carrier_pair.carrier_1.name.c_str());

		data_out.attach_dof_handler(carrier_dof_handler);
		data_out.add_data_vector(carrier_pair.carrier_1.solution,
					 postprocessor_1);

		data_out.build_patches();

		PostProcessor<dim>	postprocessor_2(sim_params,
							true,
							carrier_pair.carrier_2.name.c_str());

		data_out.add_data_vector(carrier_pair.carrier_2.solution,
					postprocessor_2);

		data_out.build_patches();
	

		std::string carrier_1_file = carrier_pair.material_name.c_str();
		carrier_1_file += Utilities::int_to_string(time_step_number,3);
		carrier_1_file += ".vtu";

		std::ofstream output(carrier_1_file.c_str());
		data_out.write_vtu(output);
		output.close();

	}

	template<int dim>
	void
	LDG<dim>::
	output_unscaled_results_on_boundary(DoFHandler<dim>	 	     & carrier_dof_handler,
					ChargeCarrierSpace::CarrierPair<dim> & carrier_pair,
					const unsigned int 		   time_step_number) const
	{
		// boundary values
		std::vector<std::string> carrier_solution_names_1;
		std::string carrier_1_current = carrier_pair.carrier_1.name.c_str();
		std::string carrier_1_density= carrier_pair.carrier_1.name.c_str();
		carrier_1_current += "_Current";		
		carrier_1_density += "_Density";

		std::vector<std::string> carrier_solution_names_2;
		std::string carrier_2_current = carrier_pair.carrier_2.name.c_str();
		std::string carrier_2_density= carrier_pair.carrier_2.name.c_str();
		carrier_2_current += "_Current";		
		carrier_2_density += "_Density";

		for(unsigned int d=0; d<dim; d++)
			carrier_solution_names_1.push_back(carrier_1_current.c_str());	
		carrier_solution_names_1.push_back(carrier_1_density.c_str());	

		for(unsigned int d=0; d<dim; d++)
			carrier_solution_names_2.push_back(carrier_2_current.c_str());	
		carrier_solution_names_2.push_back(carrier_2_density.c_str());	


		std::vector<DataComponentInterpretation::DataComponentInterpretation>
		component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);

		component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
			
		DataOutFaces<dim>	data_out;
		data_out.add_data_vector(carrier_dof_handler,
					carrier_pair.carrier_1.solution, 
					carrier_solution_names_1,
					component_interpretation);

		data_out.add_data_vector(carrier_dof_handler,
					carrier_pair.carrier_2.solution, 
					carrier_solution_names_2,
					component_interpretation);

		data_out.build_patches();;

		std::string carrier_1_file = carrier_pair.material_name.c_str();
		carrier_1_file += "-Boundary";
		carrier_1_file += Utilities::int_to_string(time_step_number,3);
		carrier_1_file += ".vtk";

		std::ofstream output(carrier_1_file.c_str());
		data_out.write_vtk(output);
		output.close();

	}
	

} // end namespace
