#include "../include/SolarCell.hpp"
#include "Grid.cpp"
#include "Assembly.cpp"
#include "MixedFEM.cpp"
#include "LDG.cpp"
#include "Generation.cpp"
#include "InitialConditions.cpp"
#include "BiasValues.cpp"
#include "CarrierPair.cpp"

namespace SOLARCELL
{
	using namespace dealii;

	template<int dim>
	SolarCellProblem<dim>::
	SolarCellProblem(const unsigned int degree,
									 ParameterHandler & param)
	:
	degree(degree),
	prm(param),
	Poisson_dof_handler(Poisson_triangulation),
	Poisson_fe(FE_RaviartThomas<dim>(degree-1), 					 1,
						 FE_DGQ<dim>(degree-1),											 1),
	Poisson_object(),
	semiconductor_dof_handler(semiconductor_triangulation),
	electron_hole_pair(),
	electrolyte_dof_handler(electrolyte_triangulation),
	redox_pair(),
	carrier_fe(FESystem<dim>(FE_DGQ<dim>(degree), dim), 1,
													 FE_DGQ<dim>(degree), 1),
	Mixed_Assembler(),
	electrons_e(),
	holes_e(),
	reductants_e(),
	oxidants_e(),
	built_in_bias(),
	applied_bias(),
	bulk_bias(),
	schottky_bias(),
	generation()
	{
		// set the parameters
		sim_params.parse_and_scale_parameters(prm);
	
		// set the voltage bias functions
		applied_bias.set_value(sim_params.scaled_applied_bias);
		built_in_bias.set_value(sim_params.scaled_built_in_bias);
		schottky_bias.set_value(sim_params.scaled_domain_height);
		schottky_bias.set_value(sim_params.scaled_schottky_bias);

		// set the charges name, charge sign, and mobility
		electron_hole_pair.carrier_1.set_name("Electrons");
		electron_hole_pair.carrier_1.charge_number = -1.0;
		electron_hole_pair.carrier_1.scaled_mobility =
						sim_params.scaled_electron_mobility;
	
		electron_hole_pair.carrier_2.set_name("Holes");
		electron_hole_pair.carrier_2.charge_number = 1.0;
		electron_hole_pair.carrier_2.scaled_mobility =
						sim_params.scaled_hole_mobility;

		// set material name
		electron_hole_pair.set_name("Semiconductor-");

		// set the material permittivities
		electron_hole_pair.material_permittivity = 
						sim_params.semiconductor_permittivity;

		redox_pair.carrier_1.set_name("Reductants");
		redox_pair.carrier_1.charge_number = -1.0;
		redox_pair.carrier_1.scaled_mobility =
						sim_params.scaled_reductant_mobility;

		redox_pair.carrier_2.set_name("Oxidants");
		redox_pair.carrier_2.charge_number = 1.0;
		redox_pair.carrier_2.scaled_mobility =
						sim_params.scaled_oxidant_mobility;

		// set material name
		redox_pair.set_name("Electrolyte-");

		// set the material permittivities
		redox_pair.material_permittivity = 
						sim_params.electrolyte_permittivity;

		// set the generation function to be on or off
		if(sim_params.illum_or_dark)
			generation.set_illuminated_params(sim_params);
		else
			generation.set_dark_params();
	}	//SolarCellProblem

	// destructor
	template<int dim>
	SolarCellProblem<dim>::
	~SolarCellProblem()
	{
		Poisson_dof_handler.clear();
		semiconductor_dof_handler.clear();
		electrolyte_dof_handler.clear();

	} // ~SolarCellProblem 

	template<int dim>
	void 
	SolarCellProblem<dim>::
	setup_dofs()
	{
		Poisson_object.setup_dofs(Poisson_fe,
															Poisson_dof_handler);
	
		electron_hole_pair.setup_dofs(carrier_fe,
																	semiconductor_dof_handler);

		if(full_system)
		{
			redox_pair.setup_dofs(carrier_fe,
														electrolyte_dof_handler);
		}
	}	// setup+dos

	/////////////////////////////////////////////////////////////////////////////
	//														PRINT INFO
	/////////////////////////////////////////////////////////////////////////////

	template<int dim>
	void
	SolarCellProblem<dim>::
	print_sim_info()
	{
		std::cout << "---------------------------------------------------------\n"
							<< "Triangulation Info\n"
							<< "---------------------------------------------------------\n"
							<< "Number of active cells : " 
							<< Poisson_triangulation.n_active_cells()
							<< std::endl
							<< "Total number of cells: " 
							<< Poisson_triangulation.n_cells() << std::endl
							<< "h_min: "
							<< GridTools::minimal_cell_diameter(Poisson_triangulation) << std::endl
							<< "h_max: "
							<< GridTools::maximal_cell_diameter(Poisson_triangulation) << std::endl
							<< std::endl;
	
		Poisson_object.print_info(Poisson_dof_handler);
		electron_hole_pair.print_info(semiconductor_dof_handler);
		if(full_system)
			redox_pair.print_info(electrolyte_dof_handler);
	}

	/////////////////////////////////////////////////////////////////////////////
	//														MAPPINGS
	/////////////////////////////////////////////////////////////////////////////

	template<int dim>
	void
	SolarCellProblem<dim>::
	setup_mappings()
	{
	
		/*----------------------------------------------------------------------*/
		/*	Semiconductor-Electrolyte Mappings																	*/
		/*----------------------------------------------------------------------*/
		if(full_system)
		{
			// make vector of interface cells and vector of the faces of those cell's
			// face that lies on the interface as well as a vector of the center point
			// of that face
			std::vector<Point<dim>>														semi_interface_centers;	
			typename DoFHandler<dim>::active_cell_iterator
			semi_cell	=	semiconductor_dof_handler.begin_active(),
			semi_endc	=	semiconductor_dof_handler.end();
		
			for(; semi_cell != semi_endc; semi_cell++)	
			{
				for(unsigned int face_no=0;
						face_no < GeometryInfo<dim>::faces_per_cell;
						face_no++)
				{
					typename DoFHandler<dim>::face_iterator	face = semi_cell->face(face_no);
				
					if(face->boundary_id() == Interface)
					{
						semi_interface_cells.push_back(
										std::pair<unsigned int, unsigned int>(semi_cell->level(),
																									 				semi_cell->index()));
						semi_interface_faces.push_back(face_no);
						semi_interface_centers.push_back(face->center());
					}
				}
			} // end semi_cell
	
			// make vector of interface cells and vector of the faces of those cell's
			// face that lies on the interface as well as a vector of the center point
			// of that face.
			std::vector<std::pair<unsigned int, unsigned int>> temp_elec_interface_cells;
			std::vector<unsigned int>													 temp_elec_interface_faces;
			std::vector<Point<dim>>														 temp_elec_interface_centers;
	
			typename DoFHandler<dim>::active_cell_iterator
			elec_cell	=	electrolyte_dof_handler.begin_active(),
			elec_endc	=	electrolyte_dof_handler.end();

			for(; elec_cell != elec_endc; elec_cell++)	
			{
				for(unsigned int face_no=0;
						face_no < GeometryInfo<dim>::faces_per_cell;
						face_no++)
				{
					typename DoFHandler<dim>::face_iterator	face = elec_cell->face(face_no);
				
					if(face->boundary_id() == Interface)
					{
						temp_elec_interface_cells.push_back(
										std::pair<unsigned int, unsigned int>(elec_cell->level(),
																									 				elec_cell->index()));
						temp_elec_interface_faces.push_back(face_no);
						temp_elec_interface_centers.push_back(face->center());
					}
				}
			} // end elec_cell
		
			elec_interface_cells.resize(semi_interface_cells.size());
			elec_interface_faces.resize(semi_interface_faces.size());

			// find the matching cells and faces that lie along the interface.
			// We copy over this data so that cells which lie across from one another
			// are ordered in the same manner
			for(unsigned int i=0; i < semi_interface_cells.size(); i++)
			{
				for(unsigned int j=0; j < temp_elec_interface_cells.size(); j++)
				{
					if(semi_interface_centers[i].distance(temp_elec_interface_centers[j]) < 1e-13)
					{
						elec_interface_cells[i] = temp_elec_interface_cells[j];
						elec_interface_faces[i] = temp_elec_interface_faces[j];
					}
				}
			} // end i


			std::pair<std::pair<unsigned int, unsigned int>, unsigned int> semi_mapping_pair;
			std::pair<std::pair<unsigned int, unsigned int>, unsigned int> elec_mapping_pair;

			// create bijective mapping so that we can get a global interface face index
			// so that when we are on a semiconductor cell's interface face we can get the
			// corresponding face on the electrolyte triangulation from elec_interface_cell
			// and elec_interface_face
			// and vice versa
			for(unsigned int i=0; i < semi_interface_cells.size(); i++)
			{
				semi_mapping_pair.first = semi_interface_cells[i];
				semi_mapping_pair.second = i;
				semi_interface_map.insert(semi_mapping_pair);
		
				elec_mapping_pair.first = elec_interface_cells[i];
				elec_mapping_pair.second = i;
				elec_interface_map.insert(elec_mapping_pair);
			}
		
		} // if full system

		/*---------------------------------------------------------------------*/
		/* Semiconductor-Poisson Mapping and Electrolyte-Poisson Mapping			 */
		/*---------------------------------------------------------------------*/

		std::vector<std::pair<unsigned int, unsigned int>> semiconductor_cells;
		std::vector<Point<dim>> 													 semiconductor_cell_centers;
		std::vector<std::pair<unsigned int, unsigned int>> Poisson_cells;
		std::vector<Point<dim>> 													 Poisson_cell_centers;
		std::vector<std::pair<unsigned int, unsigned int>> electrolyte_cells;
		std::vector<Point<dim>>														 electrolyte_cell_centers;	

		// make mapping of Poisson cells to global poisson index							 
		typename DoFHandler<dim>::active_cell_iterator
		Poisson_cell	=	Poisson_dof_handler.begin_active(),
		Poisson_endc	=	Poisson_dof_handler.end();		
		for(; Poisson_cell != Poisson_endc; Poisson_cell++)	
		{
			Poisson_cells.push_back(
										std::pair<unsigned int, unsigned int>(Poisson_cell->level(),
																									 				Poisson_cell->index()));
			Poisson_cell_centers.push_back(Poisson_cell->center());
		
		} // end Poisson_cell
	
		// make list of semiconductor cells	
		typename DoFHandler<dim>::active_cell_iterator
		semi_cell	=	semiconductor_dof_handler.begin_active(),
		semi_endc	=	semiconductor_dof_handler.end();
		for(; semi_cell != semi_endc; semi_cell++)
		{
			semiconductor_cells.push_back(
													std::pair<unsigned int, unsigned int>(semi_cell->level(),
																												 				semi_cell->index()));
			semiconductor_cell_centers.push_back(semi_cell->center());
		}

		std::pair<unsigned int, unsigned int>					semi_pair;
		std::pair<unsigned int, unsigned int>					Poisson_pair;
		
		std::pair<std::pair<unsigned int, unsigned int>, 
							std::pair<unsigned int, unsigned int>> mapping_pair;


		// find cells which are in the semiconductor and poisson and make mapping between them
		for(unsigned int i=0; i < Poisson_cells.size(); i++)
		{
			for(unsigned int j=0; j < semiconductor_cells.size(); j++)
			{
				if(Poisson_cell_centers[i].distance(semiconductor_cell_centers[j]) < 1e-13)
				{
					semi_pair = std::pair<unsigned int, unsigned int>(
																		 semiconductor_cells[j].first, // level
																		 semiconductor_cells[j].second); // index
					
					Poisson_pair = std::pair<unsigned int, unsigned int>(
																			Poisson_cells[i].first, // level
																			Poisson_cells[i].second); // index

					// mapping from semiconductor cell to Poisson cell
					mapping_pair.first  = semi_pair;				
					mapping_pair.second = Poisson_pair;				
					s_2_p_map.insert(mapping_pair);
					
				} // if
			} // end j
		} // end i

		if(full_system)
		{
			// make list of electrolyte cells
			typename DoFHandler<dim>::active_cell_iterator
			elec_cell = electrolyte_dof_handler.begin_active(),
			elec_endc	=	electrolyte_dof_handler.end();
			for(; elec_cell != elec_endc; elec_cell++)
			{
				electrolyte_cells.push_back(
													std::pair<unsigned int, unsigned int>(elec_cell->level(),
																																elec_cell->index()));
				electrolyte_cell_centers.push_back(elec_cell->center());
			}

			std::pair<unsigned int,	unsigned int>					elec_pair;

			// find cells which are in the electrolyte and poisson and make mapping between them
			for(unsigned int i=0; i < Poisson_cells.size(); i++)
			{
				for(unsigned int j=0; j < electrolyte_cells.size(); j++)
				{
					if(Poisson_cell_centers[i].distance(electrolyte_cell_centers[j]) < 1e-13)
					{
						elec_pair = std::pair<unsigned int, unsigned int>(
																		 electrolyte_cells[j].first, // level
																		 electrolyte_cells[j].second); // index
					
						Poisson_pair = std::pair<unsigned int, unsigned int>(
																			Poisson_cells[i].first,	// level
																			Poisson_cells[i].second); // index

						// mapping from electrolyte cell to Poisson cell
						mapping_pair.first  = elec_pair;				
						mapping_pair.second = Poisson_pair;				
						e_2_p_map.insert(mapping_pair);
					
					} // if
				} // end j
			} // end i
		} // end full system

	} // setup_mappings

	////////////////////////////////////////////////////////////////////////////////
	// 											MIXED METHOD ROUTINES
	////////////////////////////////////////////////////////////////////////////////
	template<int dim>
	void 
	SolarCellProblem<dim>::
	assemble_Poisson_matrix()
	{
		// this is shown in deal.II on shared memory parallelism 
		// it builds the above matrix by distributing over multi threads locally
		// and then building up the global matrix sequentially (kinda)
		WorkStream::run(Poisson_dof_handler.begin_active(),
										Poisson_dof_handler.end(),
										std_cxx11::bind(&MixedPoisson::MixedFEM<dim>::
																		assemble_local_Poisson_matrix, 
																		Mixed_Assembler, // this object object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		sim_params.semiconductor_permittivity,
																		sim_params.electrolyte_permittivity,
																		sim_params.scaled_debeye_length), 
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_Poisson_matrix,
																		this, // this object
																		std_cxx11::_1),
										Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
										Assembly::Poisson::CopyData<dim>(Poisson_fe)
										);
	} // end assemble_Poisson_matrix

	template <int dim>
	void 
	SolarCellProblem<dim>::
	copy_local_to_global_Poisson_matrix(
											const Assembly::Poisson::CopyData<dim> & data)
	{
		// distribute local matrix to global Poisson matrix
		Poisson_object.constraints.distribute_local_to_global(data.local_matrix,
																									data.local_dof_indices,
																									Poisson_object.system_matrix);
	}	// copy_local_to_global_poisson

	template <int dim>
	void 
	SolarCellProblem<dim>::
	copy_local_to_global_Poisson_rhs(
											const Assembly::Poisson::CopyData<dim> & data)
	{
		// copy the local RHS into the global RHS for Poisson
		Poisson_object.constraints.distribute_local_to_global(data.local_rhs,
																									 data.local_dof_indices,
																									 Poisson_object.system_rhs);
	}	

	template <int dim>
	void
	SolarCellProblem<dim>::
	assemble_Poisson_rhs()
	{
		// reset rhs to zero
		Poisson_object.system_rhs = 0;

		/*-------------------------------------------------------------*/
		// This is the one coupled to the semiconductor
		/*-------------------------------------------------------------*/
		WorkStream::run(semiconductor_dof_handler.begin_active(),
										semiconductor_dof_handler.end(),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		assemble_local_Poisson_rhs_for_semiconductor,
																		this, // this object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_Poisson_rhs,
																		this, // this object
																		std_cxx11::_1),
										Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
										Assembly::Poisson::CopyData<dim>(Poisson_fe)
										);
		
		if(full_system)
		{
			/*-------------------------------------------------------------*/
			// This is the one coupled to the electrolyte
			/*-------------------------------------------------------------*/
			WorkStream::run(electrolyte_dof_handler.begin_active(),
											electrolyte_dof_handler.end(),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		assemble_local_Poisson_rhs_for_electrolyte,
																		this, // this object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_Poisson_rhs,
																		this, // this object
																		std_cxx11::_1),
										Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
										Assembly::Poisson::CopyData<dim>(Poisson_fe)
										);

		} // end if full_system
	}

	template <int dim>
	void 
	SolarCellProblem<dim>::
	assemble_local_Poisson_rhs_for_semiconductor(
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


		// Get the actual values for vector field and potential from FEValues
		// Use Extractors instead of having to deal with shapefunctions directly
		const FEValuesExtractors::Vector VectorField(0); // Vector as in Vector field
		const FEValuesExtractors::Scalar Potential(dim);
		const FEValuesExtractors::Scalar Density(dim);

		// reset the local_rhs vector to be zero
		data.local_rhs=0;

		// get the Poisson cell info
		std::pair<unsigned int, unsigned int> Poisson_cell_info = 
					s_2_p_map[std::pair<unsigned int, unsigned int>(cell->level(),
																													cell->index())];

		// get the Poisson cell
		typename DoFHandler<dim>::active_cell_iterator
								Poisson_cell(&Poisson_triangulation,
													Poisson_cell_info.first,
													Poisson_cell_info.second,
													&Poisson_dof_handler);

		Poisson_cell->get_dof_indices(data.local_dof_indices);
		scratch.Poisson_fe_values.reinit(Poisson_cell);
		// get the test rhs for poisson
		// Assemble the right hand side on semiconductor side
		scratch.carrier_fe_values.reinit(cell);
			
		// get the electron density values at the previous time step
		scratch.carrier_fe_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.old_carrier_1_density_values);

		// get the hole density values at the previous time step
		scratch.carrier_fe_values[Density].get_function_values(
																				electron_hole_pair.carrier_2.solution,
																				scratch.old_carrier_2_density_values);
	

		// get doping profiles values on this cell
		electrons_e.value_list(scratch.carrier_fe_values.get_quadrature_points(),
													scratch.donor_doping_values,
													dim); // calls the density values of the donor profile

		holes_e.value_list(scratch.carrier_fe_values.get_quadrature_points(),
											scratch.acceptor_doping_values,
											dim); // calls the density values of the donor profile

		// Loop over all the quadrature points in this cell
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// copy over the test functions
			for(unsigned int k = 0; k<dofs_per_cell; k++)
				scratch.psi_i_potential[k] = scratch.Poisson_fe_values[Potential].value(k,q);

			// loop over the test function dofs for this cell
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				// i-th potential basis functions at the point q
		
				// get the local RHS values for this cell
				// = -int_{Omega_{e}} (1/lambda^{2}) v	(N_{D} - N_{A}) - (n - p) d
				data.local_rhs(i) += -scratch.psi_i_potential[i] * //-psi_i_potential *
															(
															 (scratch.donor_doping_values[q] 
													 			- 
														 		scratch.acceptor_doping_values[q])
												    	 	+
															 (electron_hole_pair.carrier_1.charge_number *
																scratch.old_carrier_1_density_values[q]	
																+
																electron_hole_pair.carrier_2.charge_number *
																scratch.old_carrier_2_density_values[q])
														  ) * scratch.Poisson_fe_values.JxW(q);
		
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
			typename DoFHandler<dim>::face_iterator 	face = Poisson_cell->face(face_no);

			// apply Dirichlet boundary conditions.. 
			// Since we are in a semicondcutor cell we know to apply the 
			// biases
			if(face->at_boundary())
			{
				if(face->boundary_id() == Dirichlet)
				{	
					// get the values of the shape functions at this boundary face
					scratch.Poisson_fe_face_values.reinit(Poisson_cell,face_no);
			
					// get the values of the dirichlet boundary conditions evaluated
					// on the quadrature points of this face
					built_in_bias.value_list(
														scratch.Poisson_fe_face_values.get_quadrature_points(),
														scratch.Poisson_bi_values);

					applied_bias.value_list(
														scratch.Poisson_fe_face_values.get_quadrature_points(),
														scratch.Poisson_bc_values);
			
					// copy over the normal vectors
					for(unsigned int k=0; k < n_face_q_points; k++)
						scratch.normals[k] = scratch.Poisson_fe_face_values.normal_vector(k);

					// loop over all the quadrature points of this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k = 0; k < dofs_per_cell; k++)
						{
							scratch.psi_i_field[k]=	
											scratch.Poisson_fe_face_values[VectorField].value(k,q);
						}

						// loop over all the test function dofs of this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// - \int_{face} p * n * (phi_{Dichlet}) dx
							data.local_rhs(i)	+=	
													 -(scratch.psi_i_field[i] *
														 scratch.normals[q] *
														(scratch.Poisson_bi_values[q] 
														 -
														 scratch.Poisson_bc_values[q]) *
														scratch.Poisson_fe_face_values.JxW(q));
						} // for i
					} // for q
				} // end if dirichlet

				if(face->boundary_id() == Schottky)
				{	
					// get the values of the shape functions at this boundary face
					scratch.Poisson_fe_face_values.reinit(Poisson_cell,face_no);
			
					// get the values of the dirichlet boundary conditions evaluated
					// on the quadrature points of this face
					built_in_bias.value_list(
														scratch.Poisson_fe_face_values.get_quadrature_points(),
														scratch.Poisson_bi_values);
					
					schottky_bias.value_list(
														scratch.Poisson_fe_face_values.get_quadrature_points(),
														scratch.Poisson_bc_values);
			
					// copy over the normal vectors
					for(unsigned int k=0; k < n_face_q_points; k++)
						scratch.normals[k] = scratch.Poisson_fe_face_values.normal_vector(k);

					// loop over all the quadrature points of this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k = 0; k < dofs_per_cell; k++)
						{
							scratch.psi_i_field[k]=	
											scratch.Poisson_fe_face_values[VectorField].value(k,q);
						}

						// loop over all the test function dofs of this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// - \int_{face} p * n * (phi_{Dichlet}) dx
							data.local_rhs(i)	+=	
													 -(scratch.psi_i_field[i] *
														 scratch.normals[q] *
														(scratch.Poisson_bi_values[q] 
														 -
														scratch.Poisson_bc_values[q]) *
														scratch.Poisson_fe_face_values.JxW(q));
						} // for i
					} // for q
				} // end if schottky
			} // if boundary
		} // end for face_no
	} // assemble_local_Poisson_rhs_for_semiconductor

	template <int dim>
	void 
	SolarCellProblem<dim>::
	assemble_local_Poisson_rhs_for_electrolyte(
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

		// Get the actual values for vector field and potential from FEValues
		// Use Extractors instead of having to deal with shapefunctions directly
		const FEValuesExtractors::Vector VectorField(0); // Vector as in Vector field
		const FEValuesExtractors::Scalar Potential(dim);
		const FEValuesExtractors::Scalar Density(dim);

		// reset the local_rhs vector to be zero
		data.local_rhs=0;

		// get the Poisson cell info
		std::pair<unsigned int, unsigned int> Poisson_cell_info = 
					e_2_p_map[std::pair<unsigned int, unsigned int>(cell->level(),
																													cell->index())];

		// get the Poisson cell
		typename DoFHandler<dim>::active_cell_iterator
								Poisson_cell(&Poisson_triangulation,
													Poisson_cell_info.first,
													Poisson_cell_info.second,
													&Poisson_dof_handler);


		Poisson_cell->get_dof_indices(data.local_dof_indices);
		scratch.Poisson_fe_values.reinit(Poisson_cell);
		// get the test rhs for poisson
		// Assemble the right hand side on semiconductor side
		scratch.carrier_fe_values.reinit(cell);

		// get the reductant density values at the previous time step
		scratch.carrier_fe_values[Density].get_function_values(
																				redox_pair.carrier_1.solution,
																				scratch.old_carrier_1_density_values);

		// get the oxidant density values at the previous time step
		scratch.carrier_fe_values[Density].get_function_values(
																				redox_pair.carrier_2.solution,
																				scratch.old_carrier_2_density_values);
	
		// get background charge profile values on this cell
//		reductants_e.value_list(scratch.carrier_fe_values.get_quadrature_points(),
//													scratch.donor_doping_values,
//													dim); // calls the density values of the donor profile

//		oxidants_e.value_list(scratch.carrier_fe_values.get_quadrature_points(),
//											scratch.acceptor_doping_values,
//											dim); // calls the density values of the donor profile
	

		// Loop over all the quadrature points in this cell
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// copy over the test functions
			for(unsigned int k = 0; k < dofs_per_cell; k++)
				scratch.psi_i_potential[k] = scratch.Poisson_fe_values[Potential].value(k,q);

			// loop over the test function dofs for this cell
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				// i-th potential basis functions at the point q
		
				// get the local RHS values for this cell
				// = - 1/\lambda^{2} \int v (-\rho_{r} + \rho_{o} )dx
				data.local_rhs(i) += -scratch.psi_i_potential[i] * 
																(
//															 (scratch.donor_doping_values[q] 
//													 			- 
//														 		scratch.acceptor_doping_values[q])
//												    	 	+
															(redox_pair.carrier_1.charge_number *
																scratch.old_carrier_1_density_values[q]	
																+
												    	 	redox_pair.carrier_2.charge_number *
																scratch.old_carrier_2_density_values[q])
														  ) * scratch.Poisson_fe_values.JxW(q);
			} // for i
		} // for q
	// No need to apply Dirichlet boundary conditions as the potential is zero
		
		// loop over all the faces of this cell to calculate the vector
		// from the dirichlet boundary conditions if the face is on the 
		// Dirichlet portion of the boundary
		for(unsigned int face_no=0;
				face_no<GeometryInfo<dim>::faces_per_cell;
				face_no++)
		{
			// obtain the face_no-th face of this cell
			typename DoFHandler<dim>::face_iterator 	face = Poisson_cell->face(face_no);

			// apply Dirichlet boundary conditions.. 
			// Since we are in a semicondcutor cell we know to apply the 
			// biases
			if((face->at_boundary()) && (face->boundary_id() == Dirichlet))
			{	
				// get the values of the shape functions at this boundary face
				scratch.Poisson_fe_face_values.reinit(Poisson_cell,face_no);
			
				// get the values of the dirichlet boundary conditions evaluated
				// on the quadrature points of this face
				bulk_bias.value_list(
														scratch.Poisson_fe_face_values.get_quadrature_points(),
														scratch.Poisson_bc_values);

				// copy over normal vectors
				for(unsigned int k=0; k < n_face_q_points; k++)
					scratch.normals[k] = scratch.Poisson_fe_face_values.normal_vector(k);

				// loop over all the quadrature points of this face
				for(unsigned int q=0; q<n_face_q_points; q++)
				{
					// copy over the test functions
					for(unsigned int k = 0; k < dofs_per_cell; k++)
						scratch.psi_i_field[k] = scratch.Poisson_fe_face_values[VectorField].value(k,q);
					
					// loop over all the test function dofs of this face
					for(unsigned int i=0; i<dofs_per_cell; i++)
					{
						// - \int_{face} p * n * (phi_{Dichlet}) dx
														//scratch.Poisson_fe_face_values[VectorField].value(i,q) *
						data.local_rhs(i)	+=	
													  -scratch.psi_i_field[i] * 
														scratch.normals[q] *
														scratch.Poisson_bc_values[q] *
														scratch.Poisson_fe_face_values.JxW(q);
						} // for i
					} // for q
				} // end if
			} // end for face_no
	} // end assemble_local_Poisson_rhs_for_electrolyte


	/////////////////////////////////////////////////////////////////////////////////
	// 								LOCAL DISCONTINUOUS GALERKIN ROUTINES
	/////////////////////////////////////////////////////////////////////////////////

	template <int dim>
	void 
	SolarCellProblem<dim>::
	assemble_LDG_system(const double & transient_or_steady)
	{	
		// assemble mass matrix for the electron_hole pair
		WorkStream::run(semiconductor_dof_handler.begin_active(),
										semiconductor_dof_handler.end(),
										std_cxx11::bind(&LDG_System::LDG<dim>::
																		assemble_local_LDG_mass_matrix,
																		LDG_Assembler, // the Assembler object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		delta_t),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_semiconductor_mass_matrix,
																		this, 			// tthis object
																		std_cxx11::_1),
										Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
										Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
										);

		// system matrix for the electron_hole pair
		WorkStream::run(semiconductor_dof_handler.begin_active(),
										semiconductor_dof_handler.end(),
										std_cxx11::bind(&LDG_System::LDG<dim>::
																		assemble_local_LDG_cell_and_bc_terms,
																		LDG_Assembler, // the Assembler object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		electron_hole_pair.carrier_1.scaled_mobility,
																		electron_hole_pair.carrier_2.scaled_mobility,
																		delta_t,
																		transient_or_steady,
																		electron_hole_pair.penalty),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_semiconductor_system_matrix,
																		this, 			// tthis object
																		std_cxx11::_1),
										Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
										Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
										);
		// LDG FLLUXES
		LDG_Assembler.assemble_flux_terms(semiconductor_dof_handler,
												electron_hole_pair,
												Poisson_fe,
												carrier_fe);	

		if(full_system)
		{
			// assemble the mass matrix for redox pair
			WorkStream::run(electrolyte_dof_handler.begin_active(),
											electrolyte_dof_handler.end(),
											std_cxx11::bind(&LDG_System::LDG<dim>::
																		assemble_local_LDG_mass_matrix,
																		LDG_Assembler, // the Assembler object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		delta_t),
											std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_electrolyte_mass_matrix,
																		this, 			// tthis object
																		std_cxx11::_1),
											Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
											Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
											);
	
			// assemble the mass matrix for redox pair
			WorkStream::run(electrolyte_dof_handler.begin_active(),
											electrolyte_dof_handler.end(),
											std_cxx11::bind(&LDG_System::LDG<dim>::
																		assemble_local_LDG_cell_and_bc_terms,
																		LDG_Assembler, // the Assembler object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		redox_pair.carrier_1.scaled_mobility,
																		redox_pair.carrier_2.scaled_mobility,
																		delta_t,
																		transient_or_steady,
																		redox_pair.penalty),
											std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_electrolyte_system_matrix,
																		this, 			// tthis object
																		std_cxx11::_1),
											Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
											Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
											);

			// LDG FLUXES
			LDG_Assembler.assemble_flux_terms(electrolyte_dof_handler,
																				redox_pair,
																				Poisson_fe,
																				carrier_fe);
		} // end if
	}

	template<int dim>
	void
	SolarCellProblem<dim>::
	copy_local_to_global_semiconductor_mass_matrix(
				const Assembly::DriftDiffusion::CopyData<dim>			& data)
	{
		// copy local mass matrix into the global mass matrix.. 
		electron_hole_pair.constraints.distribute_local_to_global(
																					data.local_mass_matrix,
																					data.local_dof_indices,
																					electron_hole_pair.mass_matrix);
	}

	template<int dim>
	void
	SolarCellProblem<dim>::
	copy_local_to_global_electrolyte_mass_matrix(
				const Assembly::DriftDiffusion::CopyData<dim>			& data)
	{
		// copy local mass matrix into the global mass matrix.. 
		redox_pair.constraints.distribute_local_to_global(
																					data.local_mass_matrix,
																					data.local_dof_indices,
																					redox_pair.mass_matrix);
	}


	template<int dim>
	void
	SolarCellProblem<dim>::
	copy_local_to_global_semiconductor_system_matrix(
				const Assembly::DriftDiffusion::CopyData<dim>			& data)
	{
		// copy local system matrix into the global system matrix.. 
		electron_hole_pair.constraints.distribute_local_to_global(
																			data.local_matrix_1,
																			data.local_dof_indices,
																			electron_hole_pair.carrier_1.system_matrix);

		// copy local system matrix into the global system matrix.. 
		electron_hole_pair.constraints.distribute_local_to_global(
																			data.local_matrix_2,
																			data.local_dof_indices,
																			electron_hole_pair.carrier_2.system_matrix);
	}

	template<int dim>
	void
	SolarCellProblem<dim>::
	copy_local_to_global_electrolyte_system_matrix(
				const Assembly::DriftDiffusion::CopyData<dim>			& data)
	{
		// copy local system matrix into the global system matrix.. 
		redox_pair.constraints.distribute_local_to_global(
																		data.local_matrix_1,
																		data.local_dof_indices,
																		redox_pair.carrier_1.system_matrix);

		// copy local syste, matrix into the global system matrix.. 
		redox_pair.constraints.distribute_local_to_global(
																		data.local_matrix_2,
																		data.local_dof_indices,
																		redox_pair.carrier_2.system_matrix);
	}

	template<int dim>
	void
	SolarCellProblem<dim>::
	copy_local_to_global_semiconductor_system_rhs(						
						const Assembly::DriftDiffusion::CopyData<dim>			& data)
	{
		// copy local system rhs into the global system rhs.. 
		electron_hole_pair.constraints.distribute_local_to_global(
																		data.local_carrier_1_rhs,
																		data.local_dof_indices,
																		electron_hole_pair.carrier_1.system_rhs);

		// copy local syste, matrix into the global system matrix.. 
		electron_hole_pair.constraints.distribute_local_to_global(
																		data.local_carrier_2_rhs,
																		data.local_dof_indices,
																		electron_hole_pair.carrier_2.system_rhs);
	}

	template<int dim>
	void
	SolarCellProblem<dim>::
	copy_local_to_global_electrolyte_system_rhs(						
				const Assembly::DriftDiffusion::CopyData<dim>			& data)
	{
		// copy local system rhs into the global system rhs.. 
		redox_pair.constraints.distribute_local_to_global(
																		data.local_carrier_1_rhs,
																		data.local_dof_indices,
																		redox_pair.carrier_1.system_rhs);

		// copy local syste, matrix into the global system matrix.. 
		redox_pair.constraints.distribute_local_to_global(
																		data.local_carrier_2_rhs,
																		data.local_dof_indices,
																		redox_pair.carrier_2.system_rhs);
	}

	template<int dim>
	void
	SolarCellProblem<dim>::
	assemble_semiconductor_rhs()
	{
		// set carrier_system_rhs = M * u^{n-1}
		electron_hole_pair.mass_matrix.vmult(electron_hole_pair.carrier_1.system_rhs,
																				 electron_hole_pair.carrier_1.solution);

		electron_hole_pair.mass_matrix.vmult(electron_hole_pair.carrier_2.system_rhs,
																				 electron_hole_pair.carrier_2.solution);

		// now run through the semiconductor cells and assemble rhs for the 
		// LDG equations.
		WorkStream::run(semiconductor_dof_handler.begin_active(),
										semiconductor_dof_handler.end(),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		assemble_local_semiconductor_rhs,
																		this, // this object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		electron_hole_pair.penalty),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_semiconductor_system_rhs,
																		this, // this object
																		std_cxx11::_1),
										Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
										Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
										);

	}		

	template<int dim>
	void 
	SolarCellProblem<dim>::
	assemble_local_semiconductor_rhs(
								const typename DoFHandler<dim>::active_cell_iterator & cell,
								Assembly::AssemblyScratch<dim>											 & scratch,
								Assembly::DriftDiffusion::CopyData<dim>							 & data,
								const double 																				 & penalty)									
	{
		/*-------------------------------------------------------------------*/
		// This is for the case of a coupled Poisson's equation and
		// reactive interface
		/*-------------------------------------------------------------------*/

		// this assembles the drift term in the ldg formulation.  it uses the 
		// electric field at the current iteration and the density of the 
		// carrier at the previous time step
		const unsigned int dofs_per_cell			 = 
															scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points					 = 
																scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points		 = 
																scratch.carrier_fe_face_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		double h = cell->diameter();

		// get the Poisson cell info
		std::pair<unsigned int, unsigned int> Poisson_cell_info = 
					s_2_p_map[std::pair<unsigned int, unsigned int>(cell->level(),
																													cell->index())];

		// get the Poisson cell
		typename DoFHandler<dim>::active_cell_iterator
								Poisson_cell(&Poisson_triangulation,
													Poisson_cell_info.first,
													Poisson_cell_info.second,
													&Poisson_dof_handler);

		// reinitialize the fe_values 
		scratch.carrier_fe_values.reinit(cell);
		scratch.Poisson_fe_values.reinit(Poisson_cell);
		
		// reset the local_rhs to be zero
		data.local_carrier_1_rhs=0;
		data.local_carrier_2_rhs=0;
	
		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);

		const FEValuesExtractors::Vector ElectricField(0);

		// get the values of carrier_1 and carrier_2 densities at the pevious time step
		scratch.carrier_fe_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.old_carrier_1_density_values);

		scratch.carrier_fe_values[Density].get_function_values(
																				electron_hole_pair.carrier_2.solution,
																				scratch.old_carrier_2_density_values);

		generation.value_list(scratch.carrier_fe_values.get_quadrature_points(),
													scratch.generation_values);

		// get the electric field values at the previous time step
		scratch.Poisson_fe_values[ElectricField].get_function_values(
																									Poisson_object.solution,
																									scratch.electric_field_values);

		const double inverse_perm			=		1.0/sim_params.semiconductor_permittivity;

		// loop over all the quadrature points in this cell and compute body integrals
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// copy over the test functions
			for(unsigned int k=0; k<dofs_per_cell; k++)
				scratch.psi_i_density[k] = scratch.carrier_fe_values[Density].value(k,q);

			for(unsigned int k=0; k<dofs_per_cell; k++)
				scratch.psi_i_current[k] = scratch.carrier_fe_values[Current].value(k,q);

			// loop over all the test function dofs and get the test functions
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				// contribution from RHS function + Drift
				// int_{Omega} v * R dx
				data.local_carrier_1_rhs(i) += ( 
													(scratch.psi_i_density[i] * scratch.generation_values[q])
													+
													(scratch.psi_i_density[i] *
													SRH_Recombination(
															scratch.old_carrier_1_density_values[q],
															scratch.old_carrier_2_density_values[q],
															sim_params))
													+
														electron_hole_pair.carrier_1.charge_number *
														scratch.psi_i_current[i] *
														inverse_perm *
														scratch.electric_field_values[q] *
														scratch.old_carrier_1_density_values[q]
														) *
														scratch.carrier_fe_values.JxW(q);

				data.local_carrier_2_rhs(i) += ( 
													(scratch.psi_i_density[i] * scratch.generation_values[q])
													+
													(scratch.psi_i_density[i] *
													SRH_Recombination(
															scratch.old_carrier_1_density_values[q],
															scratch.old_carrier_2_density_values[q],
															sim_params))
														+
														electron_hole_pair.carrier_2.charge_number *
														scratch.psi_i_current[i] *
														inverse_perm *	
														scratch.electric_field_values[q] *
														scratch.old_carrier_2_density_values[q]
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
					// Get the doping profile values for the boundary conditions
					electrons_e.value_list(
														scratch.carrier_fe_face_values.get_quadrature_points(),
														scratch.carrier_1_bc_values,
														dim); // calls the density values of the donor profile
																	// not the current ones
					holes_e.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_2_bc_values,
													dim); // calls the density values of the donor profile
																// not the current ones
					// copy over normal vectors
					for(unsigned int k=0; k<n_face_q_points; k++)
						scratch.normals[k] = scratch.carrier_fe_face_values.normal_vector(k);
	
					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k=0; k<dofs_per_cell; k++)
							scratch.psi_i_density[k] = scratch.carrier_fe_face_values[Density].value(k,q);

						for(unsigned int k=0; k<dofs_per_cell; k++)
							scratch.psi_i_current[k] = scratch.carrier_fe_face_values[Current].value(k,q);

						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// int_{\Gamma_{D}} ( -p^{-} n^{-} + penalty/h * v^{-}) * u_{D} ds
							data.local_carrier_1_rhs(i) += 
																(-1.0 * scratch.psi_i_current[i] *
													 			 scratch.normals[q]
																 + 
																 (penalty/h) * scratch.psi_i_density[i]) *
																 scratch.carrier_1_bc_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q);				
						
							data.local_carrier_2_rhs(i) +=  
																(-1.0 * scratch.psi_i_current[i] *
																scratch.normals[q] 
																 + 
																 (penalty/h) * scratch.psi_i_density[i]) *
																 scratch.carrier_2_bc_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q);				
						} // for i
					}	// for q
				} // end Dirichlet
				else if(face->boundary_id() == Interface)
				{
					// get the doping profile values for the boundary conditions
					electrons_e.value_list(
														scratch.carrier_fe_face_values.get_quadrature_points(),
														scratch.carrier_1_bc_values,
														dim); // calls the density values of the donor profile
																	// not the current ones
					holes_e.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_2_bc_values,
													dim); // calls the density values of the donor profile
																// not the current ones
					
	
					// get the electrona and hole densities at the pevious time step
					scratch.carrier_fe_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.electron_interface_values);

					scratch.carrier_fe_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_2.solution,
																				scratch.hole_interface_values);

					// get neighboring selectrolyte cell
					unsigned interface_index = semi_interface_map[
												std::pair<unsigned int, unsigned>(cell->level(),
																													cell->index())];
					typename DoFHandler<dim>::active_cell_iterator
										neighbor_cell(&electrolyte_triangulation,
																	elec_interface_cells[interface_index].first, // level
																	elec_interface_cells[interface_index].second, // index
																	&electrolyte_dof_handler);

					// update fe_neigbors_face_values for this face
					scratch.carrier_fe_neighbor_face_values.reinit(neighbor_cell,
																												 elec_interface_faces[interface_index]);

					// get the reductant and oxidant densities at the pevious time step
					scratch.carrier_fe_neighbor_face_values[Density].get_function_values(
																				redox_pair.carrier_1.solution,
																				scratch.reductant_interface_values);

					scratch.carrier_fe_neighbor_face_values[Density].get_function_values(
																				redox_pair.carrier_2.solution,
																				scratch.oxidant_interface_values);
		
					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k=0; k<dofs_per_cell; k++)
							scratch.psi_i_density[k] = scratch.carrier_fe_face_values[Density].value(k,q);

						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// int_{\sigma} -v^{-} ket (rho_n - rho_n^e) rho_o ds
							data.local_carrier_1_rhs(i) += 
																-1.0 * scratch.psi_i_density[i] *
																 sim_params.scaled_k_et *
																 (scratch.electron_interface_values[q]
																	-
																	scratch.carrier_1_bc_values[q]) * 
																 scratch.oxidant_interface_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q);				
					
							// int_{\Sigma}  v^{-} kht (rho_p - rho_p^e) rho_r ds
							data.local_carrier_2_rhs(i) +=  
																	+1.0 * scratch.psi_i_density[i] *
																 sim_params.scaled_k_ht *
																 (scratch.hole_interface_values[q]
																	-
																 scratch.carrier_2_bc_values[q]) *
																 scratch.reductant_interface_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q);						
						} // for i
					} // for j

				}
				else if(face->boundary_id() == Schottky)
				{
					// Get the doping profile values for the boundary conditions
					electrons_e.value_list(
														scratch.carrier_fe_face_values.get_quadrature_points(),
														scratch.carrier_1_bc_values,
														dim); // calls the density values of the donor profile
																	// not the current ones
					holes_e.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_2_bc_values,
													dim); // calls the density values of the donor profile
																// not the current ones
					
	
					// get the electrona and hole densities at the pevious time step
					scratch.carrier_fe_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.electron_interface_values);

					scratch.carrier_fe_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_2.solution,
																				scratch.hole_interface_values);

	
					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k=0; k<dofs_per_cell; k++)
							scratch.psi_i_density[k] = scratch.carrier_fe_face_values[Density].value(k,q);

						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// int_{\Sigma} -v^{-} ket (rho_n - rho_n^e) rho_o ds
							data.local_carrier_1_rhs(i) += 
																-1.0 * scratch.psi_i_density[i] *
																 sim_params.scaled_electron_recombo_v*
																 (scratch.electron_interface_values[q]
																	-
																 scratch.carrier_1_bc_values[q] ) * 
															 	 scratch.carrier_fe_face_values.JxW(q);				
					
							// int_{\Sigma}  v^{-} kht (rho_p - rho_p^e) rho_r ds
							data.local_carrier_2_rhs(i) +=  
																	+1.0 * scratch.psi_i_density[i] *
																 sim_params.scaled_hole_recombo_v *
																 (scratch.hole_interface_values[q]
																	-
																 scratch.carrier_2_bc_values[q] ) *
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
		} // end for face_no
	
	} // end assemble_local_semiconductor_rhs


	template<int dim>
	void
	SolarCellProblem<dim>::
	assemble_electrolyte_rhs()
	{
		// set carrier_system_rhs = M * u^{n-1}
		redox_pair.mass_matrix.vmult(redox_pair.carrier_1.system_rhs,
																 redox_pair.carrier_1.solution);

		redox_pair.mass_matrix.vmult(redox_pair.carrier_2.system_rhs,
																 redox_pair.carrier_2.solution);

		WorkStream::run(electrolyte_dof_handler.begin_active(),
										electrolyte_dof_handler.end(),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		assemble_local_electrolyte_rhs,
																		this, // this object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		redox_pair.penalty),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_electrolyte_system_rhs,
																		this, // this object
																		std_cxx11::_1),
										Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
										Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
										);

	}

	template<int dim>
	void 
	SolarCellProblem<dim>::assemble_local_electrolyte_rhs(
								const typename DoFHandler<dim>::active_cell_iterator & cell,
								Assembly::AssemblyScratch<dim>											 & scratch,
								Assembly::DriftDiffusion::CopyData<dim>							 & data,
								const double																				 & penalty)									
	{	
		/*-------------------------------------------------------------------*/
		// This is for the case of a coupled Poisson's equation and
		// reactive interface
		/*-------------------------------------------------------------------*/

		// this assembles the drift term in the ldg formulation.  it uses the 
		// electric field at the current iteration and the density of the 
		// carrier at the previous time step
		const unsigned int dofs_per_cell			 = 
															scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points					 = 
																scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points		 = 
																scratch.carrier_fe_face_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		// get the Poisson cell info
		std::pair<unsigned int, unsigned int> Poisson_cell_info = 
					e_2_p_map[std::pair<unsigned int, unsigned int>(cell->level(),
																													cell->index())];

		// get the Poisson cell
		typename DoFHandler<dim>::active_cell_iterator
								Poisson_cell(&Poisson_triangulation,
													Poisson_cell_info.first,
													Poisson_cell_info.second,
													&Poisson_dof_handler);

		// reinitialize the fe_values 
		scratch.carrier_fe_values.reinit(cell);
		scratch.Poisson_fe_values.reinit(Poisson_cell);
	
		double h =  cell->diameter();
	
		// reset the local_rhs to be zero
		data.local_carrier_1_rhs=0;
		data.local_carrier_2_rhs=0;
	
		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);

		const FEValuesExtractors::Vector ElectricField(0);

		// get the values of carrier_1 and carrier_2 densities at the pevious time step
		scratch.carrier_fe_values[Density].get_function_values(
																				redox_pair.carrier_1.solution,
																				scratch.old_carrier_1_density_values);

		scratch.carrier_fe_values[Density].get_function_values(
																				redox_pair.carrier_2.solution,
																				scratch.old_carrier_2_density_values);

		// get the electric field values at the previous time step
		scratch.Poisson_fe_values[ElectricField].get_function_values(
																									Poisson_object.solution,
																									scratch.electric_field_values);

		const double inverse_perm			=		1.0/sim_params.electrolyte_permittivity;

		// loop over all the quadrature points in this cell and compute body integrals
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// copy over the test functions
			for(unsigned int k=0; k<dofs_per_cell; k++)
				scratch.psi_i_density[k] = scratch.carrier_fe_values[Density].value(k,q);

			for(unsigned int k=0; k<dofs_per_cell; k++)
				scratch.psi_i_current[k] = scratch.carrier_fe_values[Current].value(k,q);

			// loop over all the test function dofs and get the test functions
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{

				// contribution from RHS function + Drift
				// int_{Omega} v * R dx
				data.local_carrier_1_rhs(i) += ( 
														redox_pair.carrier_1.charge_number *
														scratch.psi_i_current[i] *
														inverse_perm *	
														scratch.electric_field_values[q] *
														scratch.old_carrier_1_density_values[q]
														) *
														scratch.carrier_fe_values.JxW(q);

				data.local_carrier_2_rhs(i) += ( 
														redox_pair.carrier_2.charge_number *
														scratch.psi_i_current[i] *
														inverse_perm *
														scratch.electric_field_values[q] *
														scratch.old_carrier_2_density_values[q]
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
					// Get the doping profile values for the boundary conditions
					reductants_e.value_list(
														scratch.carrier_fe_face_values.get_quadrature_points(),
														scratch.carrier_1_bc_values,
														dim); // calls the density values of the donor profile
																	// not the current ones
					oxidants_e.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_2_bc_values,
													dim); // calls the density values of the donor profile
																// not the current ones
						// copy over normal vectors
					for(unsigned int k=0; k<n_face_q_points; k++)
						scratch.normals[k] = scratch.carrier_fe_face_values.normal_vector(k);
		
					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k=0; k<dofs_per_cell; k++)
							scratch.psi_i_density[k] = scratch.carrier_fe_face_values[Density].value(k,q);

						for(unsigned int k=0; k<dofs_per_cell; k++)
							scratch.psi_i_current[k] = scratch.carrier_fe_face_values[Current].value(k,q);

						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{

							// int_{\Gamma_{D}} (-p^{-} n^{-} + penalty * v^{-}) *u_{D} ds
							data.local_carrier_1_rhs(i) += 
																(-1.0 * scratch.psi_i_current[i] *
																scratch.normals[q]
																 +
																 (penalty/h) * scratch.psi_i_density[i]) *
																 scratch.carrier_1_bc_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q);				
					
							data.local_carrier_2_rhs(i) +=  
																(-1.0 * scratch.psi_i_current[i] *
																scratch.normals[q]
	  														 +
																 (penalty/h) * scratch.psi_i_density[i]) *
																 scratch.carrier_2_bc_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q);				


						} // for i
					}	// for q
				} // end Dirichlet
				else if(face->boundary_id() == Interface)
				{

					// get the reductant and oxidant densities at the pevious time step
					scratch.carrier_fe_face_values[Density].get_function_values(
																				redox_pair.carrier_1.solution,
																				scratch.reductant_interface_values);

					scratch.carrier_fe_face_values[Density].get_function_values(
																				redox_pair.carrier_2.solution,
																				scratch.oxidant_interface_values);

				  // get neighboring semiconductor cell
					unsigned interface_index = elec_interface_map[
												std::pair<unsigned int, unsigned>(cell->level(),
																													cell->index())];
					typename DoFHandler<dim>::active_cell_iterator
										neighbor_cell(&semiconductor_triangulation,
																	semi_interface_cells[interface_index].first, // level
																	semi_interface_cells[interface_index].second, // index
																	&semiconductor_dof_handler);

					// update fe_neigbors_face_values for this face
					scratch.carrier_fe_neighbor_face_values.reinit(neighbor_cell,
																												 semi_interface_faces[interface_index]);

					// Get the doping profile values for the interface conditions
					electrons_e.value_list(
														scratch.carrier_fe_neighbor_face_values.get_quadrature_points(),
														scratch.carrier_1_bc_values,
														dim); // calls the density values of the donor profile
																	// not the current ones
					holes_e.value_list(
													scratch.carrier_fe_neighbor_face_values.get_quadrature_points(),
													scratch.carrier_2_bc_values,
													dim); // calls the density values of the donor profile
																// not the current ones
					
	
					// get the electrona and hole densities at the pevious time step
					scratch.carrier_fe_neighbor_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.electron_interface_values);

					scratch.carrier_fe_neighbor_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_2.solution,
																				scratch.hole_interface_values);

					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k=0; k<dofs_per_cell; k++)
							scratch.psi_i_density[k] = scratch.carrier_fe_face_values[Density].value(k,q);
		
						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// get the test function
							const double current = -1.0 * scratch.psi_i_density[i] *
																 sim_params.scaled_k_et *
																 (scratch.electron_interface_values[q]
																	-
																	scratch.carrier_1_bc_values[q]) * 
																 scratch.oxidant_interface_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q) 
																 +1.0 * scratch.psi_i_density[i] *
																 sim_params.scaled_k_ht *
																 (scratch.hole_interface_values[q]
																	-
																	scratch.carrier_2_bc_values[q]) *
																 scratch.reductant_interface_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q);				
		
	
							// int_{\Sigma} -v^{-} ket (rho_n - rho_n^e) rho_o ds 
							// 	+
							// int_{\Sigma}  v^{-} kht (rho_p - rho_p^e) rho_r ds
							data.local_carrier_1_rhs(i) += current;
							data.local_carrier_2_rhs(i) += -1.0 *current;
		
						} // for i
					} // for j
				}
				else if(face->boundary_id() == Neumann)
				{
					// NOTHIN TO DO IF INSULATING
				}
				else
					Assert(false, ExcNotImplemented() );
			} // end at boundary
		} // end for face_no
	
	} ///

	template<int dim>
	void
	SolarCellProblem<dim>::
	assemble_Schottky_rhs()
	{
		// set carrier_system_rhs = M * u^{n-1}
		electron_hole_pair.mass_matrix.vmult(electron_hole_pair.carrier_1.system_rhs,
																				 electron_hole_pair.carrier_1.solution);

		electron_hole_pair.mass_matrix.vmult(electron_hole_pair.carrier_2.system_rhs,
																				 electron_hole_pair.carrier_2.solution);

		// now run through the semiconductor cells and assemble rhs for the 
		// LDG equations.
		WorkStream::run(semiconductor_dof_handler.begin_active(),
										semiconductor_dof_handler.end(),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		assemble_local_Schottky_rhs,
																		this, // this object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		electron_hole_pair.penalty),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_semiconductor_system_rhs,
																		this, // this object
																		std_cxx11::_1),
										Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
										Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
										);

	}		// end assmble Schottky

	template<int dim>
	void 
	SolarCellProblem<dim>::
	assemble_local_Schottky_rhs(
								const typename DoFHandler<dim>::active_cell_iterator & cell,
								Assembly::AssemblyScratch<dim>											 & scratch,
								Assembly::DriftDiffusion::CopyData<dim>							 & data,
								const double 																				 & penalty)									
	{
		/*-------------------------------------------------------------------*/
		// This is for the case of a coupled Poisson's equation and
		// reactive interface
		/*-------------------------------------------------------------------*/

		// this assembles the drift term in the ldg formulation.  it uses the 
		// electric field at the current iteration and the density of the 
		// carrier at the previous time step
		const unsigned int dofs_per_cell			 = 
															scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points					 = 
																scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points		 = 
																scratch.carrier_fe_face_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		double h = cell->diameter();

		// get the Poisson cell info
		std::pair<unsigned int, unsigned int> Poisson_cell_info = 
					s_2_p_map[std::pair<unsigned int, unsigned int>(cell->level(),
																													cell->index())];

		// get the Poisson cell
		typename DoFHandler<dim>::active_cell_iterator
								Poisson_cell(&Poisson_triangulation,
													Poisson_cell_info.first,
													Poisson_cell_info.second,
													&Poisson_dof_handler);

		// reinitialize the fe_values 
		scratch.carrier_fe_values.reinit(cell);
		scratch.Poisson_fe_values.reinit(Poisson_cell);
		
		// reset the local_rhs to be zero
		data.local_carrier_1_rhs=0;
		data.local_carrier_2_rhs=0;
	
		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);

		const FEValuesExtractors::Vector ElectricField(0);

		// get the values of carrier_1 and carrier_2 densities at the pevious time step
		scratch.carrier_fe_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.old_carrier_1_density_values);

		scratch.carrier_fe_values[Density].get_function_values(
																				electron_hole_pair.carrier_2.solution,
																				scratch.old_carrier_2_density_values);

		generation.value_list(scratch.carrier_fe_values.get_quadrature_points(),
													scratch.generation_values);

		// get the electric field values at the previous time step
		scratch.Poisson_fe_values[ElectricField].get_function_values(
																									Poisson_object.solution,
																									scratch.electric_field_values);

		const double inverse_perm			=		1.0/sim_params.semiconductor_permittivity;

		// loop over all the quadrature points in this cell and compute body integrals
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// copy over the test functions
			for(unsigned int k=0; k<dofs_per_cell; k++)
				scratch.psi_i_density[k] = scratch.carrier_fe_values[Density].value(k,q);

			for(unsigned int k=0; k<dofs_per_cell; k++)
				scratch.psi_i_current[k] = scratch.carrier_fe_values[Current].value(k,q);

			// loop over all the test function dofs and get the test functions
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				// contribution from RHS function + Drift
				// int_{Omega} v * R dx
				data.local_carrier_1_rhs(i) += ( 
													(scratch.psi_i_density[i] * scratch.generation_values[q])
													+
													(scratch.psi_i_density[i] *
													SRH_Recombination(
															scratch.old_carrier_1_density_values[q],
															scratch.old_carrier_2_density_values[q],
															sim_params))
													+
														electron_hole_pair.carrier_1.charge_number *
														scratch.psi_i_current[i] *
														inverse_perm *
														scratch.electric_field_values[q] *
														scratch.old_carrier_1_density_values[q]
														) *
														scratch.carrier_fe_values.JxW(q);

				data.local_carrier_2_rhs(i) += ( 
													(scratch.psi_i_density[i] * scratch.generation_values[q])
													+
													(scratch.psi_i_density[i] *
													SRH_Recombination(
															scratch.old_carrier_1_density_values[q],
															scratch.old_carrier_2_density_values[q],
															sim_params))
														+
														electron_hole_pair.carrier_2.charge_number *
														scratch.psi_i_current[i] *
														inverse_perm *	
														scratch.electric_field_values[q] *
														scratch.old_carrier_2_density_values[q]
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
					// Get the doping profile values for the boundary conditions
					electrons_e.value_list(
														scratch.carrier_fe_face_values.get_quadrature_points(),
														scratch.carrier_1_bc_values,
														dim); // calls the density values of the donor profile
																	// not the current ones
					holes_e.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_2_bc_values,
													dim); // calls the density values of the donor profile
																// not the current ones
					// copy over normal vectors
					for(unsigned int k=0; k<n_face_q_points; k++)
						scratch.normals[k] = scratch.carrier_fe_face_values.normal_vector(k);
	
					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k=0; k<dofs_per_cell; k++)
							scratch.psi_i_density[k] = scratch.carrier_fe_face_values[Density].value(k,q);

						for(unsigned int k=0; k<dofs_per_cell; k++)
							scratch.psi_i_current[k] = scratch.carrier_fe_face_values[Current].value(k,q);

						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// int_{\Gamma_{D}} ( -p^{-} n^{-} + penalty/h * v^{-}) * u_{D} ds
							data.local_carrier_1_rhs(i) += 
																(-1.0 * scratch.psi_i_current[i] *
													 			 scratch.normals[q]
																 + 
																 (penalty/h) * scratch.psi_i_density[i]) *
																 scratch.carrier_1_bc_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q);				
						
							data.local_carrier_2_rhs(i) +=  
																(-1.0 * scratch.psi_i_current[i] *
																scratch.normals[q] 
																 + 
																 (penalty/h) * scratch.psi_i_density[i]) *
																 scratch.carrier_2_bc_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q);				
						} // for q
					} // for i
				} // if Dirichlet
				else if(face->boundary_id() == Interface)
				{
					// Get the doping profile values for the boundary conditions
					electrons_e.value_list(
														scratch.carrier_fe_face_values.get_quadrature_points(),
														scratch.carrier_1_bc_values,
														dim); // calls the density values of the donor profile
																	// not the current ones
					holes_e.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_2_bc_values,
													dim); // calls the density values of the donor profile
																// not the current ones
					
	
					// get the electrona and hole densities at the pevious time step
					scratch.carrier_fe_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.electron_interface_values);

					scratch.carrier_fe_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_2.solution,
																				scratch.hole_interface_values);

	
					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k=0; k<dofs_per_cell; k++)
							scratch.psi_i_density[k] = scratch.carrier_fe_face_values[Density].value(k,q);

						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// int_{\Sigma} -v^{-} ket (rho_n - rho_n^e) rho_o ds
							data.local_carrier_1_rhs(i) += 
																-1.0 * scratch.psi_i_density[i] *
																 sim_params.scaled_electron_recombo_v*
																 (scratch.electron_interface_values[q]
																	-
																 scratch.carrier_1_bc_values[q] ) * 
															 	 scratch.carrier_fe_face_values.JxW(q);				
					
							// int_{\Sigma}  v^{-} kht (rho_p - rho_p^e) rho_r ds
							data.local_carrier_2_rhs(i) +=  
																	+1.0 * scratch.psi_i_density[i] *
																 sim_params.scaled_hole_recombo_v *
																 (scratch.hole_interface_values[q]
																	-
																 scratch.carrier_2_bc_values[q] ) *
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
		} // end for face_no
	
	} // end assemble_local_semiconductor_rhs


		////////////////////////////////////////////////////////////////////////
	//														SOLVERS
	////////////////////////////////////////////////////////////////////////

	template <int dim>
	void 
	SolarCellProblem<dim>::
	set_solvers()
	{
		Poisson_object.set_solver();
		electron_hole_pair.carrier_1.set_solver();
		electron_hole_pair.carrier_2.set_solver();
		
		if(full_system)
		{	
			redox_pair.carrier_1.set_solver();
			redox_pair.carrier_2.set_solver();
		}
	}

	
	template <int dim>
	void 
	SolarCellProblem<dim>::
	solve_Poisson()
	{
		Poisson_object.solve();	
	}		
	
	template<int dim>
	void
	SolarCellProblem<dim>::
	solve_full_system()
	{
		Threads::TaskGroup<void> task_group;
		
		task_group += Threads::new_task(&ChargeCarrierSpace::
																		Carrier<dim>::solve,
																		electron_hole_pair.carrier_1);

		task_group += Threads::new_task(&ChargeCarrierSpace::
																		Carrier<dim>::solve,
																		electron_hole_pair.carrier_2);

		task_group += Threads::new_task(&ChargeCarrierSpace::
																		Carrier<dim>::solve,
																		redox_pair.carrier_1);

		task_group += Threads::new_task(&ChargeCarrierSpace::
																		Carrier<dim>::solve,
																		redox_pair.carrier_2);

		task_group.join_all();
	} // solve_full_system

	template<int dim>
	void
	SolarCellProblem<dim>::
	solve_semiconductor_system()
	{
		Threads::TaskGroup<void> task_group;
		
		task_group += Threads::new_task(&ChargeCarrierSpace::
																		Carrier<dim>::solve,
																		electron_hole_pair.carrier_1);

		task_group += Threads::new_task(&ChargeCarrierSpace::
																		Carrier<dim>::solve,
																		electron_hole_pair.carrier_2);
		task_group.join_all();
	} // solve_semiconductor_system

	template<int dim>
	void
	SolarCellProblem<dim>::
	solve_electrolyte_system()
	{
		Threads::TaskGroup<void> task_group;

		task_group += Threads::new_task(&ChargeCarrierSpace::
																		Carrier<dim>::solve,
																		redox_pair.carrier_1);

		task_group += Threads::new_task(&ChargeCarrierSpace::
																		Carrier<dim>::solve,
																		redox_pair.carrier_2);

		task_group.join_all();
	} // solve_electrolyte_system

	////////////////////////////////////////////////////////////////////////
	//													PRINTING
	////////////////////////////////////////////////////////////////////////

	template<int dim>
	void
	SolarCellProblem<dim>::
	print_results(unsigned int time_step_number)
	{
		// TODO: Make so that rescaled version will work in debug mode
		Threads::TaskGroup<void> task_group;
		task_group += Threads::new_task(&MixedPoisson::
																		MixedFEM<dim>::output_rescaled_results,
																		Mixed_Assembler,
																		Poisson_dof_handler,
																		Poisson_object.solution,
																		sim_params,
																		time_step_number);
			

		task_group += Threads::new_task(&LDG_System::
																		LDG<dim>::output_rescaled_results,
																		LDG_Assembler,
																		semiconductor_dof_handler,
																		electron_hole_pair,
																		sim_params,
																		time_step_number);
		if(full_system)
		{
				task_group += Threads::new_task(&LDG_System::
																			LDG<dim>::output_rescaled_results,
																			LDG_Assembler,
																			electrolyte_dof_handler,
																			redox_pair,
																			sim_params,
																 			time_step_number);
		}

		task_group.join_all();
	} // print_results


	template<int dim>
	void
	SolarCellProblem<dim>::
	print_results_on_boundary(unsigned int time_step_number)
	{
		// TODO: Make so that rescaled version will work in debug mode
		Threads::TaskGroup<void> task_group;
		task_group += Threads::new_task(&LDG_System::
																		LDG<dim>::output_unscaled_results_on_boundary,
																		LDG_Assembler,
																		semiconductor_dof_handler,
																		electron_hole_pair,
		//																sim_params,
																		time_step_number);
		if(full_system)
		{
				task_group += Threads::new_task(&LDG_System::
																			LDG<dim>::output_unscaled_results_on_boundary,
																			LDG_Assembler,
																			electrolyte_dof_handler,
																			redox_pair,
//																			sim_params,
																 			time_step_number);
		}

		task_group.join_all();
		}


	/*--------------------------------------------------------------------*/

		template<int dim>
		void 
		SolarCellProblem<dim>::
		print_convergence(double & time)
		{
			std::cout << "time = " << time << std::endl;
			diff_electrons.add(-1,electron_hole_pair.carrier_1.solution);
			diff_holes.add(-1,electron_hole_pair.carrier_2.solution);
			diff_reductants.add(-1, redox_pair.carrier_1.solution);
			diff_oxidants.add(-1, redox_pair.carrier_2.solution);
			diff_Poisson.add(-1,Poisson_object.solution);

			const ComponentSelectFunction<dim> potential_mask(dim, dim+1);

  		Vector<double> cellwise_errors(semiconductor_triangulation.n_active_cells() );

			QTrapez<1>				q_trapez;
  		QIterated<dim> 		quadrature(q_trapez, degree+1);


 			 VectorTools::integrate_difference(semiconductor_dof_handler, 
																				diff_electrons, 
																				ZeroFunction<dim>(),
                                  			cellwise_errors, 
																				quadrature,
                                    		VectorTools::L2_norm,
                                     		&potential_mask);


				std::cout << "electrons diff: " 
									<< cellwise_errors.l2_norm()
									<< std::endl;

	 			 VectorTools::integrate_difference(semiconductor_dof_handler, 
																				diff_holes, 
																				ZeroFunction<dim>(),
                                  			cellwise_errors, 
																				quadrature,
                                    		VectorTools::L2_norm,
                                     		&potential_mask);

				std::cout << "holes diff: " 
									<< cellwise_errors.l2_norm() 
									<< std::endl;

			if(full_system)
			{
	 			 VectorTools::integrate_difference(electrolyte_dof_handler, 
																				diff_reductants, 
																				ZeroFunction<dim>(),
                                  			cellwise_errors, 
																				quadrature,
                                    		VectorTools::L2_norm,
                                     		&potential_mask);

				std::cout << "reductant diff: " 
									<< cellwise_errors.l2_norm() 
									<< std::endl;

 			 VectorTools::integrate_difference(electrolyte_dof_handler, 
																				diff_oxidants, 
																				ZeroFunction<dim>(),
                                  			cellwise_errors, 
																				quadrature,
                                    		VectorTools::L2_norm,
                                     		&potential_mask);


				std::cout << "oxidants diff: " 
									<< cellwise_errors.l2_norm() 
									<< std::endl;

				diff_reductants = redox_pair.carrier_1.solution;
				diff_oxidants = redox_pair.carrier_2.solution;
			}

			diff_electrons = electron_hole_pair.carrier_1.solution;
			diff_holes = electron_hole_pair.carrier_2.solution;
			diff_Poisson = Poisson_object.solution;
	}


	/*---------------------------------------------------------------------*/
	template<int dim>
	void 
	SolarCellProblem<dim>::
	print_semiconductor_current()
	{	
		double electron_current_sum = 0;
		double hole_current_sum = 0;
	
		Assembly::AssemblyScratch<dim>	scratch(Poisson_fe,
																						carrier_fe,
																						QGauss<dim>(carrier_fe.degree+2),
																						QGauss<dim-1>(carrier_fe.degree+2));

		Assembly::DriftDiffusion::CopyData<dim>	data(carrier_fe);

		typename DoFHandler<dim>::active_cell_iterator
		cell = semiconductor_dof_handler.begin_active(),
		endc = semiconductor_dof_handler.end();

		for(; cell != endc; cell++)
		{
			const unsigned int n_face_q_points = 
														scratch.carrier_fe_face_values.n_quadrature_points;

			const FEValuesExtractors::Vector Current(0);
			const FEValuesExtractors::Scalar Density(dim);

			for(unsigned int face_no=0;
					face_no < GeometryInfo<dim>::faces_per_cell;
					face_no++)
			{
				typename DoFHandler<dim>::face_iterator	face=cell->face(face_no);
				if( face->at_boundary())
				{
					if(face->boundary_id() == Dirichlet)
					{
						scratch.carrier_fe_face_values.reinit(cell, face_no);

						// get the electrona and hole densities at the pevious time step
						scratch.carrier_fe_face_values[Current].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.carrier_1_current_values);

						scratch.carrier_fe_face_values[Current].get_function_values(
																				electron_hole_pair.carrier_2.solution,
																				scratch.carrier_2_current_values);

						for(unsigned int k=0; k<n_face_q_points; k++)
							scratch.normals[k] = scratch.carrier_fe_face_values.normal_vector(k);
		
	
						// loop over all the quadrature points on this face
						for(unsigned int q=0; q<n_face_q_points; q++)
						{
							// int_{\Sigma} 
							electron_current_sum	+= 
																electron_hole_pair.carrier_1.charge_number *
																scratch.carrier_1_current_values[q] *
																scratch.normals[q] *
																scratch.carrier_fe_face_values.JxW(q);

							hole_current_sum +=
															electron_hole_pair.carrier_2.charge_number *
															scratch.carrier_2_current_values[q] *
															scratch.normals[q] *
															scratch.carrier_fe_face_values.JxW(q);				
						} // for q
					} // Dirichlet
				}  // boundary
			} // for face_no
		} // for cell
	
		electron_current_sum *= sim_params.rescale_current;
		hole_current_sum *= sim_params.rescale_current;
		
		std::cout << "electron current = "
							<< electron_current_sum 
							<< std::endl;

		std::cout << "hole current  = "
							<< hole_current_sum 
							<< std::endl;

		std::cout << "total current = " 
							<< -1.0*(electron_current_sum + hole_current_sum)
							<< std::endl; 
	} // print_semiconductor_current


	template<int dim>
	void 
	SolarCellProblem<dim>::
	print_semiconductor_current_2()
	{	
		double electron_current_sum = 0;
		double hole_current_sum = 0;
	
		Assembly::AssemblyScratch<dim>	scratch(Poisson_fe,
																						carrier_fe,
																						QGauss<dim>(carrier_fe.degree+2),
																						QGauss<dim-1>(carrier_fe.degree+2));

		Assembly::DriftDiffusion::CopyData<dim>	data(carrier_fe);

		typename DoFHandler<dim>::active_cell_iterator
		cell = semiconductor_dof_handler.begin_active(),
		endc = semiconductor_dof_handler.end();

		for(; cell != endc; cell++)
		{
			const unsigned int n_face_q_points = 
														scratch.carrier_fe_face_values.n_quadrature_points;

			const FEValuesExtractors::Vector Current(0);
			const FEValuesExtractors::Scalar Density(dim);

			for(unsigned int face_no=0;
					face_no < GeometryInfo<dim>::faces_per_cell;
					face_no++)
			{
				typename DoFHandler<dim>::face_iterator	face=cell->face(face_no);
				if( face->at_boundary())
				{
					if(face->boundary_id() == Interface)
					{
						scratch.carrier_fe_face_values.reinit(cell, face_no);

						// Get the doping profile values for the boundary conditions
						electrons_e.value_list(
														scratch.carrier_fe_face_values.get_quadrature_points(),
														scratch.carrier_1_bc_values,
														dim); // calls the density values of the donor profile
																	// not the current ones
						holes_e.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_2_bc_values,
													dim); // calls the density values of the donor profile
																// not the current ones
					
	
						// get the electrona and hole densities at the pevious time step
						scratch.carrier_fe_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.electron_interface_values);

						scratch.carrier_fe_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_2.solution,
																				scratch.hole_interface_values);

						// get neighboring selectrolyte cell
						unsigned interface_index = semi_interface_map[
												std::pair<unsigned int, unsigned>(cell->level(),
																													cell->index())];
						typename DoFHandler<dim>::active_cell_iterator
										neighbor_cell(&electrolyte_triangulation,
																	elec_interface_cells[interface_index].first, // level
																	elec_interface_cells[interface_index].second, // index
																	&electrolyte_dof_handler);

						// update fe_neigbors_face_values for this face
						scratch.carrier_fe_neighbor_face_values.reinit(neighbor_cell,
																												 elec_interface_faces[interface_index]);

						// get the reductant and oxidant densities at the pevious time step
						scratch.carrier_fe_neighbor_face_values[Density].get_function_values(
																				redox_pair.carrier_1.solution,
																				scratch.reductant_interface_values);

						scratch.carrier_fe_neighbor_face_values[Density].get_function_values(
																				redox_pair.carrier_2.solution,
																				scratch.oxidant_interface_values);
		
						// loop over all the quadrature points on this face
						for(unsigned int q=0; q<n_face_q_points; q++)
						{
							// copy over the test functions
							electron_current_sum +=
																  sim_params.rescaled_k_et *
																 (scratch.electron_interface_values[q]
																	-
																	scratch.carrier_1_bc_values[q]) * 
																 scratch.oxidant_interface_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q);				
					
							// int_{\Sigma}  v^{-} kht (rho_p - rho_p^e) rho_r ds
							hole_current_sum   +=  
																 sim_params.rescaled_k_ht *
																 (scratch.hole_interface_values[q]
																	-
																 scratch.carrier_2_bc_values[q]) *
																 scratch.reductant_interface_values[q] *
															 	 scratch.carrier_fe_face_values.JxW(q);						
						} // for q
					} // Dirichlet
				}  // boundary
			} // for face_no
		} // for cell
	
		electron_current_sum *= sim_params.rescale_current;
		hole_current_sum *= sim_params.rescale_current;
		
		std::cout << "electron current 2 = "
							<< electron_current_sum 
							<< std::endl;

		std::cout << "hole current 2 = "
							<< hole_current_sum 
							<< std::endl;

		std::cout << "total current 2 = " 
							<< -1.0*(electron_current_sum - hole_current_sum)
							<< std::endl; 
	} // print_semiconductor_current

	/*---------------------------------------------------------------------*/
	/*															RUN																		 */
	/*---------------------------------------------------------------------*/

	template<int dim>
	void
	SolarCellProblem<dim>::
	run_full_system()
	{
		full_system = true;

		TimerOutput	timer(std::cout,
										 TimerOutput::summary,
										 TimerOutput::wall_times);

		timer.enter_subsection("Make Grids");
	
		// make the triangulations
		Grid_Maker::Grid<dim> grid_maker(sim_params);

		// make the grids and set the boundary conditions
		grid_maker.make_grids(semiconductor_triangulation,
													electrolyte_triangulation,
													Poisson_triangulation,
													full_system);

//		grid_maker.print_grid(Poisson_triangulation,
//													"Grid.eps");

		timer.leave_subsection("Make Grids");

		// allocate the memory
		timer.enter_subsection("Allocate Memory");
		setup_dofs();
		timer.leave_subsection("Allocate Memory");

		
		timer.enter_subsection("Build Mappings");
		setup_mappings();
		timer.leave_subsection("Build Mappings");

		print_sim_info();

		// dont remove
		electron_hole_pair.penalty = 1.0; 
		redox_pair.penalty = 1.0;
	
		
		// assemble the global matrices
		timer.enter_subsection("Assemble Poisson Matrix");
		assemble_Poisson_matrix();
		timer.leave_subsection("Assemble Poisson Matrix");
	
		//assign delta_t
		delta_t = sim_params.delta_t;
		std::cout << "delta_t = " << delta_t << std::endl;

		// the 1.0 means we are calculating transients -> 1.0/delta_t
		timer.enter_subsection("Assemble LDG Matrices");
		assemble_LDG_system(1.0);
		timer.leave_subsection("Assemble LDG Matrices");

		// factor the matrices
		timer.enter_subsection("Factor Matrices");
		set_solvers();
		timer.leave_subsection("Factor Matrices");
	


		// make time stepping stuff... 
		unsigned int counter 					= 0;
		unsigned int period						= 1000;
		unsigned int time_step_number;
		double time;

		unsigned int number_outputs = sim_params.time_stamps;
		std::vector<double>						timeStamps(number_outputs);


		// if we are restarting the simulation, read in the dofs
		// from the end of hte last sim as the initial conditions of this one
		if(sim_params.restart_status)
		{
			electron_hole_pair.read_dofs();
			redox_pair.read_dofs();

			// make the time stamps
			for(unsigned int i=0; i<number_outputs; i++)
			{
				timeStamps[i] = sim_params.t_end 
												+ ( (i+1) * (sim_params.t_end_2 - sim_params.t_end)
												 / number_outputs);	
			}
			time_step_number	= number_outputs;
			time 							= sim_params.t_end;
			std::cout << "running until t = " << sim_params.t_end_2 << std::endl;
		}
		else // use defined initial condition functions
		{
			// Get the initial conditions
			VectorTools::project(semiconductor_dof_handler,
													 electron_hole_pair.constraints,
													 QGauss<dim>(degree+1),
													 electrons_e,
													 electron_hole_pair.carrier_1.solution);

			VectorTools::project(semiconductor_dof_handler,
													 electron_hole_pair.constraints,
													 QGauss<dim>(degree+1),
													 holes_e,
													 electron_hole_pair.carrier_2.solution);

			VectorTools::project(electrolyte_dof_handler,
													 redox_pair.constraints,
													 QGauss<dim>(degree+1),
													 reductants_e,
													 redox_pair.carrier_1.solution);

			VectorTools::project(electrolyte_dof_handler,
													 redox_pair.constraints,
													 QGauss<dim>(degree+1),
													 oxidants_e,
													 redox_pair.carrier_2.solution);

			// make the time stamps
			for(unsigned int i=0; i<number_outputs; i++)
				timeStamps[i] = (i+1) * sim_params.t_end / number_outputs;	

			time_step_number	= 0;
			time							= 0.0;
			std::cout << "running until t = " << sim_params.t_end << std::endl;
		}

		// get the intitial potential and electric field
		assemble_Poisson_rhs();
		solve_Poisson();
	
	
		// print the initial values
		print_results(time_step_number);
	 	time_step_number++;	

		// for testing convergence to steady state
		diff_electrons = electron_hole_pair.carrier_1.solution;
		diff_holes =  electron_hole_pair.carrier_2.solution;
		diff_reductants = redox_pair.carrier_1.solution;
		diff_oxidants =  redox_pair.carrier_2.solution;
		diff_Poisson	=	Poisson_object.solution;

		// time stepping until the semiconductor converges at time 
		// t = t_end_1
		for(unsigned int k = 0;
				k < number_outputs;
				k++)
		{				
			// while time < next print stamp time
			while(time < timeStamps[k])
			{
				timer.enter_subsection("Assemble semiconductor rhs");
				assemble_semiconductor_rhs();
				timer.leave_subsection("Assemble semiconductor rhs");

				timer.enter_subsection("Assemble electrolyte rhs");
				assemble_electrolyte_rhs();
				timer.leave_subsection("Assemble electrolyte rhs");

				timer.enter_subsection("Solve LDG Systems");
				solve_full_system();
				timer.leave_subsection("Solve LDG Systems");
	
				timer.enter_subsection("Assemble Poisson rhs");
				assemble_Poisson_rhs();
				timer.leave_subsection("Assemble Poisson rhs");

				timer.enter_subsection("Solve Poisson system");
				solve_Poisson();			
				timer.leave_subsection("Solve Poisson system");

				time += delta_t;
				counter++;

				// SEE IF CONVERGING
				if(counter % period == 0)
				{
					//print_convergence(time);
					timer.enter_subsection("compute current");
					print_convergence(time);
//					print_semiconductor_current();
					print_semiconductor_current_2();
					timer.leave_subsection("compute current");
				}

			} // while


			// print the results
			timer.enter_subsection("Printing");
			print_results(time_step_number);
			time_step_number++;
			timer.leave_subsection("Printing");
		} // end for

		// print the dofs of the end state
		electron_hole_pair.print_dofs();
		redox_pair.print_dofs();


	} // run

	template<int dim>
	void
	SolarCellProblem<dim>::
	run_Schottky_approx()
	{
		full_system = false;

		TimerOutput	timer(std::cout,
										 TimerOutput::summary,
										 TimerOutput::wall_times);

		timer.enter_subsection("Make Grids");
	
		// make the triangulations
		Grid_Maker::Grid<dim> grid_maker(sim_params);

		// make the grids and set the boundary conditions
		grid_maker.make_grids(semiconductor_triangulation,
													electrolyte_triangulation,
													Poisson_triangulation,
													full_system);

		timer.leave_subsection("Make Grids");

		// allocate the memory
		timer.enter_subsection("Allocate Memory");
		setup_dofs();
		timer.leave_subsection("Allocate Memory");
		
		timer.enter_subsection("Build Mappings");
		setup_mappings();
		timer.leave_subsection("Build Mappings");

		print_sim_info();

		// dont remove
		electron_hole_pair.penalty = 1.0; 
		redox_pair.penalty = 1.0;
	
		
		// assemble the global matrices
		timer.enter_subsection("Assemble Poisson Matrix");
		assemble_Poisson_matrix();
		timer.leave_subsection("Assemble Poisson Matrix");
	
		//assign delta_t
		delta_t = sim_params.delta_t;
		std::cout << "delta_t= " << delta_t << std::endl;

		// the 1.0 means we are calculating transients -> 1.0/delta_t
		timer.enter_subsection("Assemble LDG Matrices");
		assemble_LDG_system(1.0);
		timer.leave_subsection("Assemble LDG Matrices");

		// factor the matrices
		timer.enter_subsection("Factor Matrices");
		set_solvers();
		timer.leave_subsection("Factor Matrices");
		
		// Get the initial conditions
		VectorTools::project(semiconductor_dof_handler,
												 electron_hole_pair.constraints,
												 QGauss<dim>(degree+1),
												 electrons_e,
												 electron_hole_pair.carrier_1.solution);

		VectorTools::project(semiconductor_dof_handler,
												 electron_hole_pair.constraints,
												 QGauss<dim>(degree+1),
												 holes_e,
												 electron_hole_pair.carrier_2.solution);

		// get the intitial potential and electric field
		assemble_Poisson_rhs();

		solve_Poisson();
		
		// make time stepping stuff... 
		unsigned int time_step_number	= 0;
		double time										= 0.0;
		unsigned int counter 					= 0;
		unsigned int period						= 1000;

		unsigned int number_outputs = sim_params.time_stamps;
		std::vector<double>						timeStamps(number_outputs);
		
		for(unsigned int i=0; i<number_outputs; i++)
			timeStamps[i] = (i+1) * sim_params.t_end / number_outputs;
	
		std::cout << "running until t = " << sim_params.t_end << std::endl;

		// time stepping until the semiconductor converges at time 
		// t = t_end_1
		for(time_step_number=1;
				time_step_number < number_outputs;
				time_step_number++)
		{				
			std::cout << "time = " << time << std::endl;

			// while time < next print stamp time
			while(time < timeStamps[time_step_number])
			{
				timer.enter_subsection("Assemble semiconductor rhs");
				assemble_Schottky_rhs();
				timer.leave_subsection("Assemble semiconductor rhs");

				timer.enter_subsection("Solve LDG Systems");
				solve_semiconductor_system();
				timer.leave_subsection("Solve LDG Systems");
	
				timer.enter_subsection("Assemble Poisson rhs");
				assemble_Poisson_rhs();
				timer.leave_subsection("Assemble Poisson rhs");

				timer.enter_subsection("Solve Poisson system");
				solve_Poisson();			
				timer.leave_subsection("Solve Poisson system");

				time += delta_t;
				counter++;

				// SEE IF CONVERGING
				if(counter % period == 0)
				{
					
				}

			} // while


			// print the results
			timer.enter_subsection("Printing");
			print_results(time_step_number);
//			print_results_on_boundary(time_step_number);
			timer.leave_subsection("Printing");
		} // end for
	} // run_Schottk_approx

///////////////////////////////////////////////////////////////////////////////
//												TESTING ROUTINES
///////////////////////////////////////////////////////////////////////////////	

	/*----------------------------------------------------------------------*/
	// 										COUPLED MIXED-LDG TEST
	/*----------------------------------------------------------------------*/
	template<int dim>
	void
	SolarCellProblem<dim>::
	assemble_local_coupled_Poisson_test_rhs(
								const typename DoFHandler<dim>::active_cell_iterator & cell,
								Assembly::AssemblyScratch<dim>											 & scratch,
								Assembly::Poisson::CopyData<dim>										 & data,
								const double																				 & time)	
	{
		const unsigned int	dofs_per_cell 	=	
												scratch.Poisson_fe_values.dofs_per_cell;

		const unsigned int 	n_q_points 			= 
												scratch.Poisson_fe_values.n_quadrature_points;

		const unsigned int	n_face_q_points =	
												scratch.Poisson_fe_face_values.n_quadrature_points;

		// SET RHS TIME
		Mixed_Assembler.test_coupling_Poisson_rhs.set_time(time);

		// Get the actual values for vector field and potential from FEValues
		// Use Extractors instead of having to deal with shapefunctions directly
		const FEValuesExtractors::Vector VectorField(0); // Vector as in Vector field
		const FEValuesExtractors::Scalar Potential(dim);
		const FEValuesExtractors::Scalar Density(dim);

		// reset the local_rhs vector to be zero
		data.local_rhs=0;

		// get the Poisson cell info
		std::pair<unsigned int, unsigned int> Poisson_cell_info = 
					s_2_p_map[std::pair<unsigned int, unsigned int>(cell->level(),
																													cell->index())];

		// get the Poisson cell
		typename DoFHandler<dim>::active_cell_iterator
								Poisson_cell(&Poisson_triangulation,
													Poisson_cell_info.first,
													Poisson_cell_info.second,
													&Poisson_dof_handler);

		Poisson_cell->get_dof_indices(data.local_dof_indices);
		scratch.Poisson_fe_values.reinit(Poisson_cell);
		// get the test rhs for poisson
		// Assemble the right hand side on semiconductor side
		scratch.carrier_fe_values.reinit(cell);
			
		// get the electron density values at the previous time step
		scratch.carrier_fe_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.old_carrier_1_density_values);
		// get the test rhs for poisson
		Mixed_Assembler.test_coupling_Poisson_rhs.value_list(
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
				// = -int_{Omega_{e}} (C - u)
				data.local_rhs(i) += -psi_i_potential
															* ( scratch.Poisson_rhs_values[q]	
																- 	
																	scratch.old_carrier_1_density_values[q]
																)	
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
			typename DoFHandler<dim>::face_iterator 	face = Poisson_cell->face(face_no);

			// apply Dirichlet boundary conditions.. 
			// Since we are in a semicondcutor cell we know to apply the 
			// biases
			if((face->at_boundary()) && (face->boundary_id() == Dirichlet))
			{	
				// get the values of the shape functions at this boundary face
				scratch.Poisson_fe_face_values.reinit(Poisson_cell,face_no);
			
				// get the values of the dirichlet boundary conditions evaluated
				// on the quadrature points of this face	
				Mixed_Assembler.test_Poisson_bc.value_list(
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

	} // assemble_local_coupled_Poisson_test_rhs

	template<int dim>
	void
	SolarCellProblem<dim>::
	assemble_local_coupled_DD_test_rhs(
								const typename DoFHandler<dim>::active_cell_iterator & cell,
								Assembly::AssemblyScratch<dim>											 & scratch,
								Assembly::DriftDiffusion::CopyData<dim>							 & data,
								const double																				 & time,
								const double																				 & penalty)	
	{

		const unsigned int dofs_per_cell			 = 
															scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points					 = 
																scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points		 = 
																scratch.carrier_fe_face_values.n_quadrature_points;


		cell->get_dof_indices(data.local_dof_indices);

		// get the Poisson cell info
		std::pair<unsigned int, unsigned int> Poisson_cell_info = 
					s_2_p_map[std::pair<unsigned int, unsigned int>(cell->level(),
																													cell->index())];

		// get the Poisson cell
		typename DoFHandler<dim>::active_cell_iterator
								Poisson_cell(&Poisson_triangulation,
													Poisson_cell_info.first,
													Poisson_cell_info.second,
													&Poisson_dof_handler);

		scratch.Poisson_fe_values.reinit(Poisson_cell);

		scratch.carrier_fe_values.reinit(cell);
		
		// reset the local_rhs to be zero
		data.local_carrier_1_rhs=0;
		data.local_carrier_2_rhs=0;
	
		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);
		const FEValuesExtractors::Vector ElectricField(0);

		LDG_Assembler.test_DD_rhs.set_time(time);
		LDG_Assembler.test_DD_bc.set_time(time);

		// rhs values
		LDG_Assembler.test_DD_rhs.value_list(
												scratch.carrier_fe_values.get_quadrature_points(),
												scratch.generation_values);
	
		// get the values of carrier_1 density
		scratch.carrier_fe_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.old_carrier_1_density_values);


		// get the electric field values at the previous time step
		scratch.Poisson_fe_values[ElectricField].get_function_values(
																									Poisson_object.solution,
																									scratch.electric_field_values);


		double h = cell->diameter();
		// loop over all the quadrature points in this cell and compute body integrals
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// loop over all the test function dofs and get the test functions
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				const double					psi_i_density	= 
																	scratch.carrier_fe_values[Density].value(i,q);
		
				const Tensor<1,dim>		psi_i_field		=
																	scratch.carrier_fe_values[Current].value(i,q);
					
				// contribution from RHS function + Drift
				// int_{Omega} v * R +  p * E * u dx
				data.local_carrier_1_rhs(i) += ( 
													(psi_i_density * scratch.generation_values[q])
													-
													psi_i_field *
													scratch.electric_field_values[q] *
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
					LDG_Assembler.test_DD_bc.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_1_bc_values);

					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// get the test function
							const Tensor<1, dim>  psi_i_field	= 
													scratch.carrier_fe_face_values[Current].value(i,q);
							const double 					psi_i_density = 
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
				}
				else if(face->boundary_id() == Neumann)
				{
					// NOTHIN TO DO IF INSULATING
				}
				else
					Assert(false, ExcNotImplemented() );
			} // end at boundary
		} // for face_no
	
	} // assemble_local_coupling_DD_rhs

	/*-----------------------------------------------------------------------*/
	/*											INTERFACE ASSMBLY ROUTINES											 */
	/*-----------------------------------------------------------------------*/

	template<int dim>
	void 
	SolarCellProblem<dim>::
	assemble_local_test_semiconductor_rhs(
								const typename DoFHandler<dim>::active_cell_iterator & cell,
								Assembly::AssemblyScratch<dim>											 & scratch,
								Assembly::DriftDiffusion::CopyData<dim>							 & data,
								const double 																				 & time,
								const double 																				 & penalty)									
	{
		const unsigned int dofs_per_cell			 = 
															scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points					 = 
																scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points		 = 
																scratch.carrier_fe_face_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		// reinitialize the fe_values 
		scratch.carrier_fe_values.reinit(cell);
		
		// reset the local_rhs to be zero
		data.local_carrier_1_rhs=0;
		data.local_carrier_2_rhs=0;
	
		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);
		Tensor<1,dim>							field_values;
		field_values[0]	= 0.0; 
		field_values[1] = 0.0; 

		LDG_Assembler.test_interface_rhs.set_time(time);
		LDG_Assembler.test_interface_bc.set_time(time);
		LDG_Assembler.test_interface_function.set_time(time);

		// the test version	
		LDG_Assembler.test_interface_rhs.value_list(
												scratch.carrier_fe_values.get_quadrature_points(),
												scratch.generation_values);
	
		// get the values of carrier_1 and carrier_2 densities at the pevious time step
		scratch.carrier_fe_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.old_carrier_1_density_values);


		double h = cell->diameter();
		// loop over all the quadrature points in this cell and compute body integrals
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// loop over all the test function dofs and get the test functions
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				const double					psi_i_density	= 
																	scratch.carrier_fe_values[Density].value(i,q);
		
					const Tensor<1,dim>		psi_i_field		=
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
					LDG_Assembler.test_interface_bc.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_1_bc_values);

					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// get the test function
							const Tensor<1, dim>  psi_i_field	= 
													scratch.carrier_fe_face_values[Current].value(i,q);
							const double 					psi_i_density = 
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
					// get this carriers density values on the interface
					scratch.carrier_fe_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.electron_interface_values);

					// get the interface functions value on the interface
					LDG_Assembler.test_interface_function.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_2_bc_values);
	
					// get neighboring cell
					unsigned interface_index = semi_interface_map[
												std::pair<unsigned int, unsigned>(cell->level(),
																													cell->index())];
					typename DoFHandler<dim>::active_cell_iterator
										neighbor_cell(&electrolyte_triangulation,
																	elec_interface_cells[interface_index].first, // level
																	elec_interface_cells[interface_index].second, // index
																	&electrolyte_dof_handler);

					// update fe_neigbors_face_values for this face
					scratch.carrier_fe_neighbor_face_values.reinit(neighbor_cell,
																												 elec_interface_faces[interface_index]);

					// get the other carriers density values on this face
					scratch.carrier_fe_neighbor_face_values[Density].get_function_values(
																				redox_pair.carrier_1.solution,
																				scratch.oxidant_interface_values);	
			
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
															(
															scratch.electron_interface_values[q]  *
															scratch.oxidant_interface_values[q]
															-
															scratch.carrier_2_bc_values[q]
															) *
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
	}

	
	template<int dim>
	void 
	SolarCellProblem<dim>::
	assemble_local_test_electrolyte_rhs(
								const typename DoFHandler<dim>::active_cell_iterator & cell,
								Assembly::AssemblyScratch<dim>											 & scratch,
								Assembly::DriftDiffusion::CopyData<dim>							 & data,
								const double 																				 & time,
								const double 																				 & penalty)									
	{
		const unsigned int dofs_per_cell			 = 
															scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points					 = 
																scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points		 = 
																scratch.carrier_fe_face_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		// reinitialize the fe_values 
		scratch.carrier_fe_values.reinit(cell);
		
		// reset the local_rhs to be zero
		data.local_carrier_1_rhs=0;
		data.local_carrier_2_rhs=0;
	
		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);
		Tensor<1,dim>							field_values;
		field_values[0]	= 0.0; 
		field_values[1] = 0.0; 

		LDG_Assembler.test_interface_rhs.set_time(time);
		LDG_Assembler.test_interface_bc.set_time(time);
		LDG_Assembler.test_interface_function.set_time(time);

		// the test version	
		LDG_Assembler.test_interface_rhs.value_list(
												scratch.carrier_fe_values.get_quadrature_points(),
												scratch.generation_values);
	
		// get the values of carrier_1 and carrier_2 densities at the pevious time step
		scratch.carrier_fe_values[Density].get_function_values(
																				redox_pair.carrier_1.solution,
																				scratch.old_carrier_1_density_values);

		double h = cell->diameter();
		// loop over all the quadrature points in this cell and compute body integrals
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// loop over all the test function dofs and get the test functions
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				const double					psi_i_density	= 
																	scratch.carrier_fe_values[Density].value(i,q);
		
				const Tensor<1,dim>		psi_i_field		=
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
					LDG_Assembler.test_interface_bc.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_1_bc_values);

					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// get the test function
							const Tensor<1, dim>  psi_i_field	= 
													scratch.carrier_fe_face_values[Current].value(i,q);
							const double 					psi_i_density = 
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
					// get this carriers density value on the interface
					scratch.carrier_fe_face_values[Density].get_function_values(
																				redox_pair.carrier_1.solution,
																				scratch.oxidant_interface_values);

					// get the interface function value on the interface
					LDG_Assembler.test_interface_function.value_list(
													scratch.carrier_fe_face_values.get_quadrature_points(),
													scratch.carrier_2_bc_values);

						// get neighboring cell
					unsigned interface_index = elec_interface_map[
												std::pair<unsigned int, unsigned>(cell->level(),
																													cell->index())];
					typename DoFHandler<dim>::active_cell_iterator
										neighbor_cell(&semiconductor_triangulation,
																	semi_interface_cells[interface_index].first, // level
																	semi_interface_cells[interface_index].second, // index
																	&semiconductor_dof_handler);

					// update fe_neigbors_face_values for this face
					scratch.carrier_fe_neighbor_face_values.reinit(neighbor_cell,
																												 semi_interface_faces[interface_index]);

					// get the other carriers density values on this face
					scratch.carrier_fe_neighbor_face_values[Density].get_function_values(
																				electron_hole_pair.carrier_1.solution,
																				scratch.electron_interface_values);	
			
			
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
															(
															scratch.electron_interface_values[q] *
															scratch.oxidant_interface_values[q]
															-
															scratch.carrier_2_bc_values[q]
															) * 
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
	}




	/*---------------------------------------------------------------------*/
	/*												STEADY STATE TEST														 */
	/*---------------------------------------------------------------------*/

	template<int dim>
	void
	SolarCellProblem<dim>::
	test_steady_state(const unsigned int & n_refine,
										ConvergenceTable & Mixed_table,
										ConvergenceTable & LDG_table)
	{
		full_system = false;

		double 	primary_error, flux_error;
		sim_params.set_params_for_testing(n_refine);

		Grid_Maker::Grid<dim> grid_maker(sim_params);
		grid_maker.make_test_grid(Poisson_triangulation, n_refine);
		grid_maker.make_test_grid(semiconductor_triangulation, n_refine);
//		grid_maker.refine_test_grid(Poisson_triangulation,2);
//		grid_maker.refine_test_grid(semiconductor_triangulation, 2);

		setup_dofs();
		electron_hole_pair.set_semiconductor_for_testing(sim_params);

		assemble_Poisson_matrix();
		assemble_LDG_system(0.0);

		set_solvers();

		unsigned int n_active_cells = Poisson_triangulation.n_active_cells();
		unsigned int n_dofs					= Poisson_dof_handler.n_dofs();
		double h =  GridTools::maximal_cell_diameter(Poisson_triangulation);
		Mixed_table.add_value("h", h);		
		Mixed_table.add_value("cells", n_active_cells);
		Mixed_table.add_value("dofs", n_dofs);

		n_active_cells = semiconductor_triangulation.n_active_cells();
		n_dofs				 = semiconductor_dof_handler.n_dofs();

		LDG_table.add_value("h", h);		
		LDG_table.add_value("cells", n_active_cells);
		LDG_table.add_value("dofs", n_dofs);
		
		// assemble poisson rhs
		WorkStream::run(Poisson_dof_handler.begin_active(),
										Poisson_dof_handler.end(),
										std_cxx11::bind(&MixedPoisson::MixedFEM<dim>::
																		assemble_local_test_rhs,
																		Mixed_Assembler, // the object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_Poisson_rhs,
																		this, // this object
																		std_cxx11::_1),
										Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
										Assembly::Poisson::CopyData<dim>(Poisson_fe)
										);

		//assemble ldg rhs
		WorkStream::run(semiconductor_dof_handler.begin_active(),
										semiconductor_dof_handler.end(),
										std_cxx11::bind(&LDG_System::LDG<dim>::
																		assemble_local_test_rhs,
																		LDG_Assembler, // Assembler object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		electron_hole_pair.penalty),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_semiconductor_system_rhs,
																		this, // this object
																		std_cxx11::_1),
										Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
										Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
										);

		solve_Poisson();
		solve_semiconductor_system();

		Mixed_Assembler.compute_errors(Poisson_triangulation, 
																	 Poisson_dof_handler,
																	 Poisson_object.solution,
																	 primary_error,
																	 flux_error);
	
		Mixed_table.add_value("Phi", primary_error);
		Mixed_table.add_value("D", flux_error);

	
		LDG_Assembler.compute_errors(semiconductor_triangulation,
																 semiconductor_dof_handler,
																 electron_hole_pair.carrier_1.solution,
																 primary_error,
																 flux_error,
																 true,
																 0.0); // no time

		LDG_table.add_value("u", primary_error);
		LDG_table.add_value("J", flux_error);
		
		Mixed_Assembler.output_unscaled_results(Poisson_dof_handler,
																						Poisson_object.solution,
																						n_refine);

		LDG_Assembler.output_unscaled_results(semiconductor_dof_handler,
																					electron_hole_pair,
																					n_refine);
	}

	/*---------------------------------------------------------------------*/
	/*												LDG-IMEX TEST																 */
	/*---------------------------------------------------------------------*/
	template<int dim>
	void
	SolarCellProblem<dim>::
	test_transient(const unsigned int & n_refine,
							   ConvergenceTable		& LDG_table)
	{
		full_system = false;
		double 	primary_error, flux_error;

		sim_params.set_params_for_testing(n_refine);

		Grid_Maker::Grid<dim> grid_maker(sim_params);
		grid_maker.make_test_grid(Poisson_triangulation, n_refine);
		grid_maker.make_test_tran_grid(semiconductor_triangulation, n_refine);
	//	grid_maker.refine_test_grid(semiconductor_triangulation, 2);

		setup_dofs();

		double h = GridTools::maximal_cell_diameter(semiconductor_triangulation);

		electron_hole_pair.set_semiconductor_for_testing(sim_params);
 		
		delta_t = 1.0;
		for(unsigned int i=0; i < carrier_fe.degree+1; i++)
			delta_t *= h;

		std::cout << "dt = " << delta_t << std::endl;

		assemble_LDG_system(1.0);

		electron_hole_pair.carrier_1.set_solver();

		// Get the initial conditions
		VectorTools::project(semiconductor_dof_handler,
												 electron_hole_pair.constraints,
												 QGauss<dim>(degree+1),
												 LDG_Assembler.test_interface_initial,
												 electron_hole_pair.carrier_1.solution);

		double t_end = 1.0;
		double time  = 0.0;
		
		while(time < t_end)
		{
			// set carrier_system_rhs = M * u^{n-1}
			electron_hole_pair.mass_matrix.vmult(electron_hole_pair.carrier_1.system_rhs,
																				 electron_hole_pair.carrier_1.solution);
			// set rhs for time step n-1
			WorkStream::run(semiconductor_dof_handler.begin_active(),
											semiconductor_dof_handler.end(),
											std_cxx11::bind(&LDG_System::LDG<dim>::
																		assemble_local_test_transient_rhs,
																		LDG_Assembler, // Assembler object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		electron_hole_pair.carrier_1.solution,
																		time,
																		electron_hole_pair.penalty),
											std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_semiconductor_system_rhs,
																		this, // this object
																		std_cxx11::_1),
											Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
											Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
											);
			electron_hole_pair.carrier_1.solve();
	
			// UPDATE TIME
			time += delta_t;
		}	// end while																		


		// LDG ERRORS
		LDG_Assembler.compute_errors(semiconductor_triangulation,
																 semiconductor_dof_handler,
																 electron_hole_pair.carrier_1.solution,
																 primary_error,
																 flux_error,
																 false,
																 time);

		unsigned int n_active_cells = semiconductor_triangulation.n_active_cells();
		unsigned int n_dofs				 	= semiconductor_dof_handler.n_dofs();

	
		LDG_table.add_value("h",h);
		LDG_table.add_value("cells", n_active_cells);
		LDG_table.add_value("dofs", n_dofs);

		LDG_table.add_value("u", primary_error);
		LDG_table.add_value("J", flux_error);
		
		LDG_Assembler.output_unscaled_results(semiconductor_dof_handler,
																					electron_hole_pair,
																					n_refine);

	} // end test_transient


	/*---------------------------------------------------------------------*/
	/*												LDG-MIXED TEST															 */
	/*---------------------------------------------------------------------*/

	template<int dim>
	void
	SolarCellProblem<dim>::
	test_DD_Poisson(const unsigned int & n_refine,
							   ConvergenceTable		& Mixed_table,
							   ConvergenceTable		& LDG_table)
	{
		full_system = false;
		double 	primary_error, flux_error;

		sim_params.set_params_for_testing(n_refine);

		Grid_Maker::Grid<dim> grid_maker(sim_params);
		grid_maker.make_DD_Poisson_grid(Poisson_triangulation, n_refine);
		grid_maker.make_DD_Poisson_grid(semiconductor_triangulation, n_refine);
	//	grid_maker.refine_test_grid(semiconductor_triangulation, 2);

		setup_dofs();
		setup_mappings();

		double h = GridTools::maximal_cell_diameter(semiconductor_triangulation);

		electron_hole_pair.set_semiconductor_for_testing(sim_params);
 		
		delta_t = 1.0;
		for(unsigned int i=0; i < carrier_fe.degree+1; i++)
			delta_t *= h;

		std::cout << "dt = " << delta_t << std::endl;

		assemble_Poisson_matrix();
		assemble_LDG_system(1.0);

		electron_hole_pair.carrier_1.set_solver();

		// Get the initial conditions
		VectorTools::project(semiconductor_dof_handler,
												 electron_hole_pair.constraints,
												 QGauss<dim>(degree+1),
												 LDG_Assembler.test_interface_initial,
												 electron_hole_pair.carrier_1.solution);
		
		set_solvers();

		double t_end = 1.0;
		double time  = 0.0;
	
		while(time < t_end)
		{
			///////////////////////////////////////
			// UPDATE POISSON RHS AND SOLVE
			/////////////////////////////////////////
			Poisson_object.system_rhs = 0;

			WorkStream::run(semiconductor_dof_handler.begin_active(),
											semiconductor_dof_handler.end(),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		assemble_local_coupled_Poisson_test_rhs,
																		this, // this object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		time),
										std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_Poisson_rhs,
																		this, // this object
																		std_cxx11::_1),
										Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
										Assembly::Poisson::CopyData<dim>(Poisson_fe)
										);
	
			solve_Poisson();

			///////////////////////////////////////
			// UPDATE DD RHS AND SOLVE
			/////////////////////////////////////////
			// set carrier_system_rhs = M * u^{n-1}
			electron_hole_pair.mass_matrix.vmult(electron_hole_pair.carrier_1.system_rhs,
																				 electron_hole_pair.carrier_1.solution);
		
			// set rhs for time step n-1
			WorkStream::run(semiconductor_dof_handler.begin_active(),
											semiconductor_dof_handler.end(),
											std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		assemble_local_coupled_DD_test_rhs,		
																		this, // this object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		time,
																		electron_hole_pair.penalty),
											std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_semiconductor_system_rhs,
																		this, // this object
																		std_cxx11::_1),
											Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
											Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
											);

			electron_hole_pair.carrier_1.solve();

			// UPDATE TIME
			time += delta_t;
		}	// end while																		

		// LDG ERRORS
		LDG_Assembler.compute_coupled_errors(semiconductor_triangulation,
																 semiconductor_dof_handler,
																 electron_hole_pair.carrier_1.solution,
																 primary_error,
																 flux_error,
																 time);

		unsigned int n_active_cells = semiconductor_triangulation.n_active_cells();
		unsigned int n_dofs				 = semiconductor_dof_handler.n_dofs();

		LDG_table.add_value("h",h);
		LDG_table.add_value("cells", n_active_cells);
		LDG_table.add_value("dofs", n_dofs);

		LDG_table.add_value("u", primary_error);
		LDG_table.add_value("J", flux_error);
		
		// MIXED ERRORS
		Mixed_Assembler.compute_errors(Poisson_triangulation, 
																	 Poisson_dof_handler,
																	 Poisson_object.solution,
																	 primary_error,
																	 flux_error);
		
		n_active_cells = Poisson_triangulation.n_active_cells();
		n_dofs				 = Poisson_dof_handler.n_dofs();

		Mixed_table.add_value("h",h);
		Mixed_table.add_value("cells", n_active_cells);
		Mixed_table.add_value("dofs", n_dofs);

		Mixed_table.add_value("Phi", primary_error);
		Mixed_table.add_value("D", flux_error);
	}

	/*---------------------------------------------------------------------*/
	/*												LDG-INTERFACE TEST													 */
	/*---------------------------------------------------------------------*/

	template<int dim>
	void
	SolarCellProblem<dim>::
	test_interface_coupling(const unsigned int & n_refine,
							 						ConvergenceTable	 & LDG_table)
	{
		full_system = true;

		double 	primary_error, flux_error;
		sim_params.set_params_for_testing(n_refine);

		// make the triangulations
		Grid_Maker::Grid<dim> grid_maker(sim_params);

		// make the grids
		grid_maker.make_grids(semiconductor_triangulation,
													electrolyte_triangulation,
													Poisson_triangulation,
													true);

		setup_mappings();
						
		setup_dofs();

		double h = GridTools::maximal_cell_diameter(semiconductor_triangulation);

		electron_hole_pair.set_semiconductor_for_testing(sim_params);
		redox_pair.set_electrolyte_for_testing(sim_params);
 		
		delta_t = 1.0;
		for(unsigned int i=0; i < carrier_fe.degree+1; i++)
			delta_t *= h;

		std::cout << "dt = " << delta_t << std::endl;

		assemble_LDG_system(1.0);

		electron_hole_pair.carrier_1.set_solver();
		redox_pair.carrier_1.set_solver();

		// Get the initial conditions
		VectorTools::project(semiconductor_dof_handler,
												 electron_hole_pair.constraints,
												 QGauss<dim>(degree+1),
												 LDG_Assembler.test_interface_initial,
												 electron_hole_pair.carrier_1.solution);

		// Get the initial conditions
		VectorTools::project(electrolyte_dof_handler,
												 redox_pair.constraints,
												 QGauss<dim>(degree+1),
												 LDG_Assembler.test_interface_initial,
												 redox_pair.carrier_1.solution);

		double t_end = 1.0;
		double time  = 0.0;
	
		while(time < t_end)
		{
			// set carrier_system_rhs = M * u^{n-1}
			electron_hole_pair.mass_matrix.vmult(electron_hole_pair.carrier_1.system_rhs,
																				 electron_hole_pair.carrier_1.solution);
			// set rhs for time step n-1
			WorkStream::run(semiconductor_dof_handler.begin_active(),
											semiconductor_dof_handler.end(),
											std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		assemble_local_test_semiconductor_rhs,
																		this, // Assembler object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		time,
																		electron_hole_pair.penalty),
											std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_semiconductor_system_rhs,
																		this, // this object
																		std_cxx11::_1),
											Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
											Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
											);

			// set carrier_system_rhs = M * u^{n-1}
			redox_pair.mass_matrix.vmult(redox_pair.carrier_1.system_rhs,
																	redox_pair.carrier_1.solution);
			// set rhs for time step n-1
			WorkStream::run(electrolyte_dof_handler.begin_active(),
											electrolyte_dof_handler.end(),
											std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		assemble_local_test_electrolyte_rhs,
																		this, // Assembler object
																		std_cxx11::_1,
																		std_cxx11::_2,
																		std_cxx11::_3,
																		time,
																		redox_pair.penalty),
											std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
																		copy_local_to_global_electrolyte_system_rhs,
																		this, // this object
																		std_cxx11::_1),
											Assembly::AssemblyScratch<dim>(Poisson_fe,
																									 carrier_fe,
																									 QGauss<dim>(degree+2),
																									 QGauss<dim-1>(degree+2)),
											Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
											);
		
			Threads::TaskGroup<void> task_group;

			task_group += Threads::new_task(&ChargeCarrierSpace::
																		Carrier<dim>::solve,
																		electron_hole_pair.carrier_1);


			task_group += Threads::new_task(&ChargeCarrierSpace::
																			Carrier<dim>::solve,
																			redox_pair.carrier_1);

			task_group.join_all();

			time += delta_t;
		}	// end while																		

		//LDG_Assembler.compute_interface_errors(semiconductor_triangulation,
		LDG_Assembler.compute_interface_errors(semiconductor_triangulation,
																 semiconductor_dof_handler,
																 electron_hole_pair.carrier_1.solution,
																 primary_error,
																 flux_error,
																 time);


		unsigned int n_active_cells = semiconductor_triangulation.n_active_cells();
		unsigned int n_dofs				 = semiconductor_dof_handler.n_dofs();

		LDG_table.add_value("h",h);
		LDG_table.add_value("cells", n_active_cells);
		LDG_table.add_value("dofs", n_dofs);

		LDG_table.add_value("u1", primary_error);
		LDG_table.add_value("J1", flux_error);
	
		LDG_Assembler.compute_interface_errors(electrolyte_triangulation,
																 electrolyte_dof_handler,
																 redox_pair.carrier_1.solution,
																 primary_error,
																 flux_error,
																 time);

		LDG_table.add_value("u2", primary_error);
		LDG_table.add_value("J2", flux_error);

		print_results(n_refine);
	
	} // end test_interface_coupling


} // namespace SOLARCELL
 // namespace SOLARCELL
