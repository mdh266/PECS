#include "../include/Grid.hpp"


/////////////////////////////////////////////////////////////////////////////////
// GRID AND BOUNDARIES
/////////////////////////////////////////////////////////////////////////////////
namespace Grid_Maker
{
	using namespace dealii;

	template<int dim>
	Grid<dim>::
	Grid(const ParameterSpace::Parameters & sim_params)
	{
		scaled_domain_height	= sim_params.scaled_domain_height;
		scaled_domain_length	= sim_params.scaled_domain_length;
		scaled_radius_one		= sim_params.scaled_radius_one;
		scaled_radius_two		= sim_params.scaled_radius_two;
		scaled_boundary_layer	= sim_params.scaled_boundary_layer;
		n_global_refine			= sim_params.n_global_refine;	
		n_local_refine			= sim_params.n_local_refine;
		insulated				= sim_params.insulated;
		schottky				= sim_params.schottky_status;
		
		if(n_local_refine == 0)
		{
			use_boundary_layer = false;
			scaled_boundary_layer = 0.0;
		}
		else if((n_local_refine > 0)
						&&
			 	(scaled_boundary_layer > 0) )
		{
			use_boundary_layer = true;
		}
		else
		{
			std::cerr << "Boundary layer & n_local_refine need to be >= 0\n";	
		}

	}
	

	template<int dim>
	void 
	Grid<dim>::
	make_grids(Triangulation<dim> & semiconductor_triang,
               Triangulation<dim> & electrolyte_triang,
               Triangulation<dim> & Poisson_triang,
               const bool		  & full_system)
	{
		// make grids
		make_semiconductor_grid(semiconductor_triang);
		make_electrolyte_grid(electrolyte_triang);

		if(full_system)
		{
			make_merged_grid(semiconductor_triang,
							 electrolyte_triang,
							 Poisson_triang);
		}
		else
			make_semiconductor_grid(Poisson_triang);


		// globally refine now
		semiconductor_triang.refine_global(n_global_refine);
		electrolyte_triang.refine_global(n_global_refine);
		Poisson_triang.refine_global(n_global_refine);
	
		// locally refine now
		for(unsigned int refine_num=0; refine_num < n_local_refine; refine_num++)
		{
			// semiconductor
			typename Triangulation<dim>::active_cell_iterator
						cell = semiconductor_triang.begin_active(),
						endc = semiconductor_triang.end();
			// loop over all the cells and mark for refinement
			for(; cell != endc; cell++)
				if(cell->material_id() == semi_boundary_layer_id)
					cell->set_refine_flag();
		
			// electrolyte
			cell = electrolyte_triang.begin_active(),
			endc = electrolyte_triang.end();
	
			// loop over all the cells and mark for refinement
			for(; cell != endc; cell++)
				if(cell->material_id() == elec_boundary_layer_id)
					cell->set_refine_flag();

			// Poisson
			cell = Poisson_triang.begin_active(),
			endc = Poisson_triang.end();
	
			// loop over all the cells and mark for refinement
			for(; cell != endc; cell++)
			{
				if((cell->material_id() == semi_boundary_layer_id)
						||
					 (cell->material_id() == elec_boundary_layer_id) )
				{
					cell->set_refine_flag();
				}
			}

			semiconductor_triang.execute_coarsening_and_refinement();
			electrolyte_triang.execute_coarsening_and_refinement();
			Poisson_triang.execute_coarsening_and_refinement();

		} // refine_num
		
		// set Dirichlet boundary conditions
		make_Dirichlet_boundaries(semiconductor_triang);
		make_Dirichlet_boundaries(electrolyte_triang);
		make_Dirichlet_boundaries(Poisson_triang);

		// make insulating conditions.  Handled in input file
		if(insulated)
		{
			make_Neumann_boundaries(Poisson_triang);
			make_Neumann_boundaries(semiconductor_triang);
			make_Neumann_boundaries(electrolyte_triang);
		}

		// make Schottky conditions. 
		if(schottky)
		{
			make_Schottky_boundaries(semiconductor_triang);
			make_Schottky_boundaries(Poisson_triang);
		}

	} // make grids

	template<int dim>
	void
	Grid<dim>::
	make_semiconductor_grid(Triangulation<dim> 	& triangulation)
	{
//		const unsigned int dim = 2;

		// bulk layer vetrices
		static const Point<2> vertices_1[]
			= { Point<2>(0,0),
					Point<2>(scaled_radius_two - scaled_boundary_layer, 0),
					Point<2>(0,scaled_domain_height),
					Point<2>(scaled_radius_one - scaled_boundary_layer, scaled_domain_height)
				};

		const unsigned int n_vertices_1 = sizeof(vertices_1)/sizeof(vertices_1[0]);

		const std::vector<Point<dim>> vertices_list_1(&vertices_1[0],
													  &vertices_1[n_vertices_1]);

		static const int cell_vertices_1[][GeometryInfo<dim>::vertices_per_cell]
			= {  {0,1,2,3} };

		const unsigned int n_cells_1 = sizeof(cell_vertices_1)/sizeof(cell_vertices_1[0]);

		// BULK layer
		std::vector<CellData<dim> > cells_1(n_cells_1, CellData<dim>() );
		for(unsigned int i=0; i<n_cells_1; i++)
		{
			for(unsigned int j=0;
					j<GeometryInfo<dim>::vertices_per_cell;
					j++)
			{
				cells_1[i].vertices[j] = cell_vertices_1[i][j];
			}
			cells_1[i].material_id = semiconductor_id;
		}

	// boundary layer vetrices
	static const Point<2> vertices_2[]
			= { Point<2>(scaled_radius_two - scaled_boundary_layer, 0),
					Point<2>(scaled_radius_two, 0),
					Point<2>(scaled_radius_one - scaled_boundary_layer, scaled_domain_height),
					Point<2>(scaled_radius_one, scaled_domain_height)
				};

		const unsigned int n_vertices_2 = sizeof(vertices_2)/sizeof(vertices_2[0]);

		const std::vector<Point<dim>> vertices_list_2(&vertices_2[0],
													  &vertices_2[n_vertices_2]);
		static const int cell_vertices_2[][GeometryInfo<dim>::vertices_per_cell]
			= {  {0,1,2,3} };

		const unsigned int n_cells_2 = sizeof(cell_vertices_2)/sizeof(cell_vertices_2[0]);
	
			// boundary layer
		std::vector<CellData<dim> > cells_2(n_cells_2, CellData<dim>() );
		for(unsigned int i=0; i<n_cells_2; i++)
		{
			for(unsigned int j=0;
					j<GeometryInfo<dim>::vertices_per_cell;
					j++)
			{
				cells_2[i].vertices[j] = cell_vertices_2[i][j];
			}
			cells_2[i].material_id = semi_boundary_layer_id;
		}

		if(use_boundary_layer)
		{
			// create a temporary bulk layer triangulation
			Triangulation<dim>	temp_bulk_triangulation;
			temp_bulk_triangulation.create_triangulation(vertices_list_1,
														 cells_1,
														 SubCellData());
		
			// create a temporary boundary layer triangulation
			Triangulation<dim>	temp_boundary_layer_triangulation;
			temp_boundary_layer_triangulation.create_triangulation(vertices_list_2,
																														cells_2,
																														SubCellData());

			GridGenerator::merge_triangulations(temp_bulk_triangulation,
												temp_boundary_layer_triangulation,
												triangulation);
		}
		else
			triangulation.create_triangulation(vertices_list_1,
											   cells_1,
											   SubCellData());

	}

	template<int dim>
	void
	Grid<dim>::
	make_electrolyte_grid(Triangulation<dim> 	& triangulation)
	{
//		const unsigned int dim = 2;

		//	boundary layer vertices
		static const Point<2> vertices_1[]
		= {		Point<2>(scaled_radius_two, 0),
					Point<2>(scaled_radius_two + scaled_boundary_layer, 0),
					Point<2>(scaled_radius_one, scaled_domain_height),
					Point<2>(scaled_radius_one + scaled_boundary_layer, scaled_domain_height)
				};

		const unsigned int n_vertices_1 = sizeof(vertices_1)/sizeof(vertices_1[0]);

		const std::vector<Point<dim>> vertices_list_1(&vertices_1[0],
													  &vertices_1[n_vertices_1]);

		static const int cell_vertices_1[][GeometryInfo<dim>::vertices_per_cell]
			= {  {0,1,2,3}};

		const unsigned int n_cells_1 = sizeof(cell_vertices_1)/sizeof(cell_vertices_1[0]);

		// boundary layer 		
		std::vector<CellData<dim> > cells_1(n_cells_1, CellData<dim>() );
		for(unsigned int i=0; i<n_cells_1; i++)
		{
			for(unsigned int j=0;
					j<GeometryInfo<dim>::vertices_per_cell;
					j++)
			{
				cells_1[i].vertices[j] = cell_vertices_1[i][j];
			}
			cells_1[i].material_id = elec_boundary_layer_id;
		}

		// bulk layer vetrices
		static const Point<2> vertices_2[]
		= {		Point<2>(scaled_radius_two + scaled_boundary_layer, 0),
					Point<2>(scaled_domain_length, 0),
					Point<2>(scaled_radius_one + scaled_boundary_layer, scaled_domain_height),
					Point<2>(scaled_domain_length, scaled_domain_height)
				};

		const unsigned int n_vertices_2 = sizeof(vertices_2)/sizeof(vertices_2[0]);


		const std::vector<Point<dim>> vertices_list_2(&vertices_2[0],
													  &vertices_2[n_vertices_2]);


		static const int cell_vertices_2[][GeometryInfo<dim>::vertices_per_cell]
			= {  {0,1,2,3}};

		const unsigned int n_cells_2 = sizeof(cell_vertices_2)/sizeof(cell_vertices_2[0]);

		std::vector<CellData<dim> > cells_2(n_cells_2, CellData<dim>() );
		for(unsigned int i=0; i<n_cells_2; i++)
		{
			for(unsigned int j=0;
					j<GeometryInfo<dim>::vertices_per_cell;
					j++)
			{
				cells_2[i].vertices[j] = cell_vertices_2[i][j];
			}
			cells_2[i].material_id = electrolyte_id;
		}

		if(use_boundary_layer)
		{
			// create boundary layer triangulation
			Triangulation<dim>	temp_boundary_layer_triangulation;
			temp_boundary_layer_triangulation.create_triangulation(vertices_list_1,
																   cells_1,
																   SubCellData());

			// create bulk layer triangulation
			Triangulation<dim>	temp_bulk_triangulation;
			temp_bulk_triangulation.create_triangulation(vertices_list_2,
														 cells_2,
														 SubCellData());
	
			GridGenerator::merge_triangulations(temp_boundary_layer_triangulation,
												temp_bulk_triangulation,
												triangulation);
		}	
		else
			triangulation.create_triangulation(vertices_list_2,
											   cells_2,
											   SubCellData());
	}

	template <int dim>
	void 
	Grid<dim>::
	make_merged_grid(const Triangulation<dim>				& semiconductor_triang,
					 const Triangulation<dim>				& electrolyte_triang,
					 Triangulation<dim>					    & merged_triangulation)
	{
		// merges the two triangulations in to one
		GridGenerator::merge_triangulations(semiconductor_triang,
											electrolyte_triang,
											merged_triangulation);
		// note matrial id's of  merged_triangulation are inherited from
		// semiconductor_triang and electrolyte_triang
	}



	template<int dim> 
	void 
	Grid<dim>::
	make_Dirichlet_boundaries(Triangulation<dim> & triangulation)
	{
		typename Triangulation<dim>::active_cell_iterator
										cell = triangulation.begin_active(),
										endc = triangulation.end();
		// loop over all the cells
		for(; cell != endc; cell++)
		{
			// loop over all the faces of the cell and find which are on the boundary
			for(unsigned int face_no=0;
					face_no < GeometryInfo<dim>::faces_per_cell;
					face_no++)
			{
				if(cell->face(face_no)->at_boundary() )
				{
					//NOTE: Default is 0, which implies Interface

					//sets only the non-interface boundary to be Dirichlet	
					if((cell->face(face_no)->center()[0] == 0.0) ||
						(cell->face(face_no)->center()[0] == scaled_domain_length) ||
						(cell->face(face_no)->center()[1] == 0.0) ||
						(cell->face(face_no)->center()[1] == scaled_domain_height) )
						{
							cell->face(face_no)->set_boundary_id(Dirichlet);
						}
				} // end if on boundary
			} // for face_no
		} // for cell

	} 

	template<int dim> 
	void 
	Grid<dim>::
	make_Neumann_boundaries(Triangulation<dim> & triangulation)
	{
		typename Triangulation<dim>::active_cell_iterator
											cell = triangulation.begin_active(),
											endc = triangulation.end();
		// loop over all the cells
		for(; cell != endc; cell++)
		{
			// loop over all the faces of the cell and find which are on the bondary
			for(unsigned int face_no=0;
					face_no < GeometryInfo<dim>::faces_per_cell;
					face_no++)
			{
				// see if the face is on the boundary or not
				if(cell->face(face_no)->at_boundary() )
				{
					// top of the domain
					if(cell->face(face_no)->center()[1] == scaled_domain_height)
					{
						// electrolyte portion 
						if(cell->face(face_no)->center()[0] > scaled_radius_one)
						{ 
							cell->face(face_no)->set_boundary_id(Neumann);
						}
						else // in semiconductor
						{
							cell->face(face_no)->set_boundary_id(Dirichlet);
						}
					} // if top

					// bottom of the domain
					if(cell->face(face_no)->center()[1] == 0.0)	
					{
						// electrolyte portion 
						if(cell->face(face_no)->center()[0] > scaled_radius_one)
						{ 
							cell->face(face_no)->set_boundary_id(Neumann);
						}
						else // in semiconductor
						{
							cell->face(face_no)->set_boundary_id(Dirichlet);
						}
					} // if bottom 

					// left side of the domain
					if(cell->face(face_no)->center()[0] == 0.0)
							cell->face(face_no)->set_boundary_id(Neumann);
		
				} // end if on boundary
			} // for face_no
		} // for cell

	}  

	template<int dim> 
	void 
	Grid<dim>::
	make_Schottky_boundaries(Triangulation<dim> & triangulation)
	{
		typename Triangulation<dim>::active_cell_iterator
											cell = triangulation.begin_active(),
											endc = triangulation.end();
		// loop over all the cells
		for(; cell != endc; cell++)
		{
			// loop over all the faces of the cell and find which are on the boundary
			for(unsigned int face_no=0;
					face_no < GeometryInfo<dim>::faces_per_cell;
					face_no++)
			{
				if(cell->face(face_no)->at_boundary() )
				{
					//NOTE: Default is 0, which implies Interface

					//sets only the non-interface boundary to be Dirichlet	
					if(// (cell->face(face_no)->center()[1] == 0) ||
							(cell->face(face_no)->center()[1] == scaled_domain_height) )
						{
							cell->face(face_no)->set_boundary_id(Schottky);
						}
				} // end if on boundary
			} // for face_no
		} // for cell

	} 


	template <int dim>
	void 
	Grid<dim>::
	make_test_grid(Triangulation<dim> 						& triangulation,
				   const int				 				& n_global_refine)		
	{
		// make the triangulation and refine globally n_refine times
		GridGenerator::hyper_cube(triangulation,0,1);
		triangulation.refine_global(n_global_refine);
	
		// set the boundaries to be Dirichlet
		typename Triangulation<dim>::active_cell_iterator
											cell = triangulation.begin_active(),
											endc = triangulation.end();
		// loop over all the cells
		for(; cell != endc; cell++)
		{
			// loop over all the faces of the cell and find which are on the boundary
			for(unsigned int face_no=0;
					face_no < GeometryInfo<dim>::faces_per_cell;
					face_no++)
			{
				if(cell->face(face_no)->at_boundary() )
				{
					if((cell->face(face_no)->center()[1] == 0) ||
						 (cell->face(face_no)->center()[1] == 1.0) )
					{
						// set it to be Neumann boundary condition by setting the boundary 
						// indicator to be 1.  NOTE: Default is 0, which implies Dirichlet
						cell->face(face_no)->set_boundary_id(Neumann);
					}
					else
					{
						//NOTE: Default is 0, which implies Interface
						cell->face(face_no)->set_boundary_id(Dirichlet);
					}
				} // end if on boundary
			} // for face_no
		} // for cell
	} // make test_grid

	template <int dim>
	void
	Grid<dim>::
	refine_test_grid(Triangulation<dim> & triangulation,
             		 const unsigned int & local_refine)
	{
		Point<dim>	p;
		p[0]	=	0.5;
		p[1]	=	0.5;
	
    // make the triangulation and refine globally n_refine times
    for(unsigned int i =0; i <local_refine; i++)
    {
      typename Triangulation<dim>::active_cell_iterator 
										cell = triangulation.begin_active(),
                                       endc = triangulation.end();		
  		for(; cell != endc; cell++)
	 		{
				if(cell->center().distance(p) < 0.2)
					cell->set_refine_flag();
			} //  loop over all the cells
  
    triangulation.execute_coarsening_and_refinement();
	 	} // local_refine 
	} 

	template <int dim>
	void 
	Grid<dim>::
	make_test_tran_grid(Triangulation<dim> 		& triangulation,
						const int				& n_global_refine)		
	{
		// make the triangulation and refine globally n_refine times
		GridGenerator::hyper_cube(triangulation,0,1);
		triangulation.refine_global(n_global_refine);
	
		// set the boundaries to be Dirichlet
		typename Triangulation<dim>::active_cell_iterator
											cell = triangulation.begin_active(),
											endc = triangulation.end();
		// loop over all the cells
		for(; cell != endc; cell++)
		{
			// loop over all the faces of the cell and find which are on the boundary
			for(unsigned int face_no=0;
					face_no < GeometryInfo<dim>::faces_per_cell;
					face_no++)
			{
				if(cell->face(face_no)->at_boundary() )
				{
					if(cell->face(face_no)->center()[0] != 1.0)
					{
						cell->face(face_no)->set_boundary_id(Dirichlet);
					} //if x != 1
					// set to be Neumann 
					if((cell->face(face_no)->center()[1] == 0) ||
						 (cell->face(face_no)->center()[1] == 1.0) )
					{
						cell->face(face_no)->set_boundary_id(Neumann);
					}
			} // if on boundary
			} // for face_no
		} // for cell
	}

	template<int dim>
	void
	Grid<dim>::
	make_DD_Poisson_grid(Triangulation<dim>		& triangulation,
						 const	int				& n_global_refine)
	{
		// make the triangulation and refine globally n_refine times
		GridGenerator::hyper_cube(triangulation,0,1);
		triangulation.refine_global(n_global_refine);
	
		// set the boundaries to be Dirichlet
		typename Triangulation<dim>::active_cell_iterator
											cell = triangulation.begin_active(),
											endc = triangulation.end();
		// loop over all the cells
		for(; cell != endc; cell++)
		{
			// loop over all the faces of the cell and find which are on the boundary
			for(unsigned int face_no=0;
					face_no < GeometryInfo<dim>::faces_per_cell;
					face_no++)
			{
				if(cell->face(face_no)->at_boundary() )
				{
					cell->face(face_no)->set_boundary_id(Dirichlet);
				} 
			} // for face_no
		} // for cell

	}

	
	template<int dim>
	void
	Grid<dim>::
	output_mesh(Triangulation<dim> & triangulation,
				const std::string  & grid_name)
	{
		std::ofstream out(grid_name.c_str());
//		std::ofstream out("grid.eps");
		GridOut grid_out;
		grid_out.write_msh(triangulation,out);
		out.close();
	
	}

	template<int dim> 
	void 
	Grid<dim>::
	print_grid(Triangulation<dim> & triangulation,
			   const std::string  & grid_name)
	{	
		std::ofstream out(grid_name.c_str());
//		std::ofstream out("grid.eps");
		GridOut grid_out;
		grid_out.write_eps(triangulation,out);
		out.close();
	}



} 
