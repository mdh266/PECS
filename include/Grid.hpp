#ifndef _GRID_H__
#define _GRID_H__
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>
//#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <iostream>
#include <fstream>
#include <string>
#include "Parameters.hpp"

/// \namespace Grid_Maker for Grid class which creates all different types of meshes.
namespace Grid_Maker
{
	using namespace dealii;

	/// \brief This object will be used to build meshes over triangulations provided.
	/** This object will be used to build meshes over triangulations provided.  Will
	* 	build meshes for the semiconductor-electrolyte interface problem as well 	* as meshes for testing.  Can also be used to print meshes to .eps form.
	*/
	template<int dim>	
	class Grid
	{
		public:
			/** Grid object constructor just initializes private data to 
			* 	which portions of the boundary are Dirichlet and Neumann */
			Grid(const ParameterSpace::Parameters & sim_parms);

			/** Makes the outline grids for the bulk semiconductor + boundary layer
 			*	and electrolyte + boundary layer.  It then merges them to make Poisson
 			*	grid.  All the grids are globally refined.  Then we loop over all the
 			*	cells in the triangulaitions and locally refine the cells in the 
 			*	boundary layer.
 			*/
			void make_grids(Triangulation<dim> & semiconductor_triang,
					Triangulation<dim> & electrolyte_triang,
					Triangulation<dim> & Poisson_triang,
					const bool	   & full_system);

			/** Makes a grid which is of the form, \image html semiconductor-grid.png
			*	The top radius is <code> radius one <code/> while the top radius is 
			*	<code> radius two <code/>. It first makes a temporary triangulation
			* refines it, and then prints it out as "temp.msh", then reads it
			* it in again and attaches it to the passed triangulation.  The reason
			* for doing this is so that I can merge this triangulation with another
			* since in order to do that, they must be on the most refined level.
			*/
			void make_semiconductor_grid(Triangulation<dim>  & triangulation);

			/** Makes a grid which is of the form, \image html electrolyte-grid.png
			*	The top radius is <code> radius one <code/> while the top radius is 
			*	<code> radius two <code/>.  It first makes a temporary triangulation
			* refines it, and then prints it out as "temp.msh", then reads it
			* it in again and attaches it to the passed triangulation.  The reason
			* for doing this is so that I can merge this triangulation with another
			* since in order to do that, they must be on the most refined level.
			*/
			void make_electrolyte_grid(Triangulation<dim> 	 & triangulation);
			/** Makes a grid which is of the form, \image html Poisson-grid.png 
			*	The top radius is <code> radius one <code/> while the top radius is 
			*	<code> radius two <code/>.   It does so by merging the two provided grids.*/

			/**	@param semiconductor_triang which looks like 
			*	\image html semiconductor-grid.png */ 

			/** @param electrolyte_triang which looks like
 				* \image html electrolyte-grid.png. */	
			void make_merged_grid(const Triangulation<dim>	& semiconductor_triang,
					      const Triangulation<dim>	& electrolyte_triang,
					      Triangulation<dim>	& merged_triangulation);			

			///\brief Subroutine that tags boundaries to be Dirichlet portions.
			/** This subroutine loops over all the cells in <code>triangulation<code/h>
			* and finds which subroutines are on the boundary.  It then
			* flags these faces on the boundary to be <code>Dirichlet<code/>
			* portions of boundary.  
	 		*/
			void make_Dirichlet_boundaries(Triangulation<dim> & triangulation);

			/// \brief Tags boundaries of the semiconductor as Neumann/Insulating portions.
			/** This subroutine loops over all the cells in <code>triangulation<code/>
			* and finds which subroutines are on the boundary.  It then
			* flags these faces to be <code>Neumann<code/> portions of boundary.
			* The choice of which boundary faces are <code>Neumann<code/> is
			* preset in this subroutines source code. 
			*
			* NOTE: BECAREFUL TO MATCH WITH CUBE END POINTS AND DIMENSION.
			* NOTE: <code>Neumann<code/> corresponds to a Neumann bounary condition
			* for Poisson and a Robin boundary condition for the drift-diffusion equation.
			*/
			void make_Neumann_boundaries(Triangulation<dim> & triangulation);

			/// \brief Tags boundaries of the semiconductor as Schottky portions.
			void make_Schottky_boundaries(Triangulation<dim> & triangulation);
	
	 		/** \brief Creates a simple cubic grid with mixed boundaries. */
			/** Takes <code>triangulation<code/> object and creates a mesh on it.
			* Globally refines the mesh using the refinement level from  
			* <code>params.n_global_refine</code>.
			* 
 			* See	SOLARCELL::SolarCellProblem<dim>::test_steady_state. 
			*/
			void make_test_grid(Triangulation<dim>  & triangulation,
					    const int		& n_global_refine);

			/** Locally refines the triangulation around the point
			* x = 1/2, y = 1/2.*/
			void refine_test_grid(Triangulation<dim> & triangulation,
           				      const unsigned int & local_refine);

			/** \brief Grid for testing the transient LDG problems.*/
			/** See 
 			* 	SOLARCELL::SolarCellProblem<dim>::test_transient and 
 			* 	SOLARCELL::SolarCellProblem<dim>::test_interface_coupling . 
 			*/
			void
			make_test_tran_grid(Triangulation<dim> 	& triangulation,
					    const int		& n_global_refine);
		
			/** Makes a cube and sets all the boundaries to be Dirichlet. */
			void
			make_DD_Poisson_grid(Triangulation<dim>	& triangulation,
					     const	int	& n_global_refine);
	
			
			/// \brief Subroutine that prints the grid into a .msh file
			void output_mesh(Triangulation<dim> & triangulation,
					const std::string  & grid_name);


			/// \brief Subroutine that prints the grid into a .eps file
			void print_grid(Triangulation<dim> & triangulation,
					const std::string  & grid_name);

		private:
			enum
			{
				Interface,
				Dirichlet,
				Neumann,
				Schottky
			};


			enum
			{
				semiconductor_id,
				semi_boundary_layer_id,
				electrolyte_id,
				elec_boundary_layer_id
			};

			/// Set in constructor.			
			double scaled_domain_height;
			/// Set in constructor.			
			double scaled_domain_length;
			/// Set in constructor.			
			double scaled_radius_one;
			/// Set in constructor.			
			double scaled_radius_two;
			/// Set in constructor.			
			double scaled_boundary_layer;
			/// Set in constructor.			
			unsigned int n_global_refine;
			/// Set in constructor.			
			unsigned int n_local_refine;
			/// Set in constructor.			
			bool insulated;
			/// Set in constructor.			
			bool use_boundary_layer;
			/// Set in constructor.			
			bool schottky;
	
	};

}

#endif 
