#ifndef _ASSEMBLY_H__
#define _ASSEMBLY_H__

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h> // Lagrange dg fe elements
//#include <deal.II/fe/fe_dgp.h> // Legendre dg fe elements
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>



#include <fstream>
#include <iostream>

/// \namespace Assembly Contains the local assembly objects and local data objects. 
/** Assembly holds the local assembly objects and local data objects. 
*	It is necessary to run the code in parallel with the 
*	<a href="https://www.threadingbuildingblocks.org/">Thread Building Blocks</a>.
*	See
* <a href="https://www.dealii.org/developer/doxygen/deal.II/group__threads.html">
* deal.ii's explanation on parallel computing with shared memory</a> to see why this 
* necessary and how it works. */
namespace Assembly
{

	using namespace dealii;

	/// \brief Temporary assembly object that holds FE's and functions.
	template<int dim>
	struct AssemblyScratch
	{
		/** This is the constructor you would call if you wanted to create
		 * 	 an AssemblyScratch object explicitly and loop over all the cells and
		 * 	 assemble the local matrices and vectors by hand.  This would be done
		 *   sequentially.
		*/ 
		AssemblyScratch(const FiniteElement<dim> & Poisson_fe,
				const FiniteElement<dim> & carrier_fe,
				const Quadrature<dim>	& quadrature,
				const Quadrature<dim-1>	& face_quadrature);

		/** This is the copy constructor that will be called in the WorkStream
	 	* 	 call that will build and distribute the local matrices and vectors
		*   in parallel using Thread Building Blocks.
		*/ 
		AssemblyScratch(const AssemblyScratch & scratch);

		/** Finite element for Poisson's equation in a cell. */	
		FEValues<dim>				Poisson_fe_values;

		/** Finite element for Poisson's equation in a face of a cell. */	
		FEFaceValues<dim>			Poisson_fe_face_values;

		/** Finite element for carrier in a cell. */	
		FEValues<dim>				carrier_fe_values;

		/** Finite element for carrier on a face of a cell. */	
		FEFaceValues<dim>			carrier_fe_face_values;

		/** Finite element for carrier on a subface of a cell. */	
		FESubfaceValues<dim>			carrier_fe_subface_values;

		/** Finite element for carrier on a neighbor cells face. */	
		FEFaceValues<dim>			carrier_fe_neighbor_face_values;

		/// Vector which holds right hand values of Poissons equation
		/// at the quadrature points of the local cell
		std::vector<double>			Poisson_rhs_values;

		/// Vector which holds boundary conditions
		/// values of Poissons equation at the quadrature points of the local cell
		std::vector<double>			Poisson_bc_values;

		/// Vector which holds built_in potential 
		/// values of Poissons equation at the quadrature points of the local cell
		std::vector<double>			Poisson_bi_values;

		/// Vector which holds the electric field values
		///	at the quadrature points of the local cell
		std::vector<Tensor<1,dim>>  		electric_field_values;

		/// Vector which holds the carrier_1 density values
		///	at the quadrature points of the local cell
		std::vector<double> 			old_carrier_1_density_values;
	
		/// Vector which holds the carrier_2 density values
		///	at the quadrature points of the local cell
		std::vector<double>			old_carrier_2_density_values;

		/// Vector which holds the carrier_1 current density values
		/// at the quadrature points of the local cell
		std::vector<Tensor<1,dim>>  		carrier_1_current_values;
	
		/// Vector which holds the carrier_2 current density values
		/// at the quadrature points of the local cell
		std::vector<Tensor<1,dim>>  		carrier_2_current_values;

		/// Vector which holds the donor doping profile density values
		///	at the quadrature points of the local cell
		std::vector<double>			donor_doping_values;	

		/// Vector which holds the acceptor doping profile density values
		///	at the quadrature points of the local cell
		std::vector<double>			acceptor_doping_values;	
	
		/// Vector which holds the generation functions values
		///	at the quadrature points of the local cell
		std::vector<double>			generation_values;

		/// Vector which holds the carrier_1 dirichlet boundary conditions
		///	at the quadrature points of the local face
		std::vector<double>			carrier_1_bc_values;

		/// Vector which holds the carrier_2 dirichlet boundary conditions
		///	at the quadrature points of the local face
		std::vector<double>			carrier_2_bc_values;


		/// Vector which holds electron interface density values
		///	at the quadrature points of the local face
		std::vector<double>			electron_interface_values;

		/// Vector which holds hole interface density values
		///	at the quadrature points of the local face
		std::vector<double>			hole_interface_values;

		/// Vector which holds reductant interface density values
		///	at the quadrature points of the local face
		std::vector<double>			reductant_interface_values;

		/// Vector which holds oxidant interface density values
		///	at the quadrature points of the local face
		std::vector<double>			oxidant_interface_values;

		/// Vector which holds the potential basis function values at
		/// the quadrature points of the local cell
		std::vector<double>			psi_i_potential;
		
		/// Vector which holds the electric field basis function values at
		/// the quadrature points of the local cell.
		std::vector<Tensor<1,dim>>		psi_i_field;
		
		/// Vector which holds the charge carrier density basis function values at
		/// the quadrature points of the local cell.
		std::vector<double>			psi_i_density;
		
		/// Vector which holds the charge carrier current basis function values at
		/// the quadrature points of the local cell.
		std::vector<Tensor<1,dim>>		psi_i_current;

		/// Vector which holds the normal vector values at
		/// the quadrature points of the local face.
		std::vector<Tensor<1,dim>>		normals;

	};

	
	/// \namespace Poisson Contains all of the local CopyData for Poisson's equation.
	namespace Poisson
	{	

		/** \brief Local data for Poisson's equation */
		/** The object to store local Poisson data before being copied to the global
		*	data structures. */
		template<int dim>
		struct CopyData
		{
			/** This is the constructor you will call explicitly when looping over
			 * the cells and building the local matrices and vectors by hand sequentially.
			*/
			CopyData(const FiniteElement<dim> & Poisson_fe);
		
			/** This is the copy constructor that will be called in the WorkStream
			 * 	 call that will store the local matrices and vectors
			 *   in parallel using Thread Building Blocks.
			*/ 
			CopyData(const CopyData & data);
		
			/// Local vector for this cells contribution to Poissons equation's
			/// right hand side
			Vector<double>				local_rhs;

			/// Local matrix for this cells contribution to Poissons equation's
			/// left hand side matrix
			FullMatrix<double>			local_matrix;

			/// Vector which holds the local to global degrees of freedom info
			///  of this cell
			std::vector<types::global_dof_index>	local_dof_indices;

		}; // end struct

	} // end namespace Poisson

/// \namespace DriftDiffusion Contains all of the local CopyData for drift-diffusion equation.
	namespace DriftDiffusion
	{
		/** \brief Local data for drift-diffusion equation. */
		/** The object to store local drift-diffusion equation data 
		*  before being distributed to global data */
		template<int dim>
		struct CopyData
		{
			/** This is the constructor you will call explicitly when looping over
			 * the cells and building the local matrices and vectors by hand sequentially.
			*/
			CopyData(const FiniteElement<dim> & carrier_fe);

			/** This is the copy constructor that will be called in the WorkStream
			 * 	 call that will store the local matrices and vectors
			 *   in parallel using Thread Building Blocks.
			*/ 
			CopyData(const CopyData & data);

			/// Local vector for this cells contribution to carrier_1 transport
			/// equation's right hand side: drift, generation recombination and
			/// boundary terms
			Vector<double>						local_carrier_1_rhs;

			/// Local vector for this cells contribution to carrier_2 transport
			/// equation's right hand side: drift, generation recombination and
			/// boundary terms
			Vector<double>						local_carrier_2_rhs;

			/// Local matrix for this cells contribution to carrier_1 transport
			/// equation's left hand side matrix which involves body intergrals
			FullMatrix<double>					local_matrix_1;

			/// Local matrix for this cells contribution to carrier_2 transport 
			///	equation's left hand side matrix which involves body intergrals
			FullMatrix<double>					local_matrix_2;

			/// Local matrix for this cells contribution to drift diffusion equation's
			/// mass matrix
			FullMatrix<double>					local_mass_matrix;

			// flux matrices
			/// Matrix for flux from both test function and trial function interior
			/// to this cell
			FullMatrix<double>					vi_ui_matrix;
			/// Matrix for flux from test function interior
			/// to this cell and trial function exterior to this cell

			FullMatrix<double>					vi_ue_matrix;
			/// Matrix for flux from test function exterior
			/// to this cell and trial function interior to this cell

			FullMatrix<double>					ve_ui_matrix;

			/// Matrix for flux from both test function and trial function exterior
			/// to this cell
			FullMatrix<double>					ve_ue_matrix;

			/// Vector which holds the local to global degrees of freedom info
			///  of this cell
			std::vector<types::global_dof_index>	local_dof_indices;
			/// Vector which holds the local to global degrees of freedom info
			/// of this cell's neighbor cell
			std::vector<types::global_dof_index>	local_neighbor_dof_indices;

		}; // end struct

	} // end DriftDiffusion 

} // end namespace Assembly
#endif
