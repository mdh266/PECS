#ifndef _CARRIER_PAIR_H__
#define _CARRIER_PAIR_H__

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <string>
#include <fstream>
#include "Carrier.hpp"
#include "Parameters.hpp"

/** \namespace ChargeCarrierSpace Namespace for electron-hole or redox pairs.*/
namespace ChargeCarrierSpace
{
	using namespace dealii;

	/** \brief CarrierPair is the common data between electron-hole or redox pairs.
 	 *
 	 * It contains information like the name of the pair (this is set in the 
 	 * constructor of SOLARCELL::SolarCellProblem) as well as sparsity patterns, 
 	 * information on the mesh, mass matrix for the pair and the electron-hole or
 	 * reductant-oxidant pairs themselves.	*/
	template<int dim>
	class CarrierPair
	{
		public:
			CarrierPair();

			~CarrierPair();

			/** \brief This sets the name of the pair, it is called in constructor of 
 			* SOLARCELL::SolarCellProblem. */
			void 
			set_name(const std::string & str_name);

			/** \brief Distributes the dofs for the carrier pair. 
 			* Both carriers have the same distribution of dofs.*/
			void 
			setup_dofs(const FESystem<dim>	 & fe,
				   DoFHandler<dim>	 & dof_handler);

			/** \brief Prints the number DOFS in each of the carriers. */
			void
			print_info(DoFHandler<dim>	 & dof_handler);

			/** \brief Prints the DOFS into a files to be read in. */
			void 
			print_dofs();

			/** \brief Reads the DOFS from output of last run to be init cond for this one.*/
 			void
			read_dofs();

			/** \brief Set \f$\tilde{\tau} = 1\f$ and mobilities \f$=1\f$. */
			void 
			set_semiconductor_for_testing(
					const ParameterSpace::Parameters & sim_params);

			/** \brief Set \f$\tilde{\tau} = 1\f$ and mobilities \f$=1\f$. */
			void 
			set_electrolyte_for_testing(
					const ParameterSpace::Parameters & sim_params);
	
			/** Name of the material: semiconductor or electrolyte */
			std::string			material_name;

			/** Sparsity pattern for system matrix to be inverted 
			* foreach carrier in this pair. */
			SparsityPattern			system_sparsity_pattern;

			/** Sparsity pattern for the mass matrix for each carrier in
			*  this pair. */
			SparsityPattern			mass_sparsity_pattern;
		
			/** Constraint matrix for this pair.  This is mostly used 
			* to distributed local dofs to global dofs during construction,
			* since for DG methods there is not much to constrain. */
			ConstraintMatrix		constraints;

			/** Mass matrix for the carrier pair.  
			*  Note we only need one, this implies
 			* they use the same time step delta_t. */
			SparseMatrix<double>		mass_matrix;

			/** Carrier 1 should be electron or reductant.*/
			Carrier<dim>			carrier_1;

			/** Carrier 2 should be hole or oxidant.*/
			Carrier<dim>			carrier_2;

			/** The penalty term for fluxes in LDG.  This is the same for
			*  each carrier in a pair since they have the same mesh. */
			double 				penalty;

			/** Time step for each of the carriers.  This is the same.  
			* Any time you want to change this you need to update the mass 
			* matrix and system matrix.*/
			double				delta_t;

			/** The dieletric constant for this material.*/
			double				material_permittivity;

	};

}
#endif
