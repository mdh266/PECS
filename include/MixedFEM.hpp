#ifndef _MIXED_FEM_H__
#define _MIXED_FEM_H__

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/convergence_table.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h> // for block structuring
#include <deal.II/lac/sparse_matrix.h> // for block structuring
#include <deal.II/lac/block_vector.h> // for block structuring
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h> // dg fe elements
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <functional>

// Raviart-Thomas elements
#include <deal.II/fe/fe_raviart_thomas.h>

// Tensor valued functions
#include <deal.II/base/tensor_function.h>

// Multithreading
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/thread_management.h>

#include "Assembly.hpp"
#include "BiasValues.hpp"
#include "Grid.hpp"
#include "../source/test_functions.cpp"
#include "Poisson.hpp"
#include "Parameters.hpp"
#include "PostProcessor.hpp"

/// \namespace MixedPoisson Holds MixedFEM, the class for a mixed FEM on Poissons equation.
namespace MixedPoisson
{
	using namespace dealii;

	/** This class provides a mixed finite element method solver for Poisson's equation
	*   in the drit-diffusion-Poisson system with constant
	*		Debeye length. The class also can compute the error of the approximation
	* 	for testing. The data is stored in Poisson::PoissonData.  
	* 
	*   The general Poisson equation problem we are solving is, 
	* 	\f[ \begin{align}
	*  \epsilon_{r}^{-1} \ \textbf{D} \ + \ \nabla \Phi \ &= \ 0  && \text{in} \; \Omega \\
	*  \ \nabla \ \cdot \ \textbf{D} \  &= \frac{1}{\lambda^{2}}
	*  f(\textbf{x}) 
	*						&& \text{in} \; \Omega \\
	*		\textbf{D} \ \cdot \ \boldsymbol \eta \ &= \ 0 
	*						&& \text{on} \; \partial \Omega_{N} \\
	*		\Phi \ &=  \ \Phi_{D}
	*							 && \text{on} \; \partial \Omega_{D} 
	*   \end{align} \f]
	*
	* For \f$\lambda^{2} \ = \ \frac{\Phi^{*} \epsilon }{L^{2} q C^{*}}\f$. 	
	*	
	*	This becomes the problem in the weak formulation: 
	* 
	* Find 
	*	\f$( \ \Phi \ , \ \textbf{D} \ ) \ \in 
	*	\left( \ \text{W} \ \times \ [ \ 0 , \ T \ ] , 
	*	\ \textbf{V}^{d} \ \times \
	*  	[ \ 0 , \ T \ ] \ \right) \f$ such that: 
	*
	* 	\f[ \begin{align}
	* \ \left(  \textbf{p} \ , \ \epsilon_{r}^{-1} \ \textbf{D} \right)_{\Omega} 
	*	\  - \
	* \left( \ \boldsymbol \nabla \cdot \ \textbf{p}  \ , \  \Phi  \right)_{\Omega}  
	*  \; &= \;  
	*  - \langle \textbf{p} \ , \ \Phi_{D} \rangle_{\Gamma_{D}}   \\
	* -  \left(  v \ , \ \boldsymbol \nabla  \cdot \textbf{D} \right)_{\Omega} \;  
	*  &= \; -  
	*	\frac{1}{\lambda^{2}} \
	* \left(  v,\ f(\textbf{x})  \right)_{\Omega} 
	*	\end{align} \f]
	* 
	*
	*  For all \f$( v \  , \ \textbf{p}  ) \, \in \, W \, \times\, \textbf{V}^{d}\f$. 
	*	
	* 
	*	We obtain the electric field \f$-\boldsymbol \nabla \Phi\f$ by setting,
	*  	
	*	\f[ 
	* - \boldsymbol \nabla \Phi \ = \ \epsilon_{r}^{-1}  \ \textbf{D}
	*	\f]
	* 
	*  
	*/	

	template<int dim>
	class MixedFEM
	{
		public:
			/** \brief Simple constructor which instantiates test function. */
			MixedFEM();
	
			/** \brief Simple desctructor which clears the dof_handler. */
			~MixedFEM();

			/** \brief Assembles the matrix for local cell in the mixed 
			* finite element metho in the mixed 
			* finite element method  for Poissons equation. */
			/** This function can either be called when looping throug the cells by 
			* 	hand and assemblying the global matrix sequentially or by using
			* 	the WorkStream to assemble it in parallel.  If you use it in squential
			* 	mode, than you must have the AssemblyScratch and Poisson::CopyData
			*		instantiated before calling this function.
			* 
			*		This function loops through the quadrature points of this cell and
			* 	assembles a local matrix corresponding to the Mixed Method applied 
			* 	Poisson equation and stores it in Poisson::CopyData.
			*		
			*	 @param debeye length \f$\lambda{^2}\f$ is a constant for now.
			*	 @param scratch is the temporary scratch objects and data structures
			* 			that do the work.
			*	 @param data is the local data structures to this cell that the computed
			*			results get stored in before being passed off onto the global 
			*			data structures.
			*
			* 	The matrix will be of the form,
			*   \f[ \left[ \begin{matrix}
			*			 A & B \\
			*			 B^{T} & 0 
			*			 \end{matrix} \right]  \f] 
			*
			* where, 
			* 
			*	 \f[ A(\textbf{p},\textbf{D}) \; = \; 
			*				\int_{\Omega} \  \textbf{p} \  \cdot \  \textbf{D} \ dx \f]
			*
			*	 \f[ B(\textbf{p},\Phi) \; = \; 
			*				\int_{\Omega} \ \nabla \ \cdot \ \textbf{p} \  \Phi  \ dx \f]
			*
			*	 \f[ B^{T}(v,\textbf{D} ) \; = \; 
			*				\int_{\Omega} \ \nabla v \ \cdot \ \textbf{D} \ dx \f]
			*
			*
			*/
			void		
			assemble_local_Poisson_matrix(	
				const typename DoFHandler<dim>::active_cell_iterator & cell,
				Assembly::AssemblyScratch<dim>			     & scratch,
				Assembly::Poisson::CopyData<dim>		     & data,
				const double 	   				     & semi_permittivity,
				const double 		  			     & elec_permittivity,
				const double 	 				     & scaled_debye_length);
	
			/** Assembles the local cell's right hand side for the test problem. */
			void
			assemble_local_test_rhs(
				const typename DoFHandler<dim>::active_cell_iterator 	& cell,
				Assembly::AssemblyScratch<dim>				& scratch,
				Assembly::Poisson::CopyData<dim>			& data);

			/** \brief Computes the local error of your approximation for the 
				*  mixed method on the cell and stores the errors in
				* <code>potential_error<code/> and <code>field_error<code/>.
				*
				* @param solution is the <code>Poisson_solution<code/>
				*  vector to the Mixed Method.
				* @param potential_error is \f$ L^{2} \f$  error of the approximation 
				*					the potential.
				* @param field_error is \f$ L^{2} \f$ error of the approximation 
				*					electric field.
				*/
			void 
			compute_errors(const Triangulation<dim>	& triangulation,
					DoFHandler<dim>		& Poisson_dof_handler,
					Vector<double>		& solution,
					double 				& potential_error,
					double 				& field_error) const;

			/** Prints the electric field and potential with units put back in.*/
			/** \note:  I think this will only work when the code is compiled in
 			*						release mode. */	
			void	output_rescaled_results(DoFHandler<dim>	 & dof_handler,
						const Vector<double> 				 & solution,
						const ParameterSpace::Parameters 	 & sim_params,
						const unsigned int time_step_number) const;
		
			/** Prints the electric field and potential without the units. */	
			void	output_unscaled_results(DoFHandler<dim>		& dof_handler,
						const Vector<double> 					& solution,
						const unsigned int 						time_step_number) const;
	
		private:
			enum
			{
				/// Marker for cells's face to be on the Interface.	
				Interface,
				/// Marker for cells's face to be on the Dirichlet boundary.
				Dirichlet,
				/// Marker for cells's face to be on the Neumann boundary.
				Neumann,
				/// Markert for cells' face to be on the Schottky boundary.
				Schottky
			};


			enum
			{
				// Material id for cell to be bulk semiconductor cell.
				semiconductor_id,
				// Material id for cell to be boundary layer semiconductor cell.
				semi_boundary_layer_id,
				// Material id for cell to be bulk electrolyte cell.
				electrolyte_id,
				// Material id for cell to be boundary layer electrolyte cell.
				elec_boundary_layer_id
			};

		public:
			/** Uncoupled test problem solution. 
			 * 	See SOLARCELL::SolarCellProblem::test_steady_state*/
			const test_Poisson::TrueSolution<dim>			 test_Poisson_solution;
			/** Uncoupled test problem Dirichlet BC function.
			 * 	See SOLARCELL::SolarCellProblem::test_steady_state*/
			const test_Poisson::DirichletBoundaryValues<dim> test_Poisson_bc;
			/** Uncoupled test problem RHS function.
			 * 	See SOLARCELL::SolarCellProblem::test_steady_state*/
			const test_Poisson::RightHandSide<dim>			 test_Poisson_rhs;

			/** Coupled test problem RHS function.
			 * 	See SOLARCELL::SolarCellProblem::test_DD_Poisson.*/
			test_DD_Poisson::Poisson_RightHandSide<dim>		test_coupling_Poisson_rhs;	

	};

} 

#endif
