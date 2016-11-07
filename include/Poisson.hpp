#ifndef _Poisson_H__
#define _Poisson_H__

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

namespace Poisson
{

	using namespace dealii;
	/** \brief Data structures, functions/paramaters, and solver for Poisson
 	*  equation.*/

 /** This object holds the matrix, vectors, constrains and solver for the
 *  mixed finite element approximation to Poisson's equation.  See
 *  MixedPoisson::MixedFEM for more information on the mixed method.
 */
	template<int dim>
	class PoissonData
	{
		public:
			PoissonData();
			
			~PoissonData();
			/** \brief Distributes the dofs for the Mixed finte element method. */
			void setup_dofs(const FESystem<dim>	& fe,
					DoFHandler<dim>		& dof_handler);

			/// \brief Prints out the number of dofs for MFEM.*/
			void	print_info(DoFHandler<dim>	& dof_handler);
			/// set the direct solver.
			void set_solver();

			/// Solves the linear system and distributes the solution with constraints.
			void solve();
	
			/// This matrix holds the constraints for hanging nodes and Neumann BC.
			ConstraintMatrix			constraints;
	
			/// Sparsity pattern for the mixed method.
			SparsityPattern				sparsity_pattern;

			/** 	The matrix will be of the form,
			*   \f[ \left[ \begin{matrix}
			*			 A & B \\
			*			 B^{T} & 0 
			*			 \end{matrix} \right]  \f]
			*  See MixedPoisson::MixedFEM for more information on the mixed method.
 			*/
			SparseMatrix<double>			system_matrix;
	
			/** This will be used to store the right hand side of the mixed method,
 			* \f[
 			* = \ 	\frac{1}{\lambda^{2}} \
			* \left(  v,\ \rho_{n}^{e} \ - \ \rho_{p}^{e} \
			* -\ \left( \rho_{n} \ - \ \rho_{p} \right) \right) 
			* \f] 
			*  See MixedPoisson::MixedFEM for more information on the mixed method.
 			*/
			Vector<double>				system_rhs;

			/// solution to the linear system.
			Vector<double>				solution;

			/// Direct solver
			SparseDirectUMFPACK			solver;

		private:
			enum
			{
				Interface,
				Dirichlet,
				Neumann
			};

	};
}

#endif
