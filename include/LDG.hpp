#ifndef _LDG_H__
#define _LDG_H__

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/convergence_table.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h> // Lagrange dg fe elements
//#include <deal.II/fe/fe_dgp.h> // Legendre dg fe elements
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>

#include <fstream>
#include <iostream>

// multithreading
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/thread_management.h>

#include "Grid.hpp"
#include "Carrier.hpp"
#include "CarrierPair.hpp"
#include "Parameters.hpp"
#include "Assembly.hpp"
#include "PostProcessor.hpp"
#include "../source/test_functions.cpp"


/// \namespace LDG_System Basically just the namespace for LDG_System::LDG.  
namespace LDG_System
{
	using namespace dealii;
	/** \brief This class builds two system matrices using the LDG method for the
 	* two ChargeCarrierSpace::Carrier objects in each of the 
 	* ChargeCarrierSpace::CarrierPair objects. */
	/** This class will build portions of the LDG local cell matrix for the general
		* (non-dimensional) drift diffusion equation:
		*	
		*	\f[ \begin{align}
		* 	u_{t} \  - \  \boldsymbol \nabla \ \cdot  
		* 	\ \mu \left( s \boldsymbol \nabla \Phi u \ + \ \boldsymbol \nabla u \ \right) 
		* 	\; &= \; 
		* 	R(u) + G	&& \text{in} \;  \Omega   \\
 		* 	u \; &= \; u_{D} &&  \text{on} \;  \Omega_{D}     \\
		* 	- \mu \left(s \boldsymbol \nabla \Phi \ u \ + \ \boldsymbol \nabla u \ \right) 
		* 	\ \cdot \ \boldsymbol \eta 
		*	\;  &= \; K (u) && \text{on} \; \partial \Omega_{N}
		* 	\end{align} \f]
		*
		* 	We rewrite this in mixed form:
		*
		* \f[ \begin{align} 
		*	u_{t} \ + \ \nabla \ \textbf{q} \ &= \ R(u) \ + G && \text{in} \ \Omega \\
		*			\mu^{-1} \ \textbf{q} \ & =
		*					 \ -s \nabla \Phi \ u \ - \nabla u && \text{in} \ \Omega \\
		*			\mu^{-1} \ \textbf{q} \ \cdot \boldsymbol \eta &= 
		*							\ K(u) && \text{on} \ \partial \ \Omega_{N} \\
		*			u \ &= \ u_{D} && \text{on} \ \partial \Omega_{D} 
		*	\end{align} \f]
		* 	
		* 	The weak formulation for IMEX will be:
		*
		*	Find \f$(u, \textbf{q}) \in W \times [t^{k-1}, t^{k}] 
		*			\times \textbf{W}^{d} \times[ t^{k-1}, t^{k}] \f$ such that,
		*
		*  \f[ \begin{align}
		*	\frac{1}{\Delta t} \left(  v , u^{k}  \right)
		* 	-
 		*	\tau \langle  [[ \ v \ ]] , 
		*				 [[ u^{k} ]] \rangle_{\mathcal{E}_{h}^{0}}  
 		*	-
		*	\left( \boldsymbol \nabla v  ,  \textbf{q}^{k}  \right)
 		*	+
		*	 \langle [[ \ v \ ]] , 
		*	\{ \textbf{q}^{k} \} \rangle_{\mathcal{E}_{h}^{0} \cap \partial \Omega_{D}} 
		*	\ &= \
		*	\left( v ,  R(u^{k-1})  + G \right) - 
		*	\langle   v, K( u^{k-1})    \rangle_{\Sigma}   \nonumber \\
		*	 - \left(  \boldsymbol \nabla \cdot \textbf{p} ,   u^{k-1} \right)
		*  \ - \
		*	\langle  [[ \,  \textbf{p} \, ]] ,  
		*	\{  u^{k} \}  \rangle_{\mathcal{E}^{0}_{h} \cap \partial \Omega_{N}}  
		*	\ + \
		*	 \left( \textbf{p} , \textbf{q}^{k} \right)
		*	\ &= \
		*	 +
		*	\left( s \textbf{P}  \cdot \boldsymbol \nabla \Phi , u^{k-1} \right)
		*	- 
		*	\langle  \textbf{p}   ,  u_{D}  \rangle_{ \partial \Omega_{N} }
		*	\end{align} \f]
		*
		*  For all \f$(v,\textbf{p}) \in W \times \textbf{W}^{d})\f$.
		*
		* 	The corresponding matrix will be of the form,
		*   \f[ \left[ \begin{matrix}
		*		\mu^{-1} A & B_{1} + F_{1} \\
		*		 B_{2} + F_{2} & \frac{1}{\Delta t} M + C
		*	 \end{matrix} \right]  \f] 
		* 	This matrix will assembled once and stored in each carriers 
		*		ChargeCarrierSpace::Carrier::system_matrix.  
		*  	The corresponding right hand side vector will assembled at every time step 
		* 	and stored in ChargeCarrierSpace::Carrier::system_rhs.
		*
		* 	NOTE: We use IMEX time stepping so all non-linear terms and drift terms are 
		*				time lagged and therefore on the right hand side. While they are
		*			are built in parallel, this place takes place outside of this class and
		*				in SOLARCELL::SolarCellProblem::assemble_semiconductor_rhs and 
		*				SOLARCELL::SolarCellProblem::assemble_electrolyte_rhs.
		*
		*	NOTE: The LDG flux matrices \f$F_{1}\f$ and \f$F_{2}\f$.
		*				are built sequentially and this occurs in a loop
		*  			outside this class, but calls the function assemble_local_flux_terms 
		*				or in the case of a locally/adaptively refined mesh calls
		*				assemble_local_child_flux_terms.
		*
		*	NOTE:  This assembles the ChargeCarrierSpace::Carrier::system_matrix for both	
		*			ChargeCarrierSpace::CarrierPair::carrier_1 and 
		*			ChargeCarrierSpace::CarrierPair::carrier_2 
		*			at the same time.  It therefore deals with two equations of the 
		*			above form, but proceed to explain it on only one of the equations.
		*/		
	template<int dim>
	class LDG
	{
		public:
			/** \brief A simple constructor which instantiates test functions. */
			LDG();
			
			/** \brief A simple destructor which clears carrier_dof_handler. */
			~LDG();	
	
			/* \brief Assembles the local mass matrix for this cell. */
			/** This function can either be called when looping throug the cells by 
			* 	hand and assemblying the global matrix sequentially or by using
			* 	the WorkStream to assemble it in parallel.  If you use it in squential
			* 	mode, than you must have the Assembly::AssemblyScratch and 
			*	Assembly::DriftDiffusion::CopyData
			*		instantiated before calling this function.
			* 
			*		This function loops through the quadrature points of this cell and
			* 	assembles a local mass matrix for the LDG method applied to the 
			* 	drift-diffusion equation and stores it in Assembly::DriftDiffusion::CopyData.
			*		
			*	 @param delta_t is the fixed time step size.
			*  @param scaled_mobility is the scaled mobility constant \f$\mu\f$.
			*	 
			* 	The matrix will be of the form,
			*   \f[ \left[ \begin{matrix}
			*			 0 & 0 \\
			*			 0 & \frac{1}{\Delta t} M
			*			 \end{matrix} \right]  \f] 
			*
			* where, 
			* 
			*
			*	 \f[ M(v,u ) \; = \; 
			*				\int_{\Omega} \ v \ u \ dx \f]
			*/
			void		
			assemble_local_LDG_mass_matrix(
					const typename DoFHandler<dim>::active_cell_iterator & cell,
					Assembly::AssemblyScratch<dim>			& scratch,
					Assembly::DriftDiffusion::CopyData<dim>	& data,
					const double			 				& delta_t);
							
			/** \brief Assembles the local sytem matrix for this cell. */
			/** This function can either be called when looping throug the cells by 
			* 	hand and assemblying the global matrix sequentially or by using
			* 	the WorkStream to assemble it in parallel.  If you use it in squential
			* 	mode, than you must have the AssemblyScratch and DriftDiffusion::CopyData
			*		instantiated before calling this function.
			* 
			*		This function loops through the quadrature points of this cell and
			* 	assembles a local mass matrix for the LDG method applied to the 
			* 	drift-diffusion equation and stores it in DriftDiffusionCopyData.
			*		
			*	 @param delta_t is the fixed time step size.
			*  @param scaled_mobility is the scaled mobility constant \f$\mu\f$.
			*	 
			* 	The matrix will be of the form,
			*   \f[ \left[ \begin{matrix}
			*			 \mu^{-1} A & B_{1} \\
			*			 B_{2} & \frac{1}{\Delta t} M + C
			*			 \end{matrix} \right]  \f] 
			*
			* where, 
			* 
			*
			*	 \f[ A(\textbf{p},\textbf{q} ) \; = \; 
			*				\int_{\Omega} \ \textbf{p} \ \cdot \textbf{q} \ dx \f] 
			*
			*	 \f[ M(v,u ) \; = \; 
			*				\int_{\Omega} \ v \ u \ dx \f]
			*
			*	 \f[ B_{1}(v,\textbf{q} ) \; = \; 
			*				\int_{\Omega} \ \nabla \ v  \ \textbf{q} \ dx 
			*				\ + \
			*				\int_{\partial \Omega_{D}} v  \  \textbf{q} 
			*							\ \cdot \boldsymbol \eta \ ds	\f]
			*
			*	 \f[ B_{2}(\textbf{p},u ) \; = \; 
			*				\int_{\Omega} \ \nabla \ \cdot \ \textbf{p} \ u \ dx 
			*				\ + \
			*				\int_{\partial \Omega_{N}} \textbf{p} 
			*							\cdot \boldsymbol \eta \ u \ ds	\f]
			*
			*/
			void
			assemble_local_LDG_cell_and_bc_terms(
				const typename DoFHandler<dim>::active_cell_iterator & cell,
				Assembly::AssemblyScratch<dim>			     & scratch,
				Assembly::DriftDiffusion::CopyData<dim>		 & data,
				const double					     & scaled_mobility_1,
				const double					     & scaled_mobility_2,
				const double 					     & delta_t,
				const double 					     & transient_or_steady,
				const double					     & penalty);

			/** 	\brief Assembles the local LDG flux matrices for interior faces. */
			/** 	These are the matrices which correspond to
			*
			*	 \f[ F_{1}(v,\textbf{q} ) \; = \; 
			*				\langle [[ \ v \ ]] , 
			*			\{ \textbf{q} \} \rangle_{\mathcal{E}_{h}^{0} }  \f]
			*
			*	 \f[ F_{2}(\textbf{p},u ) \; = \; 
			*				\langle  [[ \,  \textbf{p} \, ]] ,  
			*							\{  u \}  \rangle_{\mathcal{E}^{0}_{h} }  \f]
			*
			* It will call LDG_System::LDG::assemble_local_flux_terms and
			* LDG_System::LDG::assemble_local_child_flux_terms.
			*/
			void
			assemble_flux_terms(DoFHandler<dim>		     & carrier_dof_handler,
					ChargeCarrierSpace::CarrierPair<dim> &	carrier_pair,
					FESystem<dim>			     		 & Poisson_fe,
					FESystem<dim>			     		 & carrier_fe);

			/** \brief Assemble the local LDG flux matrices with no local refinement
 			*  This corresponds to the case for two cells where they are on the same 
 			*  refinement level. */ 
			void
			assemble_local_flux_terms(
				Assembly::AssemblyScratch<dim>		    & scratch,
				Assembly::DriftDiffusion::CopyData<dim> & data,
				const double 				    		& penalty);

			/** \brief Function which distributes local fluxes to global matrices.*/
			/*
			* In this function we use the ConstraintMatrix to distribute
			* the local flux matrices to the global system matrix.  
			* Since I have to do this twice in assembling the 
			* system matrix, I made function to do it rather than have
			* repeated code.
			*/

			void 
			distribute_local_fluxes_to_global(
					ChargeCarrierSpace::CarrierPair<dim> 	& carrier_pair,
					Assembly::DriftDiffusion::CopyData<dim>	& data);

			/** Assembles the right hand sides carrier 
 				* test case locally on each cell.  Corresponds to the problem
 				* defined in SOLARCELL::SolarCellProblem::test_steady_state
 				*/
			void
			assemble_local_test_rhs(
				const typename DoFHandler<dim>::active_cell_iterator & cell,
				Assembly::AssemblyScratch<dim>			     & scratch,
				Assembly::DriftDiffusion::CopyData<dim>		 & data,
				const double					     		 & penalty);
	
			/**	Assembles the right hand side for 
 				* test case locally on each cell.  Corresponds to the problem
 				* defined in SOLARCELL::SolarCellProblem::test_transient.
 			*/
			void
			assemble_local_test_transient_rhs(
				const typename DoFHandler<dim>::active_cell_iterator & cell,
				Assembly::AssemblyScratch<dim>			 & scratch,
				Assembly::DriftDiffusion::CopyData<dim>  & data,
				const Vector<double>				     & old_solution,
				const double 					     	 & time,
				const double					    	 & penalty);
			


			/** \brief Computes the local error of your approximation for the 
			*  LDG method on the cell and stores the errors in
			* <code>potential_error<code/> and <code>field_error<code/>.
			*
			* @param solution is the <code>carrier_solution<code/>
			*  vector to the LDG on Poisson's equation!
			* @param potential_error is \f$ L^{2} \f$  error of the approximation 
			*					the density.
			* @param field_error is \f$ L^{2} \f$ error of the approximation 
			*					current.
 			* For problem  defined in SOLARCELL::SolarCellProblem::test_steady_state
 			* 	or defined in SOLARCELL::SolarCellProblem::test_transient.
			*/
			void 
			compute_errors(const Triangulation<dim>	& triangulation,
				       DoFHandler<dim>		& carrier_dof_handler,
				       Vector<double>		& solution,
				       double 				& potential_error,
				       double 				& field_error,
				       const bool			& steady_state,
				       const double 		& time); // const;
	
			/** \brief Computes the local error of your approximation for the 
			*  LDG method on the cell and stores the errors in
			* <code>potential_error<code/> and <code>field_error<code/>.
			*
			* @param solution is the <code>carrier_solution<code/>
			*  vector to the LDG on Poisson's equation!
			* @param potential_error is \f$ L^{2} \f$  error of the approximation 
			*					the density.
			* @param field_error is \f$ L^{2} \f$ error of the approximation 
			*					current.
			*
			*	For problem defined in SOLARCELL::SolarCellProblem::test_DD_Poisson.
			*/
			void 
			compute_coupled_errors(const Triangulation<dim>	& triangulation,
						DoFHandler<dim>		& carrier_dof_handler,
						Vector<double>		& solution,
						double 				& potential_error,
						double 				& field_error,
						const double 		& time); // const;

			/** \brief Computes the local error of your approximation for the 
			*  LDG method on the cell and stores the errors in
			* <code>potential_error<code/> and <code>field_error<code/>.
			*
			* @param solution is the <code>carrier_solution<code/>
			*  vector to the LDG on Poisson's equation!
			* @param potential_error is \f$ L^{2} \f$  error of the approximation 
			*					the density.
			* @param field_error is \f$ L^{2} \f$ error of the approximation 
			*					current.
			*
			*	For problem defined in SOLARCELL::SolarCellProblem::test_interface_coupling.
			*/
			void
			compute_interface_errors(
					const Triangulation<dim>	& triangulation,
					DoFHandler<dim>				& carrier_dof_handler,
					Vector<double>				& solution,
					double 						& potential_error,
					double 						& field_error,
					const double				& time);// const

			/** \brief Prints the current and density of each carrier in carrier pair without units.*/
			void
			output_unscaled_results(DoFHandler<dim>		 & carrier_dof_handler,
					ChargeCarrierSpace::CarrierPair<dim> & carrier_pair,
					const unsigned int 		      		 time_step_number) const;

			/** \brief Prints the current and density of each carrier in carrier pair with units.
 			*	\note The densities will still be without units because Paraview cant handle it otherwise.
 			*	\note I think this will only work when code is compiled in release mode, */
			void
			output_rescaled_results(DoFHandler<dim>			 & carrier_dof_handler,
					ChargeCarrierSpace::CarrierPair<dim> 	 & carrier_pair,
					const ParameterSpace::Parameters 		 & sim_params,
					const unsigned int 			  			 time_step_number) const;

			/** \brief Does the same as LDG_System::LDG::output_unscaled_results 
			but only on the boundary.*/
			void
			output_unscaled_results_on_boundary(DoFHandler<dim>	& carrier_dof_handler,
					ChargeCarrierSpace::CarrierPair<dim> 		& carrier_pair,
					const unsigned int 		   					time_step_number) const;


		private:
			enum
			{
				Interface,
				Dirichlet,
				Neumann,
				Schottky
			};

			// Testing functions	
			const test_Poisson::TrueSolution<dim>			 test_Poisson_solution;
			const test_Poisson::DirichletBoundaryValues<dim> test_Poisson_bc;
			const test_Poisson::RightHandSide<dim>			 test_Poisson_rhs;

		public:
			test_LDG_IMEX::RightHandSide<dim>				test_LDG_rhs;
			test_LDG_IMEX::DirichletBC<dim>					test_LDG_bc;
			test_LDG_IMEX::InterfaceFunction<dim>			test_LDG_interface;
			test_LDG_IMEX::TrueSolution<dim>				test_LDG_solution;						
	
			test_DD_Poisson::DD_RightHandSide<dim>			test_DD_rhs;
			test_DD_Poisson::DD_DirichletBC<dim>			test_DD_bc;
			test_DD_Poisson::DD_TrueSolution<dim>			test_DD_solution;						
		
			test_interface_problem::RightHandSide<dim>		test_interface_rhs;
			test_interface_problem::DirichletBC<dim>		test_interface_bc;
			test_interface_problem::InterfaceFunction<dim>	test_interface_function;
			test_interface_problem::TrueSolution<dim>		test_interface_solution;						
		
			const test_interface_problem::InitialConditions<dim> test_interface_initial;

	};

} 

#endif
