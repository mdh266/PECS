///
// Poisson solver using Local Discontinuous Galerkin Method
//
//
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/std_cxx11/bind.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h> // for block structuring
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h> // for block structuring
#include <deal.II/lac/sparse_direct.h>

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

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>

#include <fstream>
#include <iostream>
#include <string>

// multithreading
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/thread_management.h>

#include "Parameters.hpp"
#include "ParameterReader.hpp"
#include "Grid.hpp"
#include "Assembly.hpp"
#include "MixedFEM.hpp"
#include "CarrierPair.hpp"
#include "Poisson.hpp"
#include "InitialConditions.hpp"
#include "Generation.hpp"
#include "BiasValues.hpp"
#include "LDG.hpp"
#include "PostProcessor.hpp"

/** \namespace SOLARCELL This namespace is for everything that has to do particularly 
 *	with the solar cell model is the main namespace for everything that is
 *	constained in the source file for	SOLARCELL:SolarCellProblem. */
namespace SOLARCELL
{
	using namespace dealii;

	/** \brief The recombination model.*/
	/** Schockley-Reed-Hall Recombination functional
	* \f[ \begin{equation}
	*	R(\rho_{n}, \rho_{p})
	* \; = \;
	* \frac{\rho_{i}^{2} \, -\, \rho_{n} \, \rho_{p}}
	* {\tau_{n} \, (\rho_{n} \, + \, \rho_{i}) \, + \, \tau_{p} \, ( \rho_{p} \, + \, \rho_{i})}.
	* \end{equation}
 	* \f]
 	* 
 	* The term \f$\rho_{i}\f$ is the intrinsic electron density and can be spatially varying. 
 	* \f$\tau_{n}, \ \tau_{p} \f$ are constants called the electron and hole lifetimes.  
	* @param Holds the intrinsic density and recombination times.
	* @param electron_density \f$\rho_{n}^{k-1}\f$.
	* @param hole_density \f$\rho_{p}^{k-1}\f$.
	*/
	inline	double SRH_Recombination(const double & electron_density,
					 const double & hole_density,
					 const ParameterSpace::Parameters & params)
	{
		return /*0.0;*/ ((params.scaled_intrinsic_density *
				params.scaled_intrinsic_density - 
				electron_density * hole_density) /
				(params.scaled_hole_recombo_t * (electron_density +
				 params.scaled_intrinsic_density) + 
				(params.scaled_electron_recombo_t * (hole_density -
				 params.scaled_intrinsic_density)) ) );	

	}

///////////////////////////////////////////////////////////////////////////////
// Solar Cell Problem Class
///////////////////////////////////////////////////////////////////////////////


	template <int dim>
	class SolarCellProblem
	{
		public:
			/** The constructor for this problem.  The parameter values are set here
 			* and the generation function is also set here.*/
			SolarCellProblem(const unsigned int degree,
					 ParameterHandler	&param);

			/** Destructor just deletes the dof_handlers.*/
			~SolarCellProblem();

			/** Computes the parallel IMEX method on the semiconductor-electrolyte problem.*/
			void 
			run_full_system();

			/** \brief Solves test Poisson problem and prints convergence rates.*/
			/** Runs mixed method and LDG method sovler for a test Poisson problem
 			*	with Dirichlet and Neumann boundary conditions.  Computes solutions over a 
 			*	sequence of refined meshes and prints out the convergence rate of the
 			* methods. The problem is solves is,
 			*	\f[
			*	\begin{align}
			*	-\Delta \Phi(\textbf{x})  \; &= \; f(\textbf{x}) && \text{in} \quad [0,1] \times [0,1] \\
			*	-\nabla \Phi(\textbf{x}) \cdot \boldsymbol \eta \; 
			*	&= \; 0 && \text{on} \quad [0,1] \times \{y \, = \, 0, 1\}  \\
			*	\Phi(\textbf{x}) \; &= \; \Phi_{D}(\textbf{x)} && \text{on} 
			*	\quad \{x \, = \, 0, 1\} \times [0,1]
			*	\end{align} \f]		
			* We use the manufactured solution,
			*	\f[
			*	\Phi(x,y) \; = \; \cos(2 \pi y) - \sin(2 \pi x) - x,  
			*	\f]
			* with right hand side function.
			*	\f[
			*	f(x,y) \; = \; 4 \pi^{2} \, \left( \cos(2 \pi y) - \sin(2 \pi x) \right).
			* \f] 
			* \note This is used by the LDG and mixed method. */
			void
			test_steady_state(const unsigned int 	& n_refine,
					 ConvergenceTable 	& Mixed_table,
					 ConvergenceTable 	& LDG_table);

			/** \brief Tests the LDG method with IMEX time stepping on the drift-diffusion eq.*/
			/** Solves the drift-diffusion equation using the LDG method with IMEX time
 			* stepping.  The drift-diffusion equation has a fixed electric field and
 			* Dirichlet boundary conditions. Computes the solution over a sequence of 
 			* refined meshes and prints out the convergence rate. 
 			*
			* To test the time-dependent linear drift-diffusion equation we solve the problem,
			*  \f[ \begin{align}
			* u_{t} + \boldsymbol \nabla \cdot \left( \textbf{E} u -
			* \boldsymbol \nabla u \right) \; &= \; f(\textbf{x},t) && 
			* \text{in} \quad  \Omega \times (0,T] \\
			* \left( \textbf{E} u - \boldsymbol \nabla u \right) \cdot \boldsymbol \eta \; 
			* &= \; g_{I}(\textbf{x},t) && \text{on} \quad \Gamma_{I} \times (0,T]  \\
			* u \; &= \; g_{D}(\textbf{x},t) && \text{on} \quad \Gamma_{D} \times (0,T]  \\
			* u \; &= \; u_{0}(\textbf{x}) && \text{in} \quad \Omega \times \{t=0\}. 
			* \end{align} \f]
			*	
			* We take the domain to be \f$\Omega = [0,1] \times [0,1]\f$ with boundaries 
			* \f$\Gamma_{I} = \{x \, = \, 1 \} \times [0,1]\f$, 
			* and \f$\Gamma_{D} = \{x \, = \, 0\} \times [0,1] \cup [0,1] \times \{y \, = \, 0,1 \}\f$.  
			* For simplicity we assumed 
			* \f$\textbf{E}(x,y) = \langle 1, 0 \rangle\f$ and \f$T=1\f$.  The manufactured solution is,
			* \f[ u(x,y,t) \; = \; e^{-t} + \cos(2 \pi x) + \cos(2  \pi y). \f]
			* This has a corresponding right hand side function,
			* \f[ 
			* f(x,y,t) \; = \; -e^{-t} + 4\pi^{2} \cos(2 \pi x) + 4\pi^{2} \cos(2 \pi y) + 2 \pi \sin(2\pi x) .
			* \f]
			* The Dirichlet boundary conditions \f$g_{D}(x,y,t)\f$ are taken to be \f$u(x,y,t)\f$ on
			* \f$\Gamma_{D}\f$, while the Neumann boundary condition function is,
			* \f[ g_{N}(x,y,t) \; = \; -e^{-t} - \cos(2 \pi y) - 1.
			* \f]
			* The initial condition \f$u_{0}\f$ is the \f$L^{2}\f$ projection of the solution \f$u(x,y,0)\f$ 
			* onto the DG basis. The time step is chosen to be \f$ \Delta t = h^{k+1}\f$ so as to not 
			* ruin the spatial discritization convergence rate.*/
			void
			test_transient(const unsigned int & n_refine,
				       ConvergenceTable   & LDG_table);

			/** \brief Test the coupling of the mixed FEM to the LDG method.*/
			/** The coupling of the drift-diffusion-Poisson problem using a LDG-IMEX and 
			* a mixed finite element method is tested on the problem,
			* \f[ \begin{align}
			* -\Delta \Phi(\textbf{x},t)  \; &= \; C(\textbf{x},t) - u & \qquad \qquad
			*  & \text{in} \quad  \Omega \times (0,T]  \\
			* \Phi(\textbf{x}) \; &= \; \Phi_{D}(\textbf{x)} & & \text{on}
			* \quad \partial  \Omega \times (0,T]  \\
			* u_{t} + \boldsymbol \nabla \cdot \left( -\nabla \Phi \, u - \boldsymbol \nabla u \right) 
			* \; &= \; f(\textbf{x},t) & \qquad \qquad & \text{in} \quad  \Omega \times (0,T] \\
			* u \; &= \; g_{D}(\textbf{x},t) & & \text{on} \quad \partial \Omega \times (0,T]  \\
			* u \; &= \; u_{0}(\textbf{x}) & & \text{in} \quad \Omega \times \{t=0\}. 
			* \end{align} \f]
			*
			* We take the domain to be $\Omega = [0,1] \times [0,1]$.  The solutions are, 
			* \f[ \Phi(x,y) \; = \; \cos(2 \pi y) - \sin(2 \pi x) - x,  \f]
			*
			* and,
			* \f[
			* u(x,y,t) \; = \; e^{-t} + \cos(2 \pi x) + \cos(2  \pi y).
			* \f]
			* We take the right hand side functions to be,
			* \f[ C(x,y,t) \; = \;  4 \pi^{2} \, \left( \cos(2 \pi y) - \sin(2 \pi x) \right) 
			* +  
			* e^{-t} + \cos(2 \pi x) + \cos(2  \pi y),
			* \f]
			* 
			* and,
			* \f[
			* \begin{align}
			* f(x,y,t) \; = \;
			* & 4 \pi^{2} \left(\cos(2 \pi y) - \sin(2 \pi x) \right) 
			* \left( e^{-t} + \cos(2 \pi x) + \cos(2  \pi y) \right) \\
			* &  - 2 \pi \left( 2 \pi \cos(2 \pi x) + 1 \right)  \sin(2 \pi x)  
			* - 4 \pi^{2} \sin(2 \pi y)^{2}  \\
			* & -e^{-t} + 4\pi^{2} \cos(2 \pi x) +  4\pi^{2} \cos(2  \pi y).
			* \end{align}
			* \f]
			*
			* The Dirichlet boundary condition \f$g_{D}(x,y,t)\f$ and \f$\Phi_{D}(x,y)\f$ 
			* are taken to be the solutions \f$u(x,y,t)\f$ on and \f$\Phi(x,y,t)\f$ 
			* respectively. The initial condition \f$u_{0}\f$ is the \f$L^{2}\f$ projection 
			* of the solution \f$u(x,y,0)\f$ onto the DG basis. To perform time stepping 
			* we use a first-order IMEX time strategy.  In order to obtain the underlying
			* errors of the LDG method we take time steps \f$\Delta t = h^{k+1}\f$, 
			* integrating from \f$t=0\f$ to \f$T=1\f$ for polynomials of order \f$k\f$. 
			* We see that this choice of time step causes the super convergence results 
			* for \f$\textbf{D}_{h}\f$ to be deteriorated, however they will be restored
			* with the choice of \f$\Delta t = h^{k+2}\f$.
			*/
			void
			test_DD_Poisson(const unsigned int & n_refine,
					ConvergenceTable & Mixed_table,
					ConvergenceTable & LDG_table);

												

			/** \brief Tests the LDG method for a pair of parabolic problems coupled across the interface.*/
			/** Computes the solutions over a sequence of refined meshes and prints out 
			* the convergence rate. The problem is for
			* \f$u(\textbf{x},t) \in \Omega_{S} \times [0,T]\f$ and 
			* \f$v(\textbf{x},t) \in \Omega_{E} \times [0,T] \f$ that are coupled across 
			* the interface \f$\Sigma \f$ we solve the coupled parabolic problems,
			* \f[ \begin{align}
			*	u_{t} - \boldsymbol \nabla \cdot \left(  \boldsymbol \nabla u \right) \; 
			*	&= \; f_{1} (\textbf{x},t) && \text{in} \quad  \Omega_{S} \times (0,T] \\
			*	-\left( \boldsymbol \nabla u \right) \cdot \boldsymbol \eta_{\Sigma} \; 
			*	&= \; u(\textbf{x},t)v(\textbf{x},t) - I(\textbf{x},t) && \text{on} \quad \Sigma \times (0,T]\\
			*	u \; &= \; g_{1,D}(\textbf{x},t) && \text{on} \quad \Gamma_{S,D} \times (0,T]  \\
			* u \; &= \; u_{0}(\textbf{x}) && \text{in} \quad \Omega_{S} \times \{t=0\}. 
			*	\end{align} \f]
			*	and,
			* \f[ \begin{align}
			*	v_{t} - \boldsymbol \nabla \cdot \left( \boldsymbol \nabla v \right) \; 
			*	&= \; f_{2} (\textbf{x},t) && \text{in} \quad  \Omega_{E} \times (0,T] \\
			*	-\left(\boldsymbol \nabla v \right) \cdot \boldsymbol \eta_{\Sigma} \; 
			*	&= \; -u(\textbf{x},t)v(\textbf{x},t) + I(\textbf{x},t) 
			*	&& \text{on} \quad \Sigma \times (0,T]  \\
			*	v \; &= \; g_{2,D}(\textbf{x},t) && \text{on} \quad \Gamma_{E,D} \times (0,T]\\
			* v \; &= \; v_{0}(\textbf{x}) && \text{in} \quad \Omega_{E} \times \{t=0\}. 
			* \end{align} \f]
			*
			* We take the domain to be  \f$\Omega = [0,1] \times [0,1]\f$ with 
			* \f$\Omega_{S} = [0,1/2] \times [0,1]\f$ and 
			* \f$\Omega_{E} = [1/2,1] \times [0,1]\f$. The interface is
			* \f$\Sigma = \{x \, = \, 1/2 \} \times [0,1]\f$ and the boundaries are 
			* \f$\Gamma_{1,D} = \{x \, = \, 0\} \times [0,1] \cup [0,1/2] \times \{y \, = \, 0,1 \}\f$ 
			* and \f$\Gamma_{2,D} = \{x \, = \, 1\} \times [0,1] \cup [1/2,1] \times \{y \, = \, 0,1 \}\f$. 
			* We remark that the interface vector
			* \f$\boldsymbol \eta_{\Sigma} = \langle 1, 0 \rangle\f$
			* \newline The manufactured solutions for this problem are,
			* \f[ \begin{equation}
			*	u(x,y,t) \; = \; v(x,y,t) \; = \; e^{-t} + \cos(2 \pi x) + \cos(2  \pi y).
			* \end{equation} \f]
			* The right hand side functions are,
			* \f[ \begin{equation}
			*	f_{1} (x,y,t) \; = \; f_{1} (x,y,t)  \; = \; -e^{-t} + 4\pi^{2} \cos(2 \pi x) 
			*	+ 4\pi^{2} \cos(2 \pi y) + 2 \pi \sin(2\pi x) .
			*	\end{equation} \f]
			* The Dirichlet boundary conditions \f$g_{1,D}(x,y,t)\f$ is taken to be \f$u(x,y,t)\f$ 
			* on \f$\Gamma_{1,D}\f$ and \f$g_{2,D}(x,y,t)\f$ is taken to be \f$v(x,y,t)\f$ on 
			* \f$\Gamma_{2,D}\f$.  The interface condition function is,
			*\f[ \begin{equation}
			* I(x,y,t) \; = \; \left(e^{-t} + \cos(\pi y) - 1\right)^{2}.
			* \end{equation} \f]
			*	The initial conditions \f$ u_{0} \f$, \f$ v_{0} \f$ are the \f$ L^{2} \f$ 
			* projection of the 
			*	solutions \f$ u(x,y,0) \f$ and \f$v(x,y,0)\f$ onto the DG basis. 
			*  In order to obtain the underlying errors of the LDG method we take time 
			* steps \f$\Delta t = h^{2}\f$ for \f$k=1\f$ and \f$\Delta t = h^{3}\f$ 
			* for \f$k=2\f$.  \note \f$h\f$ will be the 
			* same value for both of the triangulations of \f$\Omega_{S}\f$ and \f$\Omega_{E}\f$ 
			* since they are basically the same mesh. */
			void
			test_interface_coupling(const unsigned int  & n_refine,
									ConvergenceTable    & Table);


		private:
			/// Degree of the polynomials used.
			const unsigned int degree;
			/// Time step.
			double delta_t;
			/**  Boolean that keeps track of whether we are solving semiconductor-Poisson 
 			* 	problem or semiconductor-electrolyte-Poisson. This is for testing. */
			bool full_system;

	
			enum
			{
				/// Interface boundary id set in Grid_Maker::Grid.
				Interface,
				/// Dirichlet boundary id set in Grid_Maker::Grid.
				Dirichlet,
				/// Neumann boundary id set in Grid_Maker::Grid.
				Neumann,
				/// Shottky boundary id set in Grid_Maker::Grid.
				Schottky
			};
	
			/** sim_parms holds all the parameter values.*/
			ParameterSpace::Parameters			sim_params;
	
			/** The parameter handlers reads in the paramter values from the input file
 			*		and assigns them to sim_params. */
			ParameterSpace::ParameterHandler		&prm;

			/*--------------------------------------------------------------*/
			/* 	Poisson Data Structures					*/
			/*--------------------------------------------------------------*/
			/// The Poisson's equation mesh.
			Triangulation<dim>				Poisson_triangulation;

			/// DofHandler for the Poisson equation.
			DoFHandler<dim>					Poisson_dof_handler;	
		
			/// This will be there FE -> Raviart-Thomas + DG
			FESystem<dim>					Poisson_fe;

			/// Holds the matrices, vectors, etc for Mixed method of Poisson.
			Poisson::PoissonData<dim>			Poisson_object;			
			
			/*--------------------------------------------------------------*/
			/* 	Carrier Data Structures					*/
			/*--------------------------------------------------------------*/
			/// Mesh for the semiconductor domain.
			Triangulation<dim>				semiconductor_triangulation;

			/// DoFHandler for electron and hole. They will have the same distribution of dofs.
			DoFHandler<dim>					semiconductor_dof_handler;

			/// Holds the matrices, vectors etc for the electron/hole pair.
			ChargeCarrierSpace::CarrierPair<dim> 		electron_hole_pair;

			/// Mesh for the electrolyte domain.
			Triangulation<dim>				electrolyte_triangulation;
			/// DoFHandler for the redox pair. They will have same distribution od dofs.

			DoFHandler<dim>					electrolyte_dof_handler;
			/// Holds the matrices, vectors etc for the reductant/oxidant pair.

			ChargeCarrierSpace::CarrierPair<dim> 		redox_pair;

			/// The finite element of the LDG Method DG^{2} and DG
			FESystem<dim>					carrier_fe;

			/** Object holds the assembling and printing routines for the LDG/
			* semiconductor/electrolyte carrier system. */
			MixedPoisson::MixedFEM<dim>			Mixed_Assembler;

			/** Object holds the assembling and printing routines for the Mixed FEM
			* Poisson equation on the whole domain.*/
			LDG_System::LDG<dim>				LDG_Assembler;	



			/*----------------------------------------------------------------*/
			/* Mappings														  */
			/*----------------------------------------------------------------*/
			/** \brief Vector of semiconductor cells on the interface: their level and index. */
			std::vector<std::pair<unsigned int, unsigned int>>  semi_interface_cells;
			/** \brief Vector of semiconductor cells' faces that are on the interface. */
			std::vector<unsigned int>  			semi_interface_faces;
	
			/** \brief Vector of semiconductor cells on the interface: their level and index. */
			std::vector<std::pair<unsigned int, unsigned int>>  elec_interface_cells;
			/** \brief Vector of semiconductor cells' faces that are on the interface. */
			std::vector<unsigned int>  			elec_interface_faces;

			/** \brief Bijective mapping from semiconductor cell to  global interface face index.*/
			std::map<std::pair<unsigned int, unsigned int>, unsigned int>	semi_interface_map;
			/** \brief Bijective mapping from electrolyte cell to  global interface face index.*/
			std::map<std::pair<unsigned int, unsigned int>, unsigned int>	elec_interface_map;


			/** \brief Bijective mapping from semiconductor cells to Poisson cells.*/
			std::map<std::pair<unsigned int, unsigned int>, 
			std::pair<unsigned int, unsigned int>> s_2_p_map;
			/** \brief Bijective mapping from electrolyte cells to Poisson cells.*/
			std::map<std::pair<unsigned int, unsigned int>, 
			std::pair<unsigned int, unsigned int>> e_2_p_map;

			/*-------------------------------------------------------------*/
			/*	Initial Conditions										   */
 			/*-------------------------------------------------------------*/
			/// \f$ \rho_{n}^{e}(\textbf{x}) \f$
			const Electrons_Equilibrium<dim>			electrons_e;
			/// \f$ \rho_{p}^{e} (\textbf{x}) \f$
			const Holes_Equilibrium<dim>				holes_e;
			/// \f$ \rho_{r}^{\infty}(\textbf{x})  \f$
			const Reductants_Equilibrium<dim>			reductants_e;
			/// \f$ \rho_{o}^{\infty} (\textbf{x}) \f$
			const Oxidants_Equilibrium<dim>				oxidants_e;

			/*-------------------------------------------------------------*/
			/* The potential functions					*/
			/*-------------------------------------------------------------*/
			/// \f$ \Phi_{\text{bi}} \f$
			Built_In_Bias<dim>					built_in_bias;
			/// \f$ \Phi_{\text{app.}} \f$
			Applied_Bias<dim>					applied_bias;
			/// \f$ \Phi^{\infty} \f$
			Bulk_Bias<dim>						bulk_bias;
			/// \f$ \Phi_{\text{Sch}} \f$ 	
			Schottky_Bias<dim>					schottky_bias;

			/*-------------------------------------------------------------*/
			/* Generation function						*/
 			/*-------------------------------------------------------------*/
			/** Generation function \f$G(\textbf{x})\f$.  The incident location
 			*		is set here.*/
			Generation<dim>						 generation;


			Vector<double>	diff_electrons;
			Vector<double>	diff_holes;
			Vector<double>	diff_reductants;
			Vector<double>	diff_oxidants;
			Vector<double>	diff_Poisson;

			
			/** Distributes the dofs and allocates memory for matrices and vectors.*/	
			void
			setup_dofs();	
	
			/** Prints the number of cells and dofs.*/
			void
			print_sim_info();

			/** \brief Create bijective mappings for assembling interface coupling terms.*/
 			/** The bijection is from cells' faces that are on the  interface to 
 			*  the corresponding global interface face index.
 			*  If we are on a semiconductor cell's interface face we can get the
			* corresponding face on the electrolyte triangulation from elec_interface_cell
			* and elec_interface_face. And vice versa. */
			void
			setup_mappings();
	
			/** Assembles the matrix corresponding to the mixed method for Poisson eq.
			* See MixedPoisson::MixedFEM for more information.*/
			void 
			assemble_Poisson_matrix();

			/** Copys the local matrices assembled my and MixedFEM object into
			* the global matrix which is in Poisson_object.  
			* Used for assembly in parallel by workstream.
			*/
			void 
			copy_local_to_global_Poisson_matrix(
						const Assembly::Poisson::CopyData<dim> 	& data);

			/** Copys the local right hand side vectors assembled my and this
			*  object into the global right hand side vector in Poisson_object. 
			*  Used for assembly in parallel by workstream.
			*/
			void 
			copy_local_to_global_Poisson_rhs(
						const Assembly::Poisson::CopyData<dim> & data);

			/** \brief Assembles the global poisson rhs for the coupled problem.*/
 			/** It first uses the semiconductor_dof_handler to run over the 
 			* cells in the semiconductor triangulation and calls 
 			* SOLARCELL::SolarCellProblem::assemble_local_Poisson_rhs_for_semiconductor 
 			* using workstream. It then uses the electrolyte_dof_handler to run over the cells
 			* in the electrolyte triangulation and calls 
 			* SOLARCELL::SolarCellProblem::assemble_local_Poisson_rhs_for_electrolyte
 			* using workstream.
 			*
 			* The general problem for the Poisson equation is, 
			* \f[ \begin{align}
			*	\epsilon_{r}^{-1} \ \textbf{D} \ + \ \nabla \Phi \ &= \ 0  && \text{in} \; \Omega \\
			*	\ \lambda^{2} \  \nabla \ \cdot \ \textbf{D} \  &= 
			*	 \ f(\textbf{x})					&& \text{in} \; \Omega \\
			*		\textbf{D} \ \cdot \ \boldsymbol \eta \ &= \ 0 
			*						&& \text{on} \; \partial \Omega_{N} \\
			*		\Phi \ &=  \ \Phi_{\text{app}} \ + \ \Phi_{\text{bi}}
			*							 && \text{on} \;  \Gamma_{S} \\
			*		\Phi \ &=  \ \Phi^{\infty}
			*							 && \text{on} \;  \Gamma_{E} 
			*   \end{align} \f]
			*
			* For \f$\lambda^{2} \ = \ \frac{\Phi^{*} \epsilon }{L^{2} q C^{*}}\f$.  .  
			* \f$f(\textbf{x})\f$ is defined as
			* \f[ 
			* f(\textbf{x}) \ = \
			* \left\{ \begin{array}{cc}
			* \left[ \ N_{D} \ - \ N_{A} \ -\  \left( \rho_{n} \ - \ \rho_{p} \right) \ \right]   
			* & \text{in} \; \Omega_{S} \\
			*  -\  \left( \rho_{r} \ - \ \rho_{o} \right)   
			* & \text{in} \; \Omega_{E}
			* \end{array} \right.
			* \f] 
			* 
			* See MixedPoisson::MixedFEM on how to build the corresponding left hand side of the
			* weak formulation.  This function essentially constructs:
			*
			* \f[ = \ -\langle \textbf{p} \ , \ \Phi_{\text{app}} + \Phi_{\text{bi}} 
			* \rangle_{\Gamma_{S}} - \langle \textbf{p} \ , \Phi^{\infty} \ \rangle_{\Gamma_{E}}   
			*	- \
			* \left(  v,\ f(\textbf{x}) \right)_{\Omega} \f]
			* 
			*
			*  For all \f$( v \  , \ \textbf{p}  ) \, \in \, W \, \times\, \textbf{V}^{d}\f$. 
			*/
			void
			assemble_Poisson_rhs();

		
			/**\brief Assembles the local poisson rhs for the coupled problem
 			*		in the semiconductor triangulation. */
			/**  This function loops through all the cells in the 
 			*	semiconductor triangulation and performs the following calculation, 
 			*	\f[ \ = \ -\langle \textbf{p} \ , \ \Phi_{\text{app}} + \Phi_{\text{bi}} 
			* \rangle_{\partial \Omega_{e} \cap \Gamma_{S}} 
			*	-
			*	\left( v, 
			* \left[ \ N_{D} \ - \ N_{A} \ -\  \left( \rho_{n} \ - \ \rho_{p} \right) \ \right]   
			* \ \right)_{\Omega_{e}} \f]
			* 
			*
			* For all \f$( v \  , \ \textbf{p}  ) \, \in \, W \, \times\, \textbf{V}^{d}\f$ and
			* all \f$\Omega_{e} \in \Omega_{S}\f$. It stores the data in 
			*	Assembly::Poisson::CopyData.
			*/
			void
			assemble_local_Poisson_rhs_for_semiconductor(
				const typename DoFHandler<dim>::active_cell_iterator & cell,
				Assembly::AssemblyScratch<dim>						 & scratch,
				Assembly::Poisson::CopyData<dim>					 & data);
	
			/** \brief Assembles the local poisson rhs for the coupled problem
 			*		in the electrolyte triangulation. */
			/**  This function loops through all the cells in the 
 			*	electrolyte triangulation and performs the following calculation, 
 			*	\f[ \ = \ -\langle \textbf{p} \ , \ \Phi^{\infty} 
			* \rangle_{\partial \Omega_{e} \cap \Gamma_{E}} 
			*	- 
			*	\left( v, 
			* \left[  -\  \left( \rho_{r} \ - \ \rho_{o} \right) \ \right]   
			* \ \right)_{\Omega_{e}} \f]
			* 
			*
			*  For all \f$( v \  , \ \textbf{p}  ) \, \in \, W \, \times\, \textbf{V}^{d}\f$ and
			*  all \f$\Omega_{e} \in \Omega_{E}\f$. It stores the data in 
			*  Assembly::Poisson::CopyData. */
			void 
			assemble_local_Poisson_rhs_for_electrolyte(
				const typename DoFHandler<dim>::active_cell_iterator & cell,
				Assembly::AssemblyScratch<dim>						 & scratch,
				Assembly::Poisson::CopyData<dim>					 & data);


			/** Assembles the LDG system that asembles the mass matrix and system matrix for the
 			* electrons, holes, reductants and oxidants using the IMEX method and LDG 
 			* discritizations.*/	
			void
			assemble_LDG_system(const double & transient_or_steady);
			
			/** Copies the local calculations into the global mass matrix for the
 			*		electron hole pair. */
			void
			copy_local_to_global_semiconductor_mass_matrix(						
					const Assembly::DriftDiffusion::CopyData<dim>	& data);

			/** Copies the local calculations into the global mass matrix for the
 			*		redox pair. */
			void
			copy_local_to_global_electrolyte_mass_matrix(						
					const Assembly::DriftDiffusion::CopyData<dim>	& data);


			/** Copies the local calculations into the global system matrix for the
 			*		electron hole pair. */
			void
			copy_local_to_global_semiconductor_system_matrix(						
					const Assembly::DriftDiffusion::CopyData<dim>	& data);

			/** Copies the local calculations into the global system matrix for the
 			*		redox pair. */
			void
			copy_local_to_global_electrolyte_system_matrix(						
					const Assembly::DriftDiffusion::CopyData<dim>	& data);

			/** Copies the local calculations into the global right hand side for the
 			*		electron hole pair. */
			void
			copy_local_to_global_semiconductor_system_rhs(						
					const Assembly::DriftDiffusion::CopyData<dim>	& data);

			/** Copies the local calculations into the global right hand side for the
 			*		redox pair. */
			void
			copy_local_to_global_electrolyte_system_rhs(						
					const Assembly::DriftDiffusion::CopyData<dim>	& data);

			/** Builds the RHSes for the electron and hole equations.*/
			/**	This function loops through all the cells in the semiconductor
 			*	triangulation. It does so by using semiconductor_dof_handler and 
 			*	calling SOLARCELL::SolarCellProblem::assemble_local_semiconductor_rhs 
 			*	uising workstream. For more information on the LDG method used in this 
 			*	caee see LDG_System::LDG.  */
			void
			assemble_semiconductor_rhs();
		
			/** Builds the local RHSes for the electron and hole equations.*/
			/** It stores the calculated data in Assembly::DriftDiffusion::CopyData
 			*	the  rhs at time \f$k\f$ is,
 			*	\f[ = 
			*	\left( v ,  R(u^{k-1})  + G(\textbf{x}) \right)_{\Omega_{e}} - 
			*	\langle   v, K( u^{k-1})    \rangle_{\partial \Omega_{e} \cap \Sigma}  
			*	 +
			*	\left( s \textbf{P}  \cdot \boldsymbol \nabla \Phi , u^{k-1} 
			*	\right)_{\Omega_{e}}
			*	- 
			*	\langle  \textbf{p} ,  u_{D} \rangle_{ \partial \Omega_{e} \cap \Gamma_{S} }
			*	\f]
			*
			*  For all \f$(v,\textbf{p}) \in W \times \textbf{W}^{d})\f$ and all 
			*  \f$\Omega_{e} \in \Omega_{S} \f$.
 			*/
			void
			assemble_local_semiconductor_rhs(
						const typename DoFHandler<dim>::active_cell_iterator & cell,
						Assembly::AssemblyScratch<dim>						 & scratch,
						Assembly::DriftDiffusion::CopyData<dim>		 		 & data,
						const double 					 					 & penalty);

			/** Builds the RHSes for the reductant and oxidant equations.*/
			/**	This function loops through all the cells in the electrolyte
 			*	triangulation. It does so by using electrolyte_dof_handler and 
 			*	calling SOLARCELL::SolarCellProblem::assemble_local_electrolyte_rhs 
 			*	uising workstream.*/
			void
			assemble_electrolyte_rhs();
		
			/** Builds the local RHSes for the electron and hole equations.*/
			/** It stores the calculated data in Assembly::DriftDiffusion::CopyData
 			*	the  rhs at time \f$k\f$ is,
 			*	\f[ = 		
			*	- \langle   v, K( u^{k-1})    \rangle_{\partial \Omega_{e} \cap \Sigma}  
			*	 +
			*	\left( s \textbf{P}  \cdot \boldsymbol \nabla \Phi , u^{k-1} 
			*	\right)_{\Omega_{e}}
			*	- 
			*	\langle  \textbf{p} ,  u_{D} \rangle_{ \partial \Omega_{e} \cap \Gamma_{E} }
			*	\f]
			*
			*  For all \f$(v,\textbf{p}) \in W \times \textbf{W}^{d})\f$ and all 
			*  \f$\Omega_{e} \in \Omega_{E} \f$.
 			*/
			void
			assemble_local_electrolyte_rhs(
					const typename DoFHandler<dim>::active_cell_iterator & cell,
					Assembly::AssemblyScratch<dim>			 & scratch,
					Assembly::DriftDiffusion::CopyData<dim>  & data,
					const double 					 		 & penalty);


			/** Factorizes all the matrices (Poisson and carriers). */	
			void
			set_solvers();

			/** Solves Poisson equation.*/
			void 
			solve_Poisson();

			/** Solve the linear systems of electron and holes using one thread for each.*/
			void
			solve_semiconductor_system();

			/** Solve the linear systems for all the carriers.  Using multi-threading so that
 			* each thread solves a linear system.*/
			void
			solve_full_system();

			/** Solve the linear systems of reductants and oxidants using one thread for each.*/
			void
			solve_electrolyte_system();

			/** Print the results into three .vtu files using multi-threading. 
 			* One thread prints poisson, one thread prints electron/holes, one thread prints
 			* reductant/oxidant.*/
			void
			print_results(unsigned int time_step_number);
	
			/** Multi threaded printing that only prints the electron-hole and reductant-oxidant
 			* values on the boundary of their respective domains. */
			void
			print_results_on_boundary(unsigned int time_step_number);

			/** Assembles the local cell rhs term for the LDG method applied to the 
			* drift-diffusion equation defined in
 			*	SOLARCELL::SolarCellProblem::test_DD_Poisson.*/
			void
			assemble_local_coupled_DD_test_rhs(
					const typename DoFHandler<dim>::active_cell_iterator & cell,
					Assembly::AssemblyScratch<dim>			 & scratch,
					Assembly::DriftDiffusion::CopyData<dim>	 & data,
					const double					 		 & time,
					const double 					 		 & penalty);								
			/** Assembles the local cell rhs term for the mixed FEM applied to the 
			* Poisson equation defined in
 			*	SOLARCELL::SolarCellProblem::test_DD_Poisson.*/
			void
			assemble_local_coupled_Poisson_test_rhs(
					const typename DoFHandler<dim>::active_cell_iterator & cell,
					Assembly::AssemblyScratch<dim>			& scratch,
					Assembly::Poisson::CopyData<dim>	 	& data,
					const double							& time);
	
			/** Assembles the local cell rhs term for the \$f u \f$ in problem defined in
 			*	SOLARCELL::SolarCellProblem::test_interface_coupling.*/
			void
			assemble_local_test_semiconductor_rhs(
					const typename DoFHandler<dim>::active_cell_iterator & cell,
					Assembly::AssemblyScratch<dim>		 				 & scratch,
					Assembly::DriftDiffusion::CopyData<dim>	 			 & data,
					const double										 & time,
					const double 				 						 & penalty);								

			/** Assembles the local cell rhs term for the \f$ v \f$ in problem defined in
 			*	SOLARCELL::SolarCellProblem::test_interface_coupling.*/
			void
			assemble_local_test_electrolyte_rhs(
					const typename DoFHandler<dim>::active_cell_iterator & cell,
					Assembly::AssemblyScratch<dim>						 & scratch,
					Assembly::DriftDiffusion::CopyData<dim>				 & data,
					const double 										 & time,
					const double 					 					 & penalty);									


	}; // END CLASS

} // END NAMESPACE
