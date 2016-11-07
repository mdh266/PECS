#ifndef _CARRIER_H__
#define _CARRIER_H__

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/lac/sparse_direct.h>

#include <string>
#include <iostream>

namespace ChargeCarrierSpace
{
	using namespace dealii;


	/** \brief Data structures, functions/parameters, and solvers for a charge carrier. */
	/** This object holds each individual carriers name, system matrix, vectors, solvers, 
	* 	and material properties.
	*
	* 	Each carriers (non-dimensional) transport equation has the form,
	*
	*	\f[ \begin{align}
	* 	u_{t} \  - \  \boldsymbol \nabla \ \cdot  
	* 	\ \mu \left( s \boldsymbol \nabla \Phi u \ + \boldsymbol \nabla u
	*	 \ \right) 
	* 	\; &= \; 
	* 	R(u) + G	&& \text{in} \;  \Omega   \\
 	* 	u \; &= \; u_{D} &&  \text{on} \;  \Omega_{D}     \\
	* 	- \mu \left(s \boldsymbol \nabla \Phi \ u \ + \ \ \boldsymbol 
	*	\nabla u \ \right) 
	* 	\ \cdot \ \boldsymbol \eta 
	*	\;  &= \; K (u) && \text{on} \; \partial \Omega_{N}
	* 	\end{align} \f]
	*
	*/
	template <int dim>
	struct Carrier
	{
		public:
			Carrier();
			
			~Carrier();
	
		/** Sets the name of the carrier. This will be called in constructor 
 		* 	of SOLARCELL::SolarCellProblem.*/		
		void set_name(const std::string & str_name);

		/// Factorizes this objects matrix.
		void set_solver();

		/// Proforms linear solve on the matrix.
		void solve();

		/// Name of the carrier.	
		std::string				name;

		/** This matrix is used to store the carrier's system matrix:
		*   \f[ \left[ \begin{matrix}
		*		\mu^{-1} A & B_{1} + F_{1} \\
		*		 B_{2} + F_{2} & \frac{1}{\Delta t} M + C
		*	 \end{matrix} \right]  \f] 
		*	This matrix is assembled once initially.  See LDG_System::LDG for
		*	 more details.
		*/
		SparseMatrix<double>			system_matrix;
	
		/**  This will be used to store the carrier's system right hand side vector:
		*	\f[
		*	= \; \left( v ,  R(u^{k-1})  + G \right) 
		*	- 
		*	\langle   v, K( u^{k})    \rangle_{\Sigma}  
		*	+
		*	\left( s \textbf{P}  \cdot \boldsymbol \nabla \Phi , u^{k-1} \right)
		*	- 
		*	\langle  \textbf{p}   ,  u_{D}  \rangle_{ \partial \Omega_{N} } 
		*	\f]
		*  	This vector will be assembled at every time step. 
		*  See LDG_System::LDG for more details.
		*/
		Vector<double>				system_rhs;

		/** The solution to this carriers system of linear equations at this time step.*/
		Vector<double>					solution;

		/** The linear solver for this carrier object. */
		SparseDirectUMFPACK			solver;
	
		/** The terms \f$\mu\f$ in the carriers transport equation. */
		double 					scaled_mobility;

		/** The terms \f$\alpha_{i}\f$ in the carriers transport equation. */
		double					charge_number;

//		double 												diff;	
		// SOMETHING ABOUT WHETHER IT HAS RECOMBINATION OR GENERATION ?/
		// OR MOBILITY MODEL???
	

	};

}

#endif
