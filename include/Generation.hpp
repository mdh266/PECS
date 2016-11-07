#ifndef _GENERATION_H__
#define _GENERATION_H__

#include <deal.II/base/function.h>
#include <math.h>

#include "Parameters.hpp"


using namespace dealii;

///////////////////////////////////////////////////////////////////////////////
// Carrier Recombination Functions
///////////////////////////////////////////////////////////////////////////////

/** \brief Function for macroscopic generation of electrons and holes.*/
/** The generation of electrons and holes is modeled as,
* \f[
* \begin{equation}
* G(\textbf{x}) \; = \; \left\{ 
* \begin{array}{lr}
* \alpha (\textbf{x}) 
* \, G_{0} \, e^{- \, \int_{0}^{s} \, \alpha 
* (\, \textbf{x}_{0} \, + \, s' \, \boldsymbol \theta_{0} \, ) \, ds'} 
* \qquad & \qquad \text{if} \quad  \textbf{x} 
* \; = \; \textbf{x}_{0} \, + \, s \, \boldsymbol \theta_{0} \\
* 0 \qquad & \qquad \text{otherwise}
* \end{array} \right. 
* \end{equation} \f]
*
* The point \f$\textbf{x}_{0}\f$ is the photon's incident location and 
* \f$\boldsymbol \theta_{0}\f$ is the incident direction.  The absorption coefficient 
* \f$\alpha(\textbf{x}) \f$ has been averaged over all 
* energy values of light that generate free carriers.  The term 
* \f$G(\textbf{x}_{0}) \; [ \, \text{cm}^{-2} \, \text{s}^{-1} \, ]\f$ 
* represents the surface photon flux at the point \f$\textbf{x}_{0}\f$. 
*
* \note The incident direction
* is always downward in the y-direction for now, so the incident location
* will always be the the top of the domain.
*/
	template <int dim>
	class Generation : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			Generation() : Function<dim>(1)
			{}
	
			/** \brief Sets all the parameter values to zero.*/
			/** This sets all the parameter values to zero so that when you call,
			*		this objects value at a point it will always return 0.*/
			void	set_dark_params();

			/** \brief Sets the objects scaled values from the parameters.*/
			void	set_illuminated_params(const ParameterSpace::Parameters & Params);
	
			/** \brief Returns the value of the generation function at this points.*/		
			virtual double value(const Point<dim> &p, 
					     const unsigned int component = 0 ) const;
			
		private:
			/** \brief The scaled photon flux.*/
			double scaled_photon_flux;
			/** \brief The scaled absorption coefficent.*/
			double scaled_absorption_coeff;
			/** \brief The scaled incident location.*/
			double scaled_incident_location;		
	};

#endif
