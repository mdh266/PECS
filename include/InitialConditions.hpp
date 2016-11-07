#ifndef _INITIALCONDITIONS_H__
#define _INITIALCONDITIONS_H__

#include<deal.II/base/function.h>

	using namespace dealii;

	/** \brief Initial/Boundary conditions for electrons, \f$ \rho_{n}^{e} \f$.*/
	template<int dim>
	class Electrons_Equilibrium : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			Electrons_Equilibrium() : Function<dim>(dim+1)
			{}

			/** \brief Returns value of \f$\rho_{n}^{e}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;
	};

	/** \brief Initial/Boundary conditions for holes, \f$ \rho_{p}^{e} \f$.*/
	template<int dim>
	class Holes_Equilibrium : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			Holes_Equilibrium() : Function<dim>(dim+1)
			{}

			/** \brief Returns value of \f$\rho_{p}^{e}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;
	};

	/** \brief Initial/Boundary conditions for reductants, \f$ \rho_{r}^{\infty} \f$.*/
	template<int dim>
	class Reductants_Equilibrium : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			Reductants_Equilibrium() : Function<dim>(dim+1)
			{}

			/** \brief Returns value of \f$\rho_{r}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;
	};

	

	/** \brief Initial/Boundary conditions for oxidants, \f$ \rho_{o}^{\infty} \f$.*/
	template<int dim>
	class Oxidants_Equilibrium : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			Oxidants_Equilibrium() : Function<dim>(dim+1)
			{}

			/** \brief Returns value of \f$\rho_{o}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;
	};


#endif
