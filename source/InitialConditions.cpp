#include "../include/InitialConditions.hpp"

	using namespace dealii;

	template <int dim>
	double
	Electrons_Equilibrium<dim>::
	value(const Point<dim> &p,
		  const unsigned int component) const
	{
			// dim+1 components
		if(component < dim)
		{	
			// set the components of the current initially to zero
			return ZeroFunction<dim>(dim+1).value(p, component);
		}
		else //(component == dim)
		{
//			if(p[0] < 0.2)
				return 2; 
//			else 
//				return 0;
		}
	}

	template <int dim>
	double
	Holes_Equilibrium<dim>::
	value(const Point<dim> &p,
		  const unsigned int component) const
	{
		// dim+1 components
		if(component < dim)
		{	
			// set the components of the current initially to zero
			return ZeroFunction<dim>(dim+1).value(p, component);
		}
		else //(component == dim)
		{
//			if(p[0] > 0.2)'
//				return 2; 
//			else
				return 0;
		}
	}

	template <int dim>
	double
	Reductants_Equilibrium<dim>::
	value(const Point<dim> &p,
		  const unsigned int component) const
	{
			// dim+1 components
		if(component < dim)
		{	
			// set the components of the current initially to zero
			return ZeroFunction<dim>(dim+1).value(p, component);
		}
		else //(component == dim)
		{
				return 30; 
		}
	}

	template <int dim>
	double
	Oxidants_Equilibrium<dim>::
	value(const Point<dim> &p,
		  const unsigned int component) const
	{
			// dim+1 components
		if(component < dim)
		{	
			// set the components of the current initially to zero
			return ZeroFunction<dim>(dim+1).value(p, component);
		}
		else //(component == dim)
		{
				return 29; 
		}
	}

