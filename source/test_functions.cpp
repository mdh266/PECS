#ifndef _TEST_FUNCTIONS_H__
#define _TEST_FUNCTIONS_H__

#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
#include <cmath>


/*---------------------------------------------------------------------*/
/* POISSON TEST							 */
/*---------------------------------------------------------------------*/
namespace test_Poisson
{
	using namespace dealii;

	template <int dim>
	class RightHandSide : public Function<dim>
	{
		public:
			RightHandSide() : Function<dim>(1)
			{}
			
			virtual double value(const Point<dim> &p, 
					    const unsigned int component = 0 ) const;
	};

	template <int dim>
	class DirichletBoundaryValues : public Function<dim>
	{
		public:
			DirichletBoundaryValues() : Function<dim>(1)
			{}
			
			virtual double value(const Point<dim> &p, 
					     const unsigned int component = 0 ) const;
	};

	template <int dim>
	double 
	RightHandSide<dim>::
	value(const Point<dim> &p, 
			 const unsigned int ) const
	{
		const double x = p[0];
		const double y = p[1];
		return 4*M_PI*M_PI*(cos(2*M_PI*y) - sin(2*M_PI*x));

	}

	template <int dim>
	double
	DirichletBoundaryValues<dim>::
	value(const Point<dim> &p, 
				const unsigned int ) const
	{
		const double x = p[0];
		const double y = p[1];
		return cos(2*M_PI*y) -sin(2*M_PI*x) - x;
	}

	template<int dim>
	class TrueSolution : public Function<dim>
	{
		public:
			TrueSolution() : Function<dim>(dim+1)
			{}

			virtual void vector_value(const Point<dim> & p,
						  Vector<double> &valuess) const;
	};

	template <int dim>
	void 
	TrueSolution<dim>::
	vector_value(const Point<dim> &p,
				 Vector<double> &values) const
	{
		Assert(values.size() == dim+1,	
					 ExcDimensionMismatch(values.size(), dim+1) );

		double x = p[0];
		double y = p[1];

		//gradient values
		values(0) = 1 + 2*M_PI*cos(2*M_PI*x);
		values(1) = 2*M_PI*sin(2*M_PI*y);

		// function values
		values(2) = cos(2*M_PI*y) - sin(2*M_PI*x) - x;
	}

} // end test_Poisson namespace

/*---------------------------------------------------------------------*/
/*		LDG-IMEX TEST		  			       */
/*---------------------------------------------------------------------*/

namespace test_LDG_IMEX
{
	using namespace dealii;
	
	template<int dim>
	class RightHandSide : public Function<dim>
	{
		public:
		RightHandSide() : Function<dim>(1)
		{}

		virtual double value(const Point<dim>	& p,
				     const unsigned int component=0) const;
	};

	template<int dim>
	class DirichletBC : public Function<dim>
	{
		public:
		DirichletBC() : Function<dim>(1)
		{}

		virtual double value(const Point<dim>	& p,
				     const unsigned int component=0) const;
	};


	template<int dim>
	class InterfaceFunction : public Function<dim>
	{
		public:
		InterfaceFunction() : Function<dim>(1)
		{}

		virtual double value(const Point<dim>	& p,
				     const unsigned int component=0) const;
	};

	template<int dim>
	class TrueSolution : public Function<dim>
	{
		public:
			TrueSolution() : Function<dim>(dim+1)
			{}

			virtual void vector_value(const Point<dim> & p,
						  Vector<double> &valuess) const;
	};

	template <int dim>
	double
	RightHandSide<dim>::
	value(const Point<dim> &p,
				const unsigned int) const
	{
		double x = p[0];
		double y = p[1];
	
		const double time = this->get_time();
		const double two_pi = 2*M_PI;

		return -exp(-time)
			+ two_pi * two_pi * cos(two_pi * x) 
			+ two_pi * two_pi * cos(two_pi * y)
			+ two_pi * sin(two_pi * x);
	}

	template <int dim>
	double
	DirichletBC<dim>::
	value(const Point<dim> &p,
				const unsigned int) const
	{
		double x = p[0];
		double y = p[1];
	
		const double time = this->get_time();
		const double two_pi = 2*M_PI;

		return exp(-time) 			
			+ cos(two_pi * x) 
			+ cos(two_pi * y);
	}

	template <int dim>
	double
	InterfaceFunction<dim>::
	value(const Point<dim> &p,
				const unsigned int) const
	{
		double y = p[1];
	
		const double time = this->get_time();
		const double two_pi = 2*M_PI;

		return -exp(-time) -cos(two_pi *y) - 1;
	}

	template <int dim>
	void 
	TrueSolution<dim>::
	vector_value(const Point<dim> &p,
		     Vector<double> &values) const
	{
		Assert(values.size() == dim+1,	
			ExcDimensionMismatch(values.size(), dim+1) );

		double x = p[0];
		double y = p[1];

		const double time = this->get_time();
		const double two_pi = 2*M_PI;

		//function values
		values(2) = exp(-time)  + cos(two_pi*x) + cos(two_pi*y);

		//current values
		values(0) =	two_pi * sin(two_pi*x) - values(2);

		values(1) = two_pi * sin(two_pi*y);
	}
} // LDG_IMEX


/*---------------------------------------------------------------------*/
/*		DD-Poisson TEST					       */
/*---------------------------------------------------------------------*/

namespace test_DD_Poisson
{
	using namespace dealii;

	template<int dim>
	class DD_RightHandSide : public Function<dim>
	{
		public:
		DD_RightHandSide() : Function<dim>(1)
		{}

		virtual double value(const Point<dim>	& p,
				     const unsigned int component=0) const;
	};

	template<int dim>
	class Poisson_RightHandSide : public Function<dim>
	{
		public:
		Poisson_RightHandSide() : Function<dim>(1)
		{}

		virtual double value(const Point<dim>	& p,
				     const unsigned int component=0) const;
	};

	template<int dim>
	class DD_DirichletBC : public Function<dim>
	{
		public:
		DD_DirichletBC() : Function<dim>(1)
		{}

		virtual double value(const Point<dim>	& p,
				     const unsigned int component=0) const;
	};

	template<int dim>
	class DD_TrueSolution : public Function<dim>
	{
		public:
			DD_TrueSolution() : Function<dim>(dim+1)
			{}

			virtual void vector_value(const Point<dim> & p,
						  Vector<double> &valuess) const;
	};

	template <int dim>
	double
	DD_RightHandSide<dim>::
	value(const Point<dim> &p,
		  const unsigned int) const
	{
		double x = p[0];
		double y = p[1];
	
		const double time = this->get_time();
		const double two_pi = 2*M_PI;

		const double u_at_time =  exp(-time) + cos(two_pi*x) + cos(two_pi*y);

		const double div_E_u = two_pi * two_pi *
				( cos(two_pi*y) - sin(two_pi*x) ) * u_at_time  
				- two_pi * (two_pi * cos(two_pi*x) + 1) * sin(two_pi*x)
				- two_pi * two_pi * sin(two_pi*y) * sin(two_pi*y);

		return -exp(-time)
			 + two_pi * two_pi * cos(two_pi * x) 
			 + two_pi * two_pi * cos(two_pi * y)
			 - div_E_u;
	}

	template <int dim>
	double
	Poisson_RightHandSide<dim>::
	value(const Point<dim> &p,
		  const unsigned int) const
	{	
		double x = p[0];
		double y = p[1];
	
		const double time = this->get_time();
		const double two_pi = 2*M_PI;

		return (4*M_PI*M_PI*(cos(2*M_PI*y) - sin(2*M_PI*x)) 
			+ exp(-time) + cos(two_pi*x) + cos(two_pi*y));
	}

	template <int dim>
	double
	DD_DirichletBC<dim>::
	value(const Point<dim> &p,
		  const unsigned int) const
	{
		double x = p[0];
		double y = p[1];
	
		const double time = this->get_time();
		const double two_pi = 2*M_PI;

		return exp(-time) 			
			+ cos(two_pi * x) 
			+ cos(two_pi * y);
	}

	template <int dim>
	void 
	DD_TrueSolution<dim>::
	vector_value(const Point<dim> &p,
				 Vector<double> &values) const
	{
		Assert(values.size() == dim+1,	
			 ExcDimensionMismatch(values.size(), dim+1) );

		double x = p[0];
		double y = p[1];

		const double time = this->get_time();
		const double two_pi = 2*M_PI;

		//function values
		values(2) = exp(-time)  + cos(two_pi*x) + cos(two_pi*y);


		// q = -nabla u - Eu
		//current values
		values(0) =	two_pi * sin(two_pi*x) - 
				(two_pi * cos(two_pi*x) + 1) * (values(2));
	
		values(1) = two_pi * sin(two_pi*y) -
				two_pi * sin(two_pi*y) * values(2);
	}

} // test_DD_Poisson


/*---------------------------------------------------------------------*/
/*	INTERFACE PROBLEM  TEST					       */
/*---------------------------------------------------------------------*/


namespace test_interface_problem
{
	using namespace dealii;

	template<int dim>
	class InitialConditions : public Function<dim>
	{
		public:
			InitialConditions() : Function<dim>(dim+1)
			{}

			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;
	};
	
	template<int dim>
	class RightHandSide : public Function<dim>
	{
		public:
		RightHandSide() : Function<dim>(1)
		{}

		virtual double value(const Point<dim>	& p,
				     const unsigned int component=0) const;
	};

	template<int dim>
	class DirichletBC : public Function<dim>
	{
		public:
		DirichletBC() : Function<dim>(1)
		{}

		virtual double value(const Point<dim>	& p,
				     const unsigned int component=0) const;
	};


	template<int dim>
	class InterfaceFunction : public Function<dim>
	{
		public:
		InterfaceFunction() : Function<dim>(1)
		{}

		virtual double value(const Point<dim>	& p,
				    const unsigned int component=0) const;
	};

	template<int dim>
	class TrueSolution : public Function<dim>
	{
		public:
			TrueSolution() : Function<dim>(dim+1)
			{}

			virtual void vector_value(const Point<dim> & p,
						  Vector<double> &valuess) const;
	};

	template <int dim>
	double
	InitialConditions<dim>::
	value(const Point<dim> &p,
		  const unsigned int component) const
	{
		double x = p[0];
		double y = p[1];
	
		// dim+1 components
		if(component < dim)
		{	
			// set the components of the current initially to zero
			return ZeroFunction<dim>(dim+1).value(p, component);
		}
		else //(component == dim)
		{
			return 1 + cos(2*M_PI*x) * cos(2*M_PI*y);
		}
	}

	template <int dim>
	double
	RightHandSide<dim>::
	value(const Point<dim> &p,
		  const unsigned int) const
	{
		double x = p[0];
		double y = p[1];
	
		const double time = this->get_time();
		const double two_pi = 2*M_PI;
		
		return -exp(-time)
			 + two_pi * two_pi * cos(two_pi * x) 
			 + two_pi * two_pi * cos(two_pi * y);
	}

	template <int dim>
	double
	DirichletBC<dim>::
	value(const Point<dim> &p,
		  const unsigned int) const
	{
		double x = p[0];
		double y = p[1];
	
		const double time = this->get_time();
		const double two_pi = 2*M_PI;

		return exp(-time) 			
			 + cos(two_pi * x) 
			 + cos(two_pi * y);
	}

	template <int dim>
	double
	InterfaceFunction<dim>::
	value(const Point<dim> &p,
		  const unsigned int) const
	{
		double y = p[1];
		double x = p[0];

		const double time = this->get_time();
		const double two_pi = 2*M_PI;
//		const double elec_field_term =  (exp(-time) + cos(two_pi*x) + sin(two_pi*y)); 
	
		double value = 
			exp(-time) + cos(two_pi *y) + cos(two_pi*x);
			// + elec_field_term;

		return value * value;
	}



	template <int dim>
	void 
	TrueSolution<dim>::
	vector_value(const Point<dim> &p,
				 Vector<double> &values) const
	{
		Assert(values.size() == dim+1,	
			 ExcDimensionMismatch(values.size(), dim+1) );

		double x = p[0];
		double y = p[1];

		const double time = this->get_time();
		const double two_pi = 2*M_PI;

		//current values
		values(0) =	two_pi * sin(two_pi*x); 
		values(1) = two_pi * sin(two_pi*y);
		//function values
		values(2) = exp(-time)  + cos(two_pi*x) + cos(two_pi*y);
	}
}


#endif
