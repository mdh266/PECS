#ifndef _BOUNDARY_VALUES_H___
#define _BOUNDARY_VALUES_H___

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

///////////////////////////////////////////////////////////////////////////////
// Poisson Boundary Functions
///////////////////////////////////////////////////////////////////////////////

/** \brief The built in bias of the semiconductor \$f\Phi_{\text{bi}}\f$.*/
template <int dim>
class Built_In_Bias : public dealii::Function<dim>
{
	public:
		/** \brief Default constructor.*/
		Built_In_Bias() : dealii::Function<dim>()
		{}
	
		void set_value(const double & bias_value);

		/** \brief Returns value of \$f\Phi_{\text{bi}}\f$ at point p.*/ 
		virtual double value(const dealii::Point<dim> &p, 
				     const unsigned int component = 0 ) const;

	private:
		double built_in_bias;
};

/** \brief The schottky bias \$f\Phi_{\text{sch}}\f$.*/
template <int dim>
class Schottky_Bias : public dealii::Function<dim>
{
	public:
		/** \brief Default constructor.*/
		Schottky_Bias() : dealii::Function<dim>()
		{}

		void set_location(const double & bias_location);

		void set_value(const double & bias_value);
		
		/** \brief Returns value of \$f\Phi_{\text{sch.}}\f$ at point p.*/ 
		virtual double value(const dealii::Point<dim> &p, 
				     const unsigned int component = 0 ) const;

	private:
		double Schottky_bias;
		double Schottky_location;

};



/** \brief The applied bias \$f\Phi_{\text{app}}\f$.*/
template <int dim>
class Applied_Bias : public dealii::Function<dim>
{
	public:
		/** \brief Default constructor.*/
		Applied_Bias() : dealii::Function<dim>()
		{}

		void set_value(const double & bias_value);
		
		/** \brief Returns value of \$f\Phi_{\text{app.}}\f$ at point p.*/ 
		virtual double value(const dealii::Point<dim> &p, 
				     const unsigned int component = 0 ) const;

	private:
		double applied_bias;
};



/** \brief The bulk bias of the electrolyte \$f\Phi^{\infty}\f$.*/
template <int dim>
class Bulk_Bias : public dealii::Function<dim>
{
	public:
		/** \brief Default constructor.*/
		Bulk_Bias() : dealii::Function<dim>()
		{}
		
		/** \brief Returns value of \$f\Phi^{\infty}}\f$ at point p.*/ 
		virtual double value(const dealii::Point<dim> &p, 
				     const unsigned int component = 0 ) const;
};


#endif
