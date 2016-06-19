#include "../include/Carrier.hpp"


namespace ChargeCarrierSpace
{
	using namespace dealii;

	template<int dim>
	Carrier<dim>::
	Carrier()
	{}
	
	template<int dim>
	Carrier<dim>::
	~Carrier()
	{}

	template<int dim>
	void
	Carrier<dim>::
	set_name(const std::string & str_name)
	{
		name = str_name;
	}

	template<int dim>
	void
	Carrier<dim>::
	set_solver() 
	{
		solver.initialize(system_matrix);
	}

	template<int dim>
	void
	Carrier<dim>::
	solve() 
	{
		solver.vmult(solution, system_rhs);
	}
	
}
