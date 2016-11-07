#include "../include/Poisson.hpp"


namespace Poisson
{
	using namespace dealii;

	template<int dim>
	PoissonData<dim>::
	PoissonData()
	{}

	template<int dim>
	PoissonData<dim>::
	~PoissonData()
	{}
	
	template<int dim>
	void
	PoissonData<dim>::
	setup_dofs(const FESystem<dim>	& fe,
		   DoFHandler<dim>	& dof_handler)
	{
		// distribute dofs
		dof_handler.distribute_dofs(fe);
	
		// redistribute dofs for [Elec Field, Potential]^T set up
		DoFRenumbering::component_wise(dof_handler);
	
		unsigned int n_dofs = dof_handler.n_dofs();

		// make hanging node constraints
		constraints.clear();
		DoFTools::make_hanging_node_constraints(dof_handler,											constraints);	
	
		// add constaints to Poisson_Neumann_constraints to constrain dofs
		// of electric field to be 0 along the Neumann boundary. 
		//
		// NOTE: The constraints must be added to the constraint matrix before
		// assembling Poisson's system matrix or right hand side.
		const FEValuesExtractors::Vector	ElectricField(0);
		ComponentMask	electric_field_mask	= fe.component_mask(ElectricField);
		
		DoFTools::make_zero_boundary_constraints(dof_handler,
							Neumann, // NEUMANN BOUNDARY INDICATOR
							constraints,
							electric_field_mask);
	
		constraints.close();

		// create temporary sparsity pattern
		DynamicSparsityPattern Poisson_dsp(n_dofs,n_dofs);

		DoFTools::make_sparsity_pattern(dof_handler,
						Poisson_dsp);

		constraints.condense(Poisson_dsp);
		
		sparsity_pattern.copy_from(Poisson_dsp);

		// allocate memory for the Poisson matrix 
		system_matrix.reinit(sparsity_pattern);

		// allocate memory for the Poisson solution vector 
		solution.reinit(n_dofs);

		// allocate memory for the Poisson RHS vector 
		system_rhs.reinit(n_dofs);
	}

	template<int dim>
	void
	PoissonData<dim>::
	print_info(DoFHandler<dim>	& dof_handler)
	{
		// Get the dnumber of dofs for each component
		std::vector<types::global_dof_index> dofs_per_component(dim+1);
		DoFTools::count_dofs_per_component(dof_handler, dofs_per_component);
		
		const unsigned int n_electric_field 	= dofs_per_component[0];
		const unsigned int n_potential		= dofs_per_component[dim];

		std::cout << "Number of DOFS Poisson: "
				  << dof_handler.n_dofs()
				  << " (" << n_electric_field 
				  << " + " << n_potential << ")"
				  << std::endl;
	}

	template<int dim>
	void
	PoissonData<dim>::
	set_solver()
	{
		solver.initialize(system_matrix);
	}

	template<int dim>
	void
	PoissonData<dim>::
	solve()
	{
		solver.vmult(solution, system_rhs);
		constraints.distribute(solution);
	}

}
