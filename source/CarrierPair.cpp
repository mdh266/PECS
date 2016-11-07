#include "../include/CarrierPair.hpp"
#include "Carrier.cpp"

namespace ChargeCarrierSpace
{
	using namespace dealii;

	template<int dim>
	CarrierPair<dim>::
	CarrierPair() 
		:
		carrier_1(),
		carrier_2()
	{}

	template<int dim>
	CarrierPair<dim>::
	~CarrierPair()
	{

	}
	
	template<int dim>
	void
	CarrierPair<dim>::
	setup_dofs(const FESystem<dim>		& fe,
		   DoFHandler<dim>		& dof_handler)
	{
		// distribute the dofs for the electron hole pair
		dof_handler.distribute_dofs(fe);
		
		// Renumber dofs for [ Current, Density ]^{T} set up
		DoFRenumbering::component_wise(dof_handler);
		unsigned int n_dofs = dof_handler.n_dofs();

		DynamicSparsityPattern carrier_system_dsp(n_dofs,n_dofs);

		// allocate memory for [0 , 0 ; 0 , M ]
		DynamicSparsityPattern carrier_mass_dsp(n_dofs,n_dofs);

 	 	// create the actual sparsity pattern and allocate memory for the matrices
		
		// carrier_1s:
  		DoFTools::make_flux_sparsity_pattern (dof_handler, carrier_system_dsp);
  		system_sparsity_pattern.copy_from(carrier_system_dsp);
 	
		carrier_1.system_matrix.reinit(system_sparsity_pattern);
		carrier_2.system_matrix.reinit(system_sparsity_pattern);

	  	DoFTools::make_sparsity_pattern (dof_handler, carrier_mass_dsp);
 	 	mass_sparsity_pattern.copy_from(carrier_mass_dsp);
	  	mass_matrix.reinit (mass_sparsity_pattern);

	  	// allocate memory for carrier_solutions
	  	carrier_1.solution.reinit (n_dofs); // [vector field, density]
	  	carrier_2.solution.reinit (n_dofs); // [vector field, density]
	
		// memeory for RHS
 	 	carrier_1.system_rhs.reinit (n_dofs); // [vector field, density]
 	 	carrier_2.system_rhs.reinit (n_dofs); // [vector field, density]
	
		constraints.clear();
		constraints.close();
	}

	template<int dim>
	void
	CarrierPair<dim>::
	print_info(DoFHandler<dim>  & dof_handler)
	{	
		std::vector<types::global_dof_index> dofs_per_component(dim+1);
		DoFTools::count_dofs_per_component(dof_handler, dofs_per_component);

		// get number of dofs in vector field components and of density
		// in each component/dimension of vector field has same number of dofs
		const unsigned int n_current = dim * dofs_per_component[0];
//		std::cout << "n_current = " << n_current << std::endl;	
	
		const unsigned int n_density = dofs_per_component[1]; 
//		std::cout << "n_density = " << n_density << std::endl;	
											 
		std::cout << "Number of DOFS carrier: "
				  << 2 * dof_handler.n_dofs()
				  << " = 2 x ( " << dof_handler.n_dofs() << " ) "
				  << " = 2 x (" << n_current << " + " << n_density << ")"
				  << std::endl;
	}

	template<int dim>
	void
	CarrierPair<dim>::
	print_dofs()
	{
		std::string dof_name_1 = carrier_1.name;
		dof_name_1 += ".dofs";
		std::string dof_name_2 = carrier_2.name;
		dof_name_2 += ".dofs";
		std::ofstream prt(dof_name_1.c_str());
		carrier_1.solution.block_write(prt);
		prt.close();
		prt.open(dof_name_2.c_str());
		carrier_2.solution.block_write(prt);
		prt.close();
	}
	
	template<int dim>
	void
	CarrierPair<dim>::
	read_dofs()
	{
		std::string dof_name_1 = carrier_1.name;
		dof_name_1 += ".dofs";
		std::string dof_name_2 = carrier_2.name;
		dof_name_2 += ".dofs";
		std::ifstream reader(dof_name_1.c_str());
		carrier_1.solution.block_read(reader);
		reader.close();
		reader.open(dof_name_2.c_str());
		carrier_2.solution.block_read(reader);
		reader.close();
	}

	template<int dim>
	void 
	CarrierPair<dim>::
	set_semiconductor_for_testing(const ParameterSpace::Parameters & sim_params)
	{
		penalty = 1.0;
		carrier_1.scaled_mobility = sim_params.scaled_electron_mobility;
		carrier_2.scaled_mobility = sim_params.scaled_hole_mobility;
	}
	
	template<int dim>
	void 
	CarrierPair<dim>::
	set_electrolyte_for_testing(const ParameterSpace::Parameters & sim_params)
	{
		penalty = 1.0;
		carrier_1.scaled_mobility = sim_params.scaled_reductant_mobility;
		carrier_2.scaled_mobility = sim_params.scaled_oxidant_mobility;	
	}

	template<int dim>
	void
	CarrierPair<dim>::
	set_name(const std::string & str_name)
	{
		material_name = str_name;
	}

}

