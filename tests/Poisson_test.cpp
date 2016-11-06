#include "../source/SolarCell.cpp"
#include <string>
#include <sstream>
#include <fstream>

int main()
{
	try{	
			using namespace dealii;
	
			deallog.depth_console(0);
			int degree   		= 1;

//			MultithreadInfo::set_thread_limit(1);

			ParameterHandler					prm;
			ParameterSpace::ParameterReader		param(prm);
			param.read_test_parameters("test_file.prm");


			unsigned int n_min = 1;
			unsigned int n_max = 5;

			// Run over number of refinements and calculate
			// L2 errors for LDG on Poissons equation and 
			// then print the run time performance of each the subroutines
			// finally a table listing the convergence values;	
			ConvergenceTable		Mixed_table;
			ConvergenceTable		LDG_table;

			std::cout << "STEADY STATE ON POISSONS EQUATION" << std::endl;
			
			for(unsigned int n_refine = n_min; n_refine < n_max; n_refine++)
			{
				SOLARCELL::SolarCellProblem<2> 	DeviceSimulation(degree,prm);
				DeviceSimulation.test_steady_state(n_refine,
												Mixed_table,
												LDG_table);
			}

			// get the convergence rates for potential and electric field
			Mixed_table.evaluate_convergence_rates("Phi", 
												ConvergenceTable::reduction_rate_log2);
	
			Mixed_table.evaluate_convergence_rates("D", 
												ConvergenceTable::reduction_rate_log2);
	
			// set to 3 significant digits as output
			Mixed_table.set_precision("Phi", 3);
			Mixed_table.set_precision("D", 3);

			// use scinetific notation for output
			Mixed_table.set_scientific("Phi", true);
			Mixed_table.set_scientific("D", true);
			// get the convergence rates for carrier_1 and carrier_2
			LDG_table.evaluate_convergence_rates("u", 
												ConvergenceTable::reduction_rate_log2);
	
			LDG_table.evaluate_convergence_rates("J", 
												ConvergenceTable::reduction_rate_log2);
	
			// set to 3 significant digits as output
			LDG_table.set_precision("u", 3);
			LDG_table.set_precision("J", 3);

			// use scinetific notation for output
			LDG_table.set_scientific("u", true);
			LDG_table.set_scientific("J", true);

			std::ostringstream oss; 
			oss << degree;

			std::cout << std::endl;
			std::cout << "Mixed Method L2 Errors For Poisson" << std::endl;
			std::cout << std::endl;
			Mixed_table.write_text(std::cout);

			std::cout << std::endl;
			std::cout << "LDG L2 Errors For Poisson" << std::endl;
			std::cout << std::endl;
			LDG_table.write_text(std::cout);

			std::string mixed_file = "Mixed_Deg_";
			mixed_file += oss.str();
			mixed_file += ".tex";
			
			std::ofstream table_file(mixed_file.c_str());
			Mixed_table.write_text(table_file);
			table_file.close();

			std::string LDG_file = "Poisson_LDG_Deg_";
			LDG_file += oss.str();
			LDG_file += ".tex";

			table_file.open(LDG_file.c_str());
			LDG_table.write_text(table_file);
			table_file.close();
		}
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }


		
	return 0;
}
