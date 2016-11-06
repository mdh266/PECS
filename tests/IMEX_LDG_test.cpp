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


			unsigned int n_min = 2;
			unsigned int n_max = 6;

			// Run over number of refinements and calculate
			// L2 errors for LDG on Poissons equation and 
			// then print the run time performance of each the subroutines
			// finally a table listing the convergence values;	
			std::cout << "\n\n\nTRANSIENT ON DRIFT-DIFFUSION EQUATION" << std::endl;
			ConvergenceTable		LDG_tran_table;

			for(unsigned int n_refine = n_min; n_refine < n_max; n_refine++)
			{
				SOLARCELL::SolarCellProblem<2> 	DeviceSimulation(degree,prm);
				DeviceSimulation.test_transient(n_refine,
												LDG_tran_table);
			}
		
			// get the convergence rates for carrier_1 and carrier_2
			LDG_tran_table.evaluate_convergence_rates("u", 
													ConvergenceTable::reduction_rate_log2);
	
			LDG_tran_table.evaluate_convergence_rates("J", 
													ConvergenceTable::reduction_rate_log2);
	
			// set to 3 significant digits as output
			LDG_tran_table.set_precision("u", 3);
			LDG_tran_table.set_precision("J", 3);

			// use scinetific notation for output
			LDG_tran_table.set_scientific("u", true);
			LDG_tran_table.set_scientific("J", true);

			std::ostringstream oss; 
			oss << degree;

			std::cout << std::endl;
			std::cout << "LDG L2 Errors For Drift-Diffusion" << std::endl;
			std::cout << std::endl;
			LDG_tran_table.write_text(std::cout);

			std::string DD_LDG_file = "DD_LDG_Deg_";
			DD_LDG_file += oss.str();
			DD_LDG_file += ".tex";

			std::ofstream table_file(DD_LDG_file.c_str());
			LDG_tran_table.write_text(table_file);
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
