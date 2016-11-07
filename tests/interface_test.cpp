#include "../source/SolarCell.cpp"
#include <string>
#include <sstream>
#include <fstream>

int main()
{
	try{	
			using namespace dealii;
	
			deallog.depth_console(0);
			int degree   		= 2;

//			MultithreadInfo::set_thread_limit(1);

			ParameterHandler			prm;
			ParameterSpace::ParameterReader		param(prm);
			param.read_test_parameters("test_file.prm");


			unsigned int n_min = 1;
			unsigned int n_max = 5;
	
			std::cout << "\n\n\nTRANSIENT COUPLING" << std::endl;

			ConvergenceTable		LDG_coupling_table;
			for(unsigned int n_refine = n_min; n_refine < n_max; n_refine++)
			{
				SOLARCELL::SolarCellProblem<2> 	DeviceSimulation(degree,prm);
				DeviceSimulation.test_interface_coupling(n_refine,
								LDG_coupling_table);
			}

			// get the convergence rates for carrier_1 and carrier_2
			LDG_coupling_table.evaluate_convergence_rates("u1", 
						ConvergenceTable::reduction_rate_log2);
	
			LDG_coupling_table.evaluate_convergence_rates("J1", 
						ConvergenceTable::reduction_rate_log2);
	
			// set to 3 significant digits as output
			LDG_coupling_table.set_precision("u1", 3);
			LDG_coupling_table.set_precision("J1", 3);

			// use scinetific notation for output
			LDG_coupling_table.set_scientific("u1", true);
			LDG_coupling_table.set_scientific("J1", true);


			// get the convergence rates for carrier_1 and carrier_2
			LDG_coupling_table.evaluate_convergence_rates("u2", 
						ConvergenceTable::reduction_rate_log2);
	
			LDG_coupling_table.evaluate_convergence_rates("J2", 
						ConvergenceTable::reduction_rate_log2);
	
			// set to 3 significant digits as output
			LDG_coupling_table.set_precision("u2", 3);
			LDG_coupling_table.set_precision("J2", 3);

			// use scinetific notation for output
			LDG_coupling_table.set_scientific("u2", true);
			LDG_coupling_table.set_scientific("J2", true);

			std::ostringstream oss; 
			oss << degree;

			std::cout << std::endl;
			std::cout << "LDG L2 Errors For Interface Coupling" << std::endl;
			std::cout << std::endl;
			LDG_coupling_table.write_text(std::cout);
			std::cout << std::endl;

			std::string file = "Interface_LDG_Deg_";
			file += oss.str();
			file += ".tex";

			std::ofstream table_file(file.c_str());
			LDG_coupling_table.write_text(table_file);
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
