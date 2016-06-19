#include "../SolarCell.cpp"

 

int main()
{
	try{	
			using namespace dealii;
	
			deallog.depth_console(0);
			int degree   		= 1;

			ParameterHandler									prm;
			ParameterSpace::ParameterReader		param(prm);
			param.read_parameters("input_file.prm");

			// Run over number of refinements and calculate
			// L2 errors for LDG on Poissons equation and 
			// then print the run time performance of each the subroutines
			// finally a table listing the convergence values;	

			::MultithreadInfo::set_thread_limit();

			SOLARCELL::SolarCellProblem<2> 	DeviceSimulation(degree,prm);
			DeviceSimulation.run_full_system();
			
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
