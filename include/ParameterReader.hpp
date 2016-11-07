#ifndef _PARAMETERREADER_H__
#define _PARAMETERREADER_H__

#include<deal.II/base/parameter_handler.h>



/** \brief This namespace holds the paramaters Parameters and input file reader
*	ParameterReader.
 */
namespace ParameterSpace
{
	using namespace dealii;

	/// \brief Reads in the input file, if one isn't provided a default is created.
	/** This object is instantiated in the main method and used as follows,
	*
	*		ParameterHandler			prm;
	*		ParameterSpace::ParameterReader		param(prm);
	*		param.read_parameters("input_file.prm");
	*
	*		SOLARCELL::SolarCellProblem<2> 	DeviceSimulation(degree,prm);
	*
	*/
	class ParameterReader
	{
		public:
			/** \brief Class constructor 
			*	The class will need a ParameterHander object to be instantiated before
			*	it is instantiated. 
			*/
			ParameterReader(ParameterHandler & param_handler);
	
			/** \brief Reads input file, if the file doesn't exist a default one is created.
			*	The information is then stored in the ParameterHandler before being
			*	called stored in parse_and_scaled_parameters() when called 
			*	from the SolarCell constructor.
			*/
			void read_parameters(const std::string);
	
			void read_test_parameters(const std::string parameter_file);
		private:
			/** \brief States what parameters are to be read in and what their default 
			* values are. 
			*/
			void	declare_parameters();

			/** \brief States which test parameters are to be read in and what their default
 			* value are.
 			*/
			void	declare_test_parameters();

			/** deal.ii object which temporarily stores the parameters to be read in.*/
			ParameterHandler & prm;
	};
}


/// @author Michael Harmon

#endif	
