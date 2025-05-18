***MAKE FILE***
1. To build the make file cmake cmake version 3.16.3 is used in linux environment.
	cmake -S. -B.
	make ---> makefile is in the build directory
	
***INPUT FILE***
2. All the input files are located inside the input folder. Change the path in DEM.cu accordingly for each case to read .par and .in file. 
	eg:- Default path is /input/$ CASE_NAME $
		 considering bounce.
         std::string path_par {"input/bounce.par"};
         std::string path_in  {"input/bounce.in"}; 
***OUTPUT FILE***	
3. All the output files are located in ./output/DEM/$ NAME SPECIFIED  IN THE INPUT PARAMETER FILE $
   i.e there are seperate folders foe each case.
   	
	eg:- Default path is /output/DEM/$ CASE_NAME $
		 considering bounce.
         std::string path_par {"input/bounce.par"};
         std::string path_in  {"input/bounce.in"}; 
         std::string out_path {"output/DEM/bounce/"};
         
         you can make it empty if you wish to generate all the outputs in the current location.
***EXECUTABLE***       
4.  DEM_SOLVER is my executable for the solver which is now in the build directory
**visualization**
paraview is used for the visualization purposes
