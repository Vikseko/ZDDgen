# ZDDgen
To compile the application, you need the installed CUDD package.  
  
Example compile command:  
  
//g++ -o ZDDgen_direct ZDDgen_direct.cpp -I /home/path_to_CUDD_directory -I /home/path_to_CUDD_directory/util -I /home/path_to_CUDD_directory/cudd /home/path_to_CUDD_directory/cudd/.libs/libcudd.a /home/path_to_CUDD_directory/dddmp/.libs/libdddmp.a
  
Run as follows:  
./app input_file_name  
The input file is CNF in DIMACS format.  
