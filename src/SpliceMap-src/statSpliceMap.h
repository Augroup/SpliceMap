

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>

using namespace std; 

int file_index(string filename);
void count_read(ifstream &junction, bool debug, int index);
void count_indep_read(ifstream &junction,bool debug, int lowerbound, int index);
void print_idep_read(ifstream &junction, int lowerbound, int index);


