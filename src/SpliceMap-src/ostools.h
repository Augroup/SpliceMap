

#include "SpliceMap_utils.h"



using namespace std;

#ifndef _OSTOOLS_H
#define  _OSTOOLS_H

class ostools  {
public:
	ostools(int os);
	
	
	bool program_exists(string name);
	bool mkdir(string dir_name);
	vector<string> list_dir(string path);
	
	vector<string> list_bowtie(string package_path,string path);
	
	const static int NIX = 1;    
	const static int WIN = 2;    
private:
	int curr_os;
};

#endif

