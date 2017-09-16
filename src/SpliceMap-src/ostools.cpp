

#include "ostools.h"



ostools::ostools(int os)
{
	curr_os = os;
}

bool ostools::program_exists(string name)
{
	int out = 1;
	if (curr_os == NIX) {
		string cmd = "which " + name + " > /dev/null 2> /dev/null";
		out = system(cmd.c_str());
	}else if (curr_os == WIN) {
		
	}
	
	return (out == 0);
}

bool ostools::mkdir(string dir_name)
{
	int out = 1;
	if (curr_os == NIX) {
		string cmd = "mkdir " + dir_name + " 2> /dev/null";
		out = system(cmd.c_str());
	}else if (curr_os == WIN) {
		
	}
	
	return (out == 0);
}


vector<string> ostools::list_dir(string path)
{
	vector<string> result;
	string cmd;
	string temp;
	if (curr_os == NIX) {
		
		FILE *fp;
		char line[1024];
		
		cmd = "ls " + path;
		
		fp = popen(cmd.c_str(), "r");
		
		if (fp == NULL) {
			cout << "ERROR: failed to open pipe" << endl;
			
			exit(2);
		}else {
			while ( fgets( line, 1024, fp))
			{
				temp = line;
				trim2(temp);
				if (temp.length() > 0) {
					result.push_back(temp);
				}
			}
		}
		pclose(fp);

	}else if (curr_os == WIN) {
		
	}
	
	return result;	
}

vector<string> ostools::list_bowtie(string package_path,string path)
{
	vector<string> result;
	string cmd;
	string temp;
	if (curr_os == NIX) {
		
		FILE *fp;
		char line[1024];
		
		cmd = package_path+"bowtie-inspect -n " + path;
		
		fp = popen(cmd.c_str(), "r");
		
		if (fp == NULL) {
			cout << "ERROR: failed to open pipe" << endl;
			
			exit(2);
		}else {
			while ( fgets( line, 1024, fp))
			{
				temp = line;
				trim2(temp);
				if (temp.length() > 0) {
					result.push_back(temp);
				}
			}
		}
		pclose(fp);
		
	}else if (curr_os == WIN) {
		
	}
	
	return result;	
}



