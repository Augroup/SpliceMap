


#include "runSpliceMap.h"


int main (int argc, char * const argv[]) {
    
	string genome_direc;  
	string bowtie_base_dir; 
	string ref_filename = "";
	char mapping_opt = '#'; 
	
	string temp_path;
	string out_path;
	vector<string> reads_filename[2];
	vector<string>::iterator readslist_it;
	string cfg_filename;
	string read_format;
	
	string temp;
	string line;
	string cmd;
	int system_status;
	bool cmd_out;
	
	ostools os(ostools::NIX);
	
	ofstream debug_out;

	vector<string> chr_file_list;
	
	pair<string,string> package_path_filename;
	pair<string,string> ref_path_filename;   
	
	map<string,reference_t> ref_map; 
	map<string,reference_t>::iterator ref_map_it;
	
	int num_pair = 0;
	int_fast32_t num_seed_mismatch = 1; 
	int_fast32_t num_read_mismatch = 2;  
	int_fast32_t max_clip_allowed = 40;  
	
	
	int fullread_length = -1;
	unsigned int head_clip_length = 0;
	bool print_sam = false; 
	bool cufflinks = false;
	
	string merge_opt = "";
	int num_threads = 2;
	string bowtie_opts = "";
	int max_intron = 400000;
	int min_intron = 20000;
	int min_extend_req = 0;
	string chromosome_wildcard = "chr*.fa";
	int num_chromosome_together = 1;
	int quality_format = 0; 
	int max_multi_hit = 10; 
	
	
	
	cout << "---== Welcome to SpliceMap " << VERSION << " ==---" << endl;
	cout << "Developed by Kin Fai Au and John C. Mu" << endl;
	cout << WEBSITE << endl;
	cout << "__________" << endl;
	
	struct timeval tv;
	struct timeval start_tv;
	struct timeval temp_start_tv;
	struct timeval map_temp_start_tv;
	
	gettimeofday(&start_tv, NULL);
	
	
	if (argc == 2) {
		
		temp = argv[0];
		package_path_filename = get_path_and_filename(temp);
		
		
		ifstream temp_infile;
		
		cfg_filename = argv[1];
		cout << "Loading configuration file... " << cfg_filename << endl;
		
		cfgfile *run_cfg = new cfgfile(cfg_filename);
		
		
		temp_path = run_cfg->getVal("temp_path");
		if (temp_path.length() == 0) {
			cout << "Warning: Temp directory location (temp_path) not specified, using default: temp/" << endl;
			temp_path = "temp/";
		}else {
			if (temp_path[temp_path.length()-1] != '/') {
				temp_path = temp_path + '/';
			}
		}
		
		

		
		
		out_path = run_cfg->getVal("out_path");
		if (out_path.length() == 0) {
			cout << "Warning: Output directory location (out_path) not specified, using default: output/" << endl;
			
			
			out_path = "output/";
		}else {
			if (out_path[out_path.length()-1] != '/') {
				out_path = out_path + '/';
			}
		}
		
		cmd_out = os.mkdir(out_path);
		if (cmd_out) {
			cout << "Created output directory" << endl;
			
		}else {
			cout << "output directory exists" << endl;
			
		}
		
		
		os.mkdir(out_path+debug_path);
		temp = out_path+debug_path+"run_debug.log";
		debug_out.open(temp.c_str(), ios::out);
		
		
		cmd_out = os.mkdir(temp_path);
		if (cmd_out) {
			cout << "Created temp directory" << endl;
			debug_out << "Created temp directory" << endl;
		}else {
			cout << "temp directory exists" << endl;
			debug_out << "temp directory exists" << endl;
		}
		
		
		genome_direc = run_cfg->getVal("genome_dir");
		if (genome_direc.length() == 0) {
			cout << "ERROR: The location of the genome files should be specified in the configuration file tag \"genome_dir\". " << endl;
			debug_out << "ERROR: The location of the genome files should be specified in the configuration file tag \"genome_dir\". " << endl;
			exit(2);
		}
		
		if (genome_direc[genome_direc.length()-1] != '/' && genome_direc[genome_direc.length()-1] != '\\') {
			genome_direc = genome_direc+'/';  
		}
		
		temp = run_cfg->getVal("chromosome_wildcard");
		if(temp.length() != 0){
			chromosome_wildcard = temp;
		}
		
		chr_file_list = os.list_dir(genome_direc+chromosome_wildcard);
		vector<string>::iterator chr_file_list_it = chr_file_list.begin();
		
		if (chr_file_list.size() == 0) {
			cout << "ERROR: I found no chromosomes, please check the genome_dir ("<< genome_direc <<") and the chromosome_wildcard (" << chromosome_wildcard <<")"<< endl;
			debug_out << "ERROR: I found no chromosomes, please check the genome_dir ("<< genome_direc <<") and the chromosome_wildcard (" << chromosome_wildcard <<")"<< endl;
			exit(2);
		}
		
		
		
		cout << "Scaning genome: " << genome_direc+chromosome_wildcard << endl;
		debug_out << "Scaning genome: " << genome_direc+chromosome_wildcard << endl;
		
		if (!read_reference_map(genome_direc+chromosome_wildcard, ref_map)) {
			cout << "ERROR: something wrong reading the genome files..." << endl;
			debug_out << "ERROR: something wrong reading the genome files..." << endl;
			exit(1);
		}
		
		if(ref_map.size() == 0){
			cout << "ERROR: No chromosomes found within the reference files" << endl;
			cout << "Please check that they are FASTA format" << endl;
			exit(1);
		}
		
		
		ofstream ref_list_file;
		temp = temp_path+reference_filename;
		ref_list_file.open(temp.c_str(), ios::out);
		
		ref_map_it = ref_map.begin();
		while (ref_map_it != ref_map.end()) {
			
			ref_list_file << ref_map_it->first << "\t" 
			<< ref_map_it->second.file_name << "\t" 
			<< ref_map_it->second.file_path << "\t" 
			<< ref_map_it->second.file_index_start << "\t"
			<< ref_map_it->second.file_index_end << "\n"; 
			
			ref_map_it++;
		}
		
		ref_list_file.close();
		ref_list_file.clear();
		
		
		cout << "List of chromosomes to be searched: " << endl;
		debug_out << "List of chromosomes to be searched: " << endl;
		
		
		ref_map_it = ref_map.begin();
		while (ref_map_it != ref_map.end()) {
			cout << ref_map_it->first << " | " << ref_map_it->second.file_path << ref_map_it->second.file_name << " | pos:" << ref_map_it->second.file_index_start << " - " << ref_map_it->second.file_index_end << "\n";
			debug_out << ref_map_it->first << " | " << ref_map_it->second.file_path << ref_map_it->second.file_name << " | pos:" << ref_map_it->second.file_index_start << " - " << ref_map_it->second.file_index_end << "\n";
			
			ref_map_it++;
		}
		cout << "Please check that these are correct... continuing in 7 s ...  < control-c >  to exit" << endl;
		cout << "If they are not correct please check chromosome_wildcard: " << chromosome_wildcard << endl;
		debug_out << flush;
		
		
		
		sleep(7);

		
		
		
		reads_filename[0] = run_cfg->getList("reads_list1");
		reads_filename[1] = run_cfg->getList("reads_list2");

		
		
		if(reads_filename[0].size() == 0){
			cout << "ERROR: There should be at least one reads file in the first list (reads_list1). " << endl;
			debug_out << "ERROR: There should be at least one reads file in the first list (reads_list1). " << endl;
			exit(2);
		}
		
		if(reads_filename[1].size() == 0){
			num_pair = 1;
		}else {
			num_pair = 2;
			
			if (reads_filename[0].size() != reads_filename[1].size()) {
				cout << "ERROR: Both read lists must be the same length. Please check reads_list1 and reads_list2" << endl;
				debug_out << "ERROR: Both read lists must be the same length. Please check reads_list1 and reads_list2" << endl;
				exit(2);
			}
		}
		
		
		temp = run_cfg->getVal("max_multi_hit");
		if (temp.length() > 0) {
			max_multi_hit = atoi(temp.c_str());
			
			if (max_multi_hit == 0) {
				cout << "ERROR: Zero max_multi_hit, we need at least one hit" << endl;
				debug_out << "ERROR: Zero max_multi_hit, we need at least one hit" << endl;
				exit(2);
			}
		}
		

		temp = run_cfg->getVal("full_read_length");
		if (temp.length() > 0) {
			fullread_length = atoi(temp.c_str());
			
			if (fullread_length == 0) {
				cout << "ERROR: Zero full_read_length" << endl;
				debug_out << "ERROR: Zero full_read_length" << endl;
				exit(2);
			}
		}
		
		
		temp = run_cfg->getVal("head_clip_length");
		if (temp.length() > 0) {
			head_clip_length = atoi(temp.c_str());
			if (head_clip_length < 0) {
				cout << "Warning: Negative head_clip_length (reseting to 0): " << head_clip_length << endl;
				debug_out << "Warning: Negative head_clip_length (reseting to 0): " << head_clip_length << endl;

				head_clip_length = 0;
			}
			if (fullread_length>0 && fullread_length-head_clip_length < 50) {
				cout << "ERROR: fullread_length-head_clip_length < 50, the minimum read length accepted by SpliceMap is 50" << endl;
				debug_out << "ERROR: fullread_length-head_clip_length < 50, the minimum read length accepted by SpliceMap is 50" << endl;
				
				exit(2);
			}
		}
		
		
		
		temp = run_cfg->getVal("seed_mismatch");
		if (temp.length() > 0) {
			num_seed_mismatch = atoi(temp.c_str());
		}
		
		
		temp = run_cfg->getVal("read_mismatch");
		if (temp.length() > 0) {
			num_read_mismatch = atoi(temp.c_str());
		}
		
		
		temp = run_cfg->getVal("max_clip_allowed");
		if (temp.length() > 0) {
			max_clip_allowed = atoi(temp.c_str());
		}
		
		
		ref_filename = run_cfg->getVal("annotations");
		ifstream temp_file;
		temp_file.open(ref_filename.c_str(), ios::in);
		if (ref_filename.length() > 0 && !temp_file.is_open()) {
			cout << "ERROR: I'm sorry, cannot open the annotation file, " << ref_filename << endl;
			debug_out << "ERROR: I'm sorry, cannot open the annotation file, " << ref_filename << endl;
			
			exit(2);
		}
		temp_file.close();
		temp_file.clear();
		
		ref_path_filename = get_path_and_filename(ref_filename);
		
		
		
		temp = run_cfg->getVal("mapper");
		if (temp.compare("bowtie") == 0) {
			mapping_opt = 'B';
		}else if(temp.compare("eland") == 0){
			mapping_opt = 'E';
			bool eland_exists = os.program_exists(package_path_filename.first+"eland_25");
			if (!eland_exists) {
				cout << "ERROR: Eland doesn't exist. Are you sure you have correctly installed Eland to your PATH?" << endl;
				cout << "Also make sure that it is called \"eland_25\" and compiled to handle 25bp reads" << endl;
				debug_out << "ERROR: Eland doesn't exist. Are you sure you have correctly installed Eland to your PATH?" << endl;
				debug_out << "Also make sure that it is called \"eland_25\" and compiled to handle 25bp reads" << endl;
				
				exit(2);
			}
			eland_exists = os.program_exists(package_path_filename.first+"squashGenome");
			if (!eland_exists) {
				cout << "ERROR: squashGenome doesn't exist. Are you sure you have correctly installed Eland to your PATH?" << endl;
				debug_out << "ERROR: squashGenome doesn't exist. Are you sure you have correctly installed Eland to your PATH?" << endl;
				
				exit(2);
			}
		}else if(temp.compare("seqmap") == 0){
			mapping_opt = 'S';
			bool seqmap_exists = os.program_exists(package_path_filename.first+"seqmap");
			if (!seqmap_exists) {
				cout << "ERROR: SeqMap doesn't exist. Are you sure you have correctly installed SeqMap to your PATH?" << endl;
				debug_out << "ERROR: squashGenome doesn't exist. Are you sure you have correctly installed Eland to your PATH?" << endl;
				
				exit(2);
			}
		}else if(temp.length() == 0){
			
		}else {
			cout << "ERROR: The choices for mapper are bowtie, eland or seqmap ... not  " << temp << endl;
			debug_out << "ERROR: The choices for mapper are bowtie, eland or seqmap ... not  " << temp << endl;
			
			exit(2);
		}
		
		if (temp.compare("bowtie") == 0) {
			
			bool bowtie_exists = os.program_exists(package_path_filename.first+"bowtie");
			if (!bowtie_exists) {
				cout << "ERROR: bowtie doesn't exist. Are you sure you have correctly installed bowtie to your PATH?" << endl;
				debug_out << "ERROR: bowtie doesn't exist. Are you sure you have correctly installed bowtie to your PATH?" << endl;
				
				exit(2);
			}
			bowtie_exists = os.program_exists(package_path_filename.first+"bowtie-build");
			if (!bowtie_exists) {
				cout << "ERROR: bowtie-build doesn't exist. Are you sure you have correctly installed bowtie-build to your PATH?" << endl;
				debug_out << "ERROR: bowtie-build doesn't exist. Are you sure you have correctly installed bowtie-build to your PATH?" << endl;
				
				exit(2);
			}
			
			
			bowtie_base_dir = run_cfg->getVal("bowtie_base_dir");
			
			bool bowtie_index_exist = true;
			
			if (bowtie_base_dir.length() == 0) {
				
				
				
				cout << "Warning: Bowtie base index not specified, building the index from the genome can take some time..." << endl;
				debug_out << "Warning: Bowtie base index not specified, building the index from the genome can take some time..." << endl;
				
				bowtie_index_exist = false;

			}else{

				
				temp = bowtie_base_dir+".1.ebwt";
				temp_infile.open(temp.c_str(), ios::in);
				if (!temp_infile.is_open()) {
					
					
					bowtie_index_exist = false;
				}
				temp_infile.close();
				temp_infile.clear();
				
				temp = bowtie_base_dir+".2.ebwt";
				temp_infile.open(temp.c_str(), ios::in);
				if (!temp_infile.is_open()) {
					
					
					bowtie_index_exist = false;
				}
				temp_infile.close();
				temp_infile.clear();
				
				temp = bowtie_base_dir+".3.ebwt";
				temp_infile.open(temp.c_str(), ios::in);
				if (!temp_infile.is_open()) {
					
					
					bowtie_index_exist = false;
				}
				temp_infile.close();
				temp_infile.clear();
				
				temp = bowtie_base_dir+".4.ebwt";
				temp_infile.open(temp.c_str(), ios::in);
				if (!temp_infile.is_open()) {
					
					
					bowtie_index_exist = false;
				}
				temp_infile.close();
				temp_infile.clear();
				
				temp = bowtie_base_dir+".rev.1.ebwt";
				temp_infile.open(temp.c_str(), ios::in);
				if (!temp_infile.is_open()) {
					
					
					bowtie_index_exist = false;
				}
				temp_infile.close();
				temp_infile.clear();
				
				temp = bowtie_base_dir+".rev.2.ebwt";
				temp_infile.open(temp.c_str(), ios::in);
				if (!temp_infile.is_open()) {
					
					
					bowtie_index_exist = false;
				}
				temp_infile.close();
				temp_infile.clear();
				
			}
			
			if (!bowtie_index_exist) {
				cout << "Warning: Bowtie index files not found, will attempt to generate them from the genome..." << endl;
				debug_out << "Warning: Bowtie index files not found, will attempt to generate them from the genome..." << endl;
				cout << "This could take a long time, I suggest you see the website for instruction on obtaining index files." << endl;
				
				cout << "continuing in 10 s ...  < control-c >  to exit" << endl;
				
				sleep(10);
				
				string chr_comma = "";
				chr_file_list_it = chr_file_list.begin();
				while (chr_file_list_it != (chr_file_list.end() - 1)) {
					chr_comma += (*chr_file_list_it) + ",";
					chr_file_list_it++;
				}
				chr_comma += (*chr_file_list_it);
				bowtie_base_dir = temp_path + "genome";
				cmd = package_path_filename.first+"bowtie-build " + chr_comma + " " + bowtie_base_dir;
				debug_out << "COMMAND: " << cmd << endl;
				system_status = system(cmd.c_str());
				
			}else {
				
				vector<string> bowtie_chr = os.list_bowtie(package_path_filename.first, bowtie_base_dir);
				
				sort(bowtie_chr.begin(), bowtie_chr.end());
				
				ref_map_it = ref_map.begin();
				while (ref_map_it != ref_map.end()) {
					
					if (!binary_search(bowtie_chr.begin(), bowtie_chr.end(), ref_map_it->first)) {
						cout << "ERROR: The following chromosome does not exist in the bowtie index" << endl;
						cout << ref_map_it->first << " | " << ref_map_it->second.file_path << ref_map_it->second.file_name << " | pos:" << ref_map_it->second.file_index_start << " - " << ref_map_it->second.file_index_end << "\n";
						cout << "Please make sure that all chromosomes are contained within the bowtie index" << endl;
						debug_out << "ERROR: The following chromosome does not exist in the bowtie index" << endl;
						debug_out << ref_map_it->first << " | " << ref_map_it->second.file_path << ref_map_it->second.file_name << " | pos:" << ref_map_it->second.file_index_start << " - " << ref_map_it->second.file_index_end << "\n";
						debug_out << "Please make sure that all chromosomes are contained within the bowtie index" << endl;
						
						exit(1);
					}
					
					ref_map_it++;
				}
			}

		}
		
		
		
		temp = run_cfg->getVal("sam_file");
		if (temp.compare("sam") == 0) {
			print_sam = true;
			merge_opt+=" -sam ";
		}else if(temp.compare("cuff") == 0){
			cufflinks = true;
			merge_opt+=" -cuff ";
		}else if(temp.length() == 0){
			print_sam = false;
			cufflinks = false;
		}else {
			cout << "ERROR: The choices for SAM output are \"sam\", \"cuff\" ... not  " << temp << endl;
			debug_out << "ERROR: The choices for SAM output are \"sam\", \"cuff\" ... not  " << temp << endl;
			
			exit(2);
		}

		
		temp = run_cfg->getVal("num_threads");
		if (temp.length() > 0) {
			num_threads = atoi(temp.c_str());
		}	
		
		
		temp = run_cfg->getVal("try_hard");
		if (temp.compare("yes") == 0) {
			bowtie_opts += " -y ";
		}else if (temp.compare("no") == 0) {
			
			
		}else if (temp.length() == 0) {
			
		}else {
			cout << "ERROR: The choices for try_hard are \"yes\" or \"no\" ... not  " << temp << endl;
			debug_out << "ERROR: The choices for try_hard are \"yes\" or \"no\" ... not  " << temp << endl;
			
			exit(2);
		}

		
		
		
		
		read_format = run_cfg->getVal("read_format");
		if (read_format.compare("FASTQ") == 0) {
			
		}else if (read_format.compare("FASTA") == 0) {
			
		}else if (read_format.compare("RAW") == 0) {
			
		}else {
			cout << "ERROR: The read format " << read_format << " is not supported. " << endl;
			cout << "Supported formats are FASTQ, FASTA and RAW" << endl;
			
			debug_out << "ERROR: The read format " << read_format << " is not supported. " << endl;
			debug_out << "Supported formats are FASTQ, FASTA and RAW" << endl;
			
			exit(2);
		}
		
		
		temp = run_cfg->getVal("quality_format");
		if (temp.compare("phred-33") == 0) {
			quality_format = 0;
		}else if (temp.compare("phred-64") == 0) {
			quality_format = 1;
		}else if (temp.compare("solexa") == 0) {
			quality_format = 2;
		}else if (temp.length() == 0){
			
		}else {
			cout << "ERROR: The quality format " << quality_format << " is not supported. " << endl;
			debug_out << "ERROR: The quality format " << quality_format << " is not supported. " << endl;
			
			cout << "Supported formats are phred-33(default), phred-64 and solexa" << endl;
			cout << "Note that \"Illumnina 1.3+\" format is the same as phred-64" << endl;
			exit(2);
		}

		
		int temp_intron_int;
		temp = run_cfg->getVal("max_intron");
		temp_intron_int = atoi(temp.c_str());
		if (temp_intron_int > 0) {
			max_intron = temp_intron_int;
		}
		
		
		temp = run_cfg->getVal("min_intron");
		temp_intron_int = atoi(temp.c_str());
		if (temp_intron_int > 50) {
			min_intron = temp_intron_int;
		}else if (temp_intron_int == 0) {
			
		}else {
			cout << "ERROR: \"min_intron\" option should be greater than 50. " << endl;
			cout << "Note that this is not the minimum search distance, SpliceMap will find small introns. :)" << endl;
			
			debug_out << "ERROR: \"min_intron\" option should be greater than 50. " << endl;
			debug_out << "Note that this is not the minimum search distance, SpliceMap will find small introns. :)" << endl;
			
			exit(2);
		}

		
		
		if (max_intron>0 && min_intron>0){
			if (min_intron > max_intron) {
				cout << "Warning: minimum intron parameter greater than maximum intron parameter... " << endl;
				debug_out << "Warning: minimum intron parameter greater than maximum intron parameter... " << endl;
				
				max_intron = min_intron;
			}
		}
		
		temp = run_cfg->getVal("min_extend_req");
		min_extend_req = atoi(temp.c_str()); 
		
		if (min_extend_req < 0 || min_extend_req > 25) {
			cout << "ERROR: valid values for min_extend_req are 0-25 not " << temp << endl;
			debug_out << "ERROR: valid values for min_extend_req are 0-25 not " << temp << endl;
			
			exit(2);
		}
		
		
		temp = run_cfg->getVal("num_chromosome_together");
		num_chromosome_together = atoi(temp.c_str()); 
		
		if (num_chromosome_together <= 0) {
			num_chromosome_together = 1;
		}
		
		delete run_cfg;
		
	} else {
		print_usage_and_exit();
	}
	
	
	
	if (cufflinks && print_sam) {
		cout << "ERROR: Please choose either -cuff or -sam. " << endl;
		debug_out << "ERROR: Please choose either -cuff or -sam. " << endl;
		
		cout << "-cuff will output a SAM file in a format compatible with cufflinks" << endl;
		cout << "-sam will output a proper SAM file" << endl;
		print_usage_and_exit();
	}
	
	
	
	

	

	
	
	
	
	
	
	cout << "__________" << endl;
	debug_out << "__________" << endl;
	
	cout << "Temp directory:   " << temp_path << endl;
	cout << "Output directory: " << out_path << endl;
	debug_out << "Temp directory:   " << temp_path << endl;
	debug_out << "Output directory: " << out_path << endl;
	
	cout << "Maximum number of multiple mapped reads allowed: " << max_multi_hit << endl;
	debug_out << "Maximum number of multiple mapped reads allowed: " << max_multi_hit << endl;
	
	cout << "Maximum number of mismatches allowed in 25-mer seed: " << num_seed_mismatch << endl;
	debug_out << "Maximum number of mismatches allowed in 25-mer seed: " << num_seed_mismatch << endl;
	
	cout << "Maximum number of mismatches allowed in full read: " << num_read_mismatch << endl;
	debug_out << "Maximum number of mismatches allowed in full read: " << num_read_mismatch << endl;
	
	
	cout << "Maximum number of bases SpliceMap is allowed to clip: " << max_clip_allowed << endl;
	debug_out << "Maximum number of bases SpliceMap is allowed to clip: " << max_clip_allowed << endl;
	
	cout << "Mapper used: ";
	debug_out << "Mapper used: ";
	if (mapping_opt == 'B') {
		cout << "bowtie" << endl;
		debug_out << "bowtie" << endl;
	}else if (mapping_opt == 'E') {
		cout << "eland" << endl;
		debug_out << "eland" << endl;
	}else {
		cout << "seqmap" << endl;
		debug_out << "seqmap" << endl;
	}

	cout << "(25th-percentile) intron size: " << min_intron << endl;
	cout << "(99th-percentile) intron size: " << max_intron << endl;
	debug_out << "(25th-percentile) intron size: " << min_intron << endl;
	debug_out << "(99th-percentile) intron size: " << max_intron << endl;
	
	
	
	cout << "Annotations path: " << ref_path_filename.first << " name: " << ref_path_filename.second << endl;
	cout << "Package path:     " << package_path_filename.first << " name: " << package_path_filename.second << endl;
	debug_out << "Annotations path: " << ref_path_filename.first << " name: " << ref_path_filename.second << endl;
	debug_out << "Package path:     " << package_path_filename.first << " name: " << package_path_filename.second << endl;
	
	cout << "Read format: " << read_format << endl;
	debug_out << "Read format: " << read_format << endl;
	
	if (read_format == "FASTQ") {
		cout << "Quality format: ";
		debug_out << "Quality format: ";
		if (quality_format == 0) {
			cout << "phred-33" << endl;
			debug_out << "phred-33" << endl;
		}else if (quality_format == 1){
			cout << "phred-64" << endl;
			debug_out << "phred-64" << endl;
		}else {
			cout << "solexa" << endl;
			debug_out << "solexa" << endl;
		}
	}
	
	cout << "Number of threads: " << num_threads << endl;
	debug_out << "Number of threads: " << num_threads << endl;
	
	cout << "Number of chromosomes to run together: " << num_chromosome_together << endl;
	debug_out << "Number of chromosomes to run together: " << num_chromosome_together << endl;
	
	
	if(print_sam){
		cout << "Will print regular SAM file" << endl;
		debug_out << "Will print regular SAM file" << endl;
	}
	if(cufflinks){
		cout << "Will print Cufflinks compatible SAM file" << endl;
		debug_out << "Will print Cufflinks compatible SAM file" << endl;
	}
	
	

	for(int j = 0;j<2;j++){
		cout << "Reads List " << (j+1) << ":" << endl;
		debug_out << "Reads List " << (j+1) << ":" << endl;
		for (unsigned int i = 0; i<reads_filename[j].size(); i++) {
			cout << reads_filename[j][i] << endl;
			debug_out << reads_filename[j][i] << endl;
		}
	}
	
	
	if (num_seed_mismatch<0 || num_seed_mismatch > 2){
		cout << "I'm sorry, SpliceMap " << VERSION << " supports a maximum of 2 mismatchs in each seed" << endl;
		debug_out << "I'm sorry, SpliceMap " << VERSION << " supports a maximum of 2 mismatchs in each seed" << endl;
		exit(1);
	}
	
	if (num_read_mismatch<0){
		cout << "ERROR: Read mismatch number negative: " << num_read_mismatch << endl;
		debug_out << "ERROR: Read mismatch number negative: " << num_read_mismatch << endl;
		exit(1);
	}
	
	if (max_clip_allowed<0){
		cout << "ERROR: Maximum clipping number negative: " << max_clip_allowed << endl;
		debug_out << "ERROR: Maximum clipping number negative: " << max_clip_allowed << endl;
		exit(1);
	}
	
	if (mapping_opt != '#') {
		
		
		cout << "Preparing the reads!..." << endl;
		debug_out << "Preparing the reads!..." << endl;
		
		
		
		
		
		
		
		
		
		
		if (read_format.compare("FASTQ") == 0) {
			for (int k = 0; k<num_pair; k++) {
				
				vector<string>::iterator reads_list_it = reads_filename[k].begin();
				while (reads_list_it != reads_filename[k].end()) {
					ifstream temp_infile;
					
					temp_infile.open(reads_list_it->c_str(), ios::in);
					
					if(!temp_infile.is_open()){
						cout << "ERROR: I'm sorry, cannot open the sequencer reads file, " << *reads_list_it << endl;
						cout << "Please check you have the correct path and are the correct working directory." << endl;
						debug_out << "ERROR: I'm sorry, cannot open the sequencer reads file, " << *reads_list_it << endl;
						debug_out << "Please check you have the correct path and are the correct working directory." << endl;
						exit(1);
					}
					
					getline(temp_infile, line);
					
					
					if (line[0] != '@') {
						cout << "Malformed FASTQ reads file: " << *reads_list_it << endl;
						cout << line << endl;
						cout << "Should begin with an @" << endl;
						debug_out << "Malformed FASTQ reads file: " << *reads_list_it << endl;
						debug_out << line << endl;
						debug_out << "Should begin with an @" << endl;
						exit(1);
					}
					
					
					temp_infile.close();
					temp_infile.clear();
					
					reads_list_it++;
				}
				
			}
		}else if (read_format.compare("FASTA") == 0) {
			for (int k = 0; k<num_pair; k++) {
				
				vector<string>::iterator reads_list_it = reads_filename[k].begin();
				while (reads_list_it != reads_filename[k].end()) {
					ifstream temp_infile;
					temp_infile.open(reads_list_it->c_str(), ios::in);
					
					if(!temp_infile.is_open()){
						cout << "ERROR: I'm sorry, cannot open the sequencer reads file, " << *reads_list_it << endl;
						cout << "Please check you have the correct path and are the correct working directory." << endl;
						debug_out << "ERROR: I'm sorry, cannot open the sequencer reads file, " << *reads_list_it << endl;
						debug_out << "Please check you have the correct path and are the correct working directory." << endl;
						exit(1);
					}
					
					getline(temp_infile, line);
					
					
					if (line[0] != '>') {
						cout << "Malformed FASTA reads file: " << *reads_list_it << endl;
						cout << line << endl;
						cout << "Should begin with a >" << endl;
						debug_out << "Malformed FASTA reads file: " << *reads_list_it << endl;
						debug_out << line << endl;
						debug_out << "Should begin with a >" << endl;
						exit(1);
					}
					
					
					temp_infile.close();
					temp_infile.clear();
					
					reads_list_it++;
				}
				
				
			}
		}else if (read_format.compare("RAW") == 0) {
			for (int k = 0; k<num_pair; k++) {
				
				vector<string>::iterator reads_list_it = reads_filename[k].begin();
				while (reads_list_it != reads_filename[k].end()) {
					ifstream temp_infile;
					temp_infile.open(reads_list_it->c_str(), ios::in);
					
					if(!temp_infile.is_open()){
						cout << "ERROR: I'm sorry, cannot open the sequencer reads file, " << *reads_list_it << endl;
						cout << "Please check you have the correct path and are in the correct working directory." << endl;
						debug_out << "ERROR: I'm sorry, cannot open the sequencer reads file, " << *reads_list_it << endl;
						debug_out << "Please check you have the correct path and are in the correct working directory." << endl;
						exit(1);
					}
					
					getline(temp_infile, line);
					
					
					if (!(line[0] == 'A' || line[0] == 'C' || line[0] == 'T' || line[0] == 'G'
						  || line[0] == 'a' || line[0] == 'c' || line[0] == 't' || line[0] == 'g'
						  || line[0] == 'N' || line[0] == 'n')) {
						cout << "Warning: Are you sure this is a RAW reads file? " << *reads_list_it << endl;
						cout << line << endl;
						debug_out << "Warning: Are you sure this is a RAW reads file? " << *reads_list_it << endl;
						debug_out << line << endl;
					}
					
					temp_infile.close();
					temp_infile.clear();
					reads_list_it++;
				}
				
			}
		}
		
		
		if(fullread_length > 0){
			cout << "Using at most bases " << (head_clip_length+1) << " to " << fullread_length << " (inclusive) of each read for mapping." << endl ;
			cout << "At most " << (fullread_length-head_clip_length) <<  " bases in total." << endl;
			debug_out << "Using at most bases " << (head_clip_length+1) << " to " << fullread_length << " (inclusive) of each read for mapping." << endl ;
			debug_out << "At most " << (fullread_length-head_clip_length) <<  " bases in total." << endl;
			
		}else { 
			cout << "Bases removed from front: " << head_clip_length << endl;
			debug_out << "Bases removed from front: " << head_clip_length << endl;
			cout << "Using as many bases as possible." << endl;
			debug_out << "Using as many bases as possible." << endl;
			
		}
		
		if (read_format.compare("FASTQ") == 0) {
			
			
			
			
			
			
			int read_count = 1;
			vector<string>::iterator reads_list_it[2];
			reads_list_it[0]= reads_filename[0].begin();
			reads_list_it[1]= reads_filename[1].begin();
			
			uint_fast32_t read_idx = 1;
			
			while (reads_list_it[0] != reads_filename[0].end()) { 
				ofstream temp_outfile[2];
				ofstream temp_name_outfile[2];
				ofstream temp_quality_outfile[2];
				
				ifstream temp_infile[2];
				
				for (int k = 0; k<num_pair; k++) {
					
					
					temp_infile[k].open(reads_list_it[k]->c_str(), ios::in);
					
					string new_filename = temp_path+"read_"+IntToStr(read_count)+"_"+IntToStr(k+1);
					temp_outfile[k].open(new_filename.c_str(), ios::out);
					
					*(reads_list_it[k]) = new_filename;
					
					
					temp = new_filename + ".names";
					temp_name_outfile[k].open(temp.c_str(), ios::out);
					temp = new_filename + ".quals";
					temp_quality_outfile[k].open(temp.c_str(), ios::out);
					
				}
				
				uint_fast32_t line_num = 1;
				
				bool line_good[2] = {false,false};
				string line_name[2] = {"",""};
				string line_seq[2] = {"",""};
				string line_qual[2] = {"",""};
				
				while (!temp_infile[0].eof()) {  
					
					bool end_reached = false;
					
					for (int k = 0; k<num_pair; k++) {
						
						
						getline(temp_infile[k], line);
						
						trim2(line);
						
						if (line.length() == 0) {
							
							end_reached = true;
							
							continue;
						}
						
						if (line_num%4 == 1) {
							line_good[k] = false; 
							
							if (line[0] != '@') {
								cout << "Malformed FASTQ reads file: " << *(reads_list_it[k]) << endl;
								cout << "Line: " << line_num << endl;
								cout << line << endl;
								cout << "Should begin with an @" << endl;
								debug_out << "Malformed FASTQ reads file: " << *(reads_list_it[k]) << endl;
								debug_out << "Line: " << line_num << endl;
								debug_out << line << endl;
								debug_out << "Should begin with an @" << endl;
								exit(1);
							}else {
								size_t name_end_loc = line.find_first_of(" \t/", 1);
								line_name[k] = line.substr(1, name_end_loc-1); 
								
							}
							
						}else if(line_num%4 == 2){
							
							if (line.length() > head_clip_length) {
								if (fullread_length > 0){
									line = line.substr(head_clip_length, fullread_length-head_clip_length);   
								}else if(head_clip_length > 0 ){
									line = line.substr(head_clip_length);   
								}else {
									
								}
								
								if (!(line[0] == 'A' || line[0] == 'C' || line[0] == 'T' || line[0] == 'G'
									  || line[0] == 'a' || line[0] == 'c' || line[0] == 't' || line[0] == 'g'
									  || line[0] == 'N' || line[0] == 'n')) {
									cout << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
									cout << "Ignoring line" << endl;
									cout << line << endl;
									debug_out << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
									debug_out << "Ignoring line" << endl;
									debug_out << line << endl;
									
									line_good[k] = false;
								}else {
									
									if (line.length() >= 50) {
										
										line_seq[k] = line;
										
										line_good[k] = true;
									}else {
										cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
										debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
										
									}
								}
								
							}else {
								cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
								debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
							}
							
						}else if (line_num%4 == 3) {
							if (line[0] != '+') {
								cout << "Malformed FASTQ reads file: " << *(reads_list_it[k]) << endl;
								cout << "Line: " << line_num << endl;
								cout << line << endl;
								cout << "Should begin with a +" << endl;
								debug_out << "Malformed FASTQ reads file: " << *(reads_list_it[k]) << endl;
								debug_out << "Line: " << line_num << endl;
								debug_out << line << endl;
								debug_out << "Should begin with a +" << endl;
								exit(1);
							}
						}else {
							
							
							if (quality_format == 0) { 
								
							}else if(quality_format == 1){
								phred642phred33(line);
							}else { 
								solexa2phred33(line);
							}
							
							if (line.length() > head_clip_length) {
								if (fullread_length > 0){
									line = line.substr(head_clip_length, fullread_length-head_clip_length);   
								}else if(head_clip_length > 0){
									line = line.substr(head_clip_length);   
								}else {
									
								}
							}
							
							line_qual[k] = line;  
							
							
						}
						
						
						
						
					}
					
					
					
					if (!end_reached && line_good[0] && (line_good[1]||num_pair==1) && line_num%4 == 0) { 
						
						for (int k = 0; k<num_pair; k++) {
							temp_quality_outfile[k] << line_qual[k] << "\n"; 
							
							temp_outfile[k] << line_seq[k] << "\n";  
							temp_name_outfile[k] << line_name[k] <<"[" << read_idx << "]\n";
						}
					}
					
					line_num++;
					read_idx++;
				}
				
				for (int k = 0; k<num_pair; k++) {
					
					temp_infile[k].close();
					temp_infile[k].clear();
					temp_outfile[k].close();
					temp_outfile[k].clear();
					temp_name_outfile[k].close();
					temp_name_outfile[k].clear();
					temp_quality_outfile[k].close();
					temp_quality_outfile[k].clear();
					
					(reads_list_it[k])++;
				}
				
				read_count++;
				
			}
			
		}else if (read_format.compare("FASTA") == 0) {
			
			
			
			
			
			
			uint_fast32_t read_count = 1;
			vector<string>::iterator reads_list_it[2];
			reads_list_it[0]= reads_filename[0].begin();
			reads_list_it[1]= reads_filename[1].begin();
			
			uint_fast32_t read_idx = 1;
			
			while (reads_list_it[0] != reads_filename[0].end()) { 
				ofstream temp_outfile[2];
				ofstream temp_name_outfile[2];
				ofstream temp_quality_outfile[2];
				
				ifstream temp_infile[2];
				
				for (int k = 0; k<num_pair; k++) {
					
					
					temp_infile[k].open(reads_list_it[k]->c_str(), ios::in);
					
					string new_filename = temp_path+"read_"+IntToStr(read_count)+"_"+IntToStr(k+1);
					temp_outfile[k].open(new_filename.c_str(), ios::out);
					
					*(reads_list_it[k]) = new_filename;
					
					
					temp = new_filename + ".names";
					temp_name_outfile[k].open(temp.c_str(), ios::out);
					temp = new_filename + ".quals";
					temp_quality_outfile[k].open(temp.c_str(), ios::out);
					
				}
				
				uint_fast32_t line_num = 1;
				
				bool line_good[2] = {false,false};
				string line_name[2] = {"",""};
				string line_seq[2] = {"",""};
				
				
				while (!temp_infile[0].eof()) {  
					
					bool end_reached = false;
					
					for (int k = 0; k<num_pair; k++) {
						
						
						getline(temp_infile[k], line);
						
						trim2(line);
						
						if (line.length() == 0) {
							end_reached = true;
							continue;
						}
						
						if (line_num%2 == 1) {
							line_good[k] = false; 
							
							if (line[0] != '>') {
								cout << "Malformed FASTA reads file: " << *(reads_list_it[k]) << endl;
								cout << "Line: " << line_num << endl;
								cout << line << endl;
								cout << "Should begin with an >" << endl;
								debug_out << "Malformed FASTA reads file: " << *(reads_list_it[k]) << endl;
								debug_out << "Line: " << line_num << endl;
								debug_out << line << endl;
								debug_out << "Should begin with an >" << endl;
								exit(1);
							}else {
								size_t name_end_loc = line.find_first_of(" \t/", 1);
								line_name[k] = line.substr(1, name_end_loc-1); 
								
							}
							
						}else if(line_num%2 == 0){
							
							if (line.length() > head_clip_length) {
								if (fullread_length > 0){
									line = line.substr(head_clip_length, fullread_length-head_clip_length);   
								}else if(head_clip_length > 0){
									line = line.substr(head_clip_length);   
								}else {
									
								}
								
								if (!(line[0] == 'A' || line[0] == 'C' || line[0] == 'T' || line[0] == 'G'
									  || line[0] == 'a' || line[0] == 'c' || line[0] == 't' || line[0] == 'g'
									  || line[0] == 'N' || line[0] == 'n')) {
									cout << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
									cout << "Ignoring line" << endl;
									cout << line << endl;
									debug_out << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
									debug_out << "Ignoring line" << endl;
									debug_out << line << endl;
									
									line_good[k] = false;
								}else {
									
									if (line.length() >= 50) {
										
										line_seq[k] = line;
										
										line_good[k] = true;
									}else {
										cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
										debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
									}
								}
								
							}else {
								cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
								debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
							}
							
						}

					}
					
					
					
					if (!end_reached && line_good[0] && (line_good[1]||num_pair==1) && line_num%2 == 0) { 
						
						for (int k = 0; k<num_pair; k++) {
							temp_quality_outfile[k] << (char)0 << "\n"; 
							
							temp_outfile[k] << line_seq[k] << "\n";  
							temp_name_outfile[k] << line_name[k] << "[" << read_idx << "]\n";
						}
					}
					
					line_num++;
					read_idx++;
				}
				
				for (int k = 0; k<num_pair; k++) {
					
					temp_infile[k].close();
					temp_infile[k].clear();
					temp_outfile[k].close();
					temp_outfile[k].clear();
					temp_name_outfile[k].close();
					temp_name_outfile[k].clear();
					temp_quality_outfile[k].close();
					temp_quality_outfile[k].clear();
					
					(reads_list_it[k])++;
				}
				
				read_count++;
				
			}
			
		}else if (read_format.compare("RAW") == 0) {
			
			
			
			
			
			
			uint_fast32_t read_count = 1;
			vector<string>::iterator reads_list_it[2];
			reads_list_it[0]= reads_filename[0].begin();
			reads_list_it[1]= reads_filename[1].begin();
			
			uint_fast32_t read_idx = 1;
			
			while (reads_list_it[0] != reads_filename[0].end()) { 
				ofstream temp_outfile[2];
				ofstream temp_name_outfile[2];
				ofstream temp_quality_outfile[2];
				
				ifstream temp_infile[2];
				
				for (int k = 0; k<num_pair; k++) {
					
					
					temp_infile[k].open(reads_list_it[k]->c_str(), ios::in);
					
					string new_filename = temp_path+"read_"+IntToStr(read_count)+"_"+IntToStr(k+1);
					temp_outfile[k].open(new_filename.c_str(), ios::out);
					
					*(reads_list_it[k]) = new_filename;
					
					
					temp = new_filename + ".names";
					temp_name_outfile[k].open(temp.c_str(), ios::out);
					temp = new_filename + ".quals";
					temp_quality_outfile[k].open(temp.c_str(), ios::out);
					
				}
				
				uint_fast32_t line_num = 1;
				
				bool line_good[2] = {false,false};
				
				string line_seq[2] = {"",""};
				
				
				while (!temp_infile[0].eof()) {  
					
					bool end_reached = false;
					
					for (int k = 0; k<num_pair; k++) {
						
						
						getline(temp_infile[k], line);
						
						trim2(line);
						
						if (line.length() == 0) {
							end_reached = true;
							continue;
						}
						
						
						
						if (line.length() > head_clip_length) {
							if (fullread_length > 0){
								line = line.substr(head_clip_length, fullread_length-head_clip_length);   
							}else if(head_clip_length > 0){
								line = line.substr(head_clip_length);   
							}else {
								
							}
							
							if (!(line[0] == 'A' || line[0] == 'C' || line[0] == 'T' || line[0] == 'G'
								  || line[0] == 'a' || line[0] == 'c' || line[0] == 't' || line[0] == 'g'
								  || line[0] == 'N' || line[0] == 'n')) {
								cout << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
								cout << "Ignoring line" << endl;
								cout << line << endl;
								debug_out << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
								debug_out << "Ignoring line" << endl;
								debug_out << line << endl;
								
								line_good[k] = false;
							}else {
								
								if (line.length() >= 50) {
									
									line_seq[k] = line;
									
									line_good[k] = true;
								}else {
									cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
									debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
								}
							}
							
						}else {
							cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
							debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
						}
						
					}
					
					
					
					if (!end_reached && line_good[0] && (line_good[1]||num_pair==1)) { 
						
						for (int k = 0; k<num_pair; k++) {
							temp_quality_outfile[k] << (char)0 << "\n"; 
							
							temp_outfile[k] << line_seq[k] << "\n";  
							temp_name_outfile[k] << "R["<< read_idx << "]\n";
						}
					}
					
					line_num++;
					read_idx++;
				}
				
				for (int k = 0; k<num_pair; k++) {
					
					temp_infile[k].close();
					temp_infile[k].clear();
					temp_outfile[k].close();
					temp_outfile[k].clear();
					temp_name_outfile[k].close();
					temp_name_outfile[k].clear();
					temp_quality_outfile[k].close();
					temp_quality_outfile[k].clear();
					
					(reads_list_it[k])++;
				}
				
				read_count++;
			}
			
		}
		
		
		
		
		
		
		
	}
	
	
	
	
	if(mapping_opt == '#'){
		int read_count = 1;
		vector<string>::iterator reads_list_it[2];
		reads_list_it[0]= reads_filename[0].begin();
		reads_list_it[1]= reads_filename[1].begin();
		
		while (reads_list_it[0] != reads_filename[0].end()) { 
			
			for (int k = 0; k<num_pair; k++) {
				
				
				string new_filename = temp_path+"read_"+IntToStr(read_count)+"_"+IntToStr(k+1);
				
				*(reads_list_it[k]) = new_filename;
				
				
			}
			
			for (int k = 0; k<num_pair; k++) {
				
				(reads_list_it[k])++;
			}
			
			read_count++;
			
		}
	}
	
	
	
	
	
	ofstream reads_list_temp_file;
	temp = temp_path+"reads_list1";
	reads_list_temp_file.open(temp.c_str(), ios::out);
	readslist_it = reads_filename[0].begin();
	while (readslist_it != reads_filename[0].end()) {
		reads_list_temp_file << *readslist_it << endl;
		readslist_it++;
	}
	reads_list_temp_file.close();
	reads_list_temp_file.clear();

	if(num_pair == 2){
		temp = temp_path+"reads_list2";
		reads_list_temp_file.open(temp.c_str(), ios::out);
		readslist_it = reads_filename[1].begin();
		while (readslist_it != reads_filename[1].end()) {
			reads_list_temp_file << *readslist_it << endl;
			readslist_it++;
		}
		reads_list_temp_file.close();
		reads_list_temp_file.clear();

	}
	
	
	
	
	
	if(mapping_opt != '#'){
		
		
		
		cout << "Extracting 25-mers... "  << endl;
		debug_out << "Extracting 25-mers... "  << endl;
		
		
		gettimeofday(&temp_start_tv, NULL);
		
		
		
		
		
		
		ofstream map_out;
		
		temp = temp_path+map_filename;
		map_out.open(temp.c_str(), ios::out);
		
		
		readslist_it = reads_filename[0].begin();
		while (readslist_it != reads_filename[0].end()) {

			output_index(*readslist_it,map_out);

			readslist_it++;
		}
		
		if(num_pair == 2){
			
			readslist_it = reads_filename[1].begin();
			while (readslist_it != reads_filename[1].end()) {
				
				output_index(*readslist_it,map_out);			

				readslist_it++;
			}
			
		}
		
		
		
		
		
		
		
		map_out << flush;
		
		map_out.close();
		map_out.clear();
		
		
		
		gettimeofday(&tv, NULL);
		cout<<"Total 25-mer extraction section time: "<<diffclock(temp_start_tv,tv)<<" s."<<endl;
		debug_out<<"Total 25-mer extraction section time: "<<diffclock(temp_start_tv,tv)<<" s."<<endl;
		gettimeofday(&temp_start_tv, NULL);
		
	}
	
	
	
	
	
	
	
	
	if (mapping_opt == 'S'){
		cout << "_________" << endl;
		cout << "Getting ready to do the SeqMap mapping!..." << endl;
		cmd= "cat "+ genome_direc + chromosome_wildcard + " > "+ temp_path +"genome.fa";
		
		debug_out << "_________" << endl;
		debug_out << "Getting ready to do the SeqMap mapping!..." << endl;
		debug_out << "COMMAND: " << cmd << endl;
		
		system_status = system(cmd.c_str());
		
	}else if (mapping_opt == 'E'){
		cout << "_________" << endl;
		cout << "Getting ready to do the Eland mapping!..." << endl;
		cmd =  package_path_filename.first+"squashGenome " + temp_path + " " +genome_direc+chromosome_wildcard;
		
		debug_out << "_________" << endl;
		debug_out << "Getting ready to do the Eland mapping!..." << endl;
		debug_out << "COMMAND: " << cmd << endl;
		
		system_status = system(cmd.c_str());
	}
	
	
	gettimeofday(&tv, NULL);
	cout<<"Pre-mapping section time: "<< diffclock(start_tv,tv) <<" s."<<endl;
	debug_out<<"Pre-mapping section time: "<< diffclock(start_tv,tv) <<" s."<<endl;
	gettimeofday(&temp_start_tv, NULL);
	
	
	if(mapping_opt != '#'){
		cout << "_________" << endl;
		cout << "Doing the initial mapping!... " << endl;
		debug_out << "_________" << endl;
		debug_out << "Doing the initial mapping!... " << endl;
	}
	
	gettimeofday(&map_temp_start_tv, NULL);
	
	string output_file_path = temp_path+map_filename + ".out";
	
	if(mapping_opt == 'S'){
		
		
		cmd = package_path_filename.first+"seqmap " + "2" + " " 
		+ temp_path+map_filename + " " + temp_path + "genome.fa"+ " "+ output_file_path +" /eland:3";
		
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		gettimeofday(&tv, NULL);
		cout<<"Total Seqmap mapping section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
		debug_out<<"Total Seqmap mapping section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
		gettimeofday(&map_temp_start_tv, NULL);

	}else if(mapping_opt == 'E'){
		
		
		cmd = package_path_filename.first+"eland_25 " + temp_path + map_filename +" " 
		+temp_path+ " " + output_file_path + " --multi="+IntToStr(max_multi_hit)+" ";
		
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		gettimeofday(&tv, NULL);
		cout<<"Total Eland mapping section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
		debug_out<<"Total Eland mapping section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
		gettimeofday(&map_temp_start_tv, NULL);
		
		
		
	}else if(mapping_opt == 'B'){
		
		
		
		temp = " -S -k "+IntToStr(max_multi_hit)+" -m "+IntToStr(max_multi_hit)+" -v 2 -r -p "+IntToStr(num_threads)+" --best --strata";
		
		cmd = package_path_filename.first+"bowtie " + bowtie_opts + temp + " "+ bowtie_base_dir + " "
		+ temp_path+map_filename +" "  + temp_path + map_filename + "_unsorted 2> " + out_path+debug_path+"bowtie_out.log";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		if (system_status != 0) {
			cout << "ERROR: Bowtie execution failed with error: " << system_status << endl;
			debug_out << "ERROR: Bowtie execution failed with error: " << system_status << endl;

			exit(system_status);
		}
		
		
		
		cmd = package_path_filename.first+"sortsam -idx "+temp_path + map_filename + "_unsorted " +output_file_path;
		
		debug_out << "COMMAND:" << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = "rm "  +temp_path + map_filename + "_unsorted " ;
		
		debug_out << "COMMAND:" << cmd << endl;
		system_status = system(cmd.c_str());
		
		
		gettimeofday(&tv, NULL);
		cout<<"Total Bowtie mapping section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
		debug_out<<"Total Bowtie mapping section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
		gettimeofday(&map_temp_start_tv, NULL);

	}
	
	cout << "Generating .t file!..." << endl;
	debug_out << "Generating .t file!..." << endl;
	
	
	if(mapping_opt == 'S' || mapping_opt == 'E'){
		
		string infile = output_file_path;
		string outfile = output_file_path+".t";
		processmultimap(infile,outfile,max_multi_hit);
		
		gettimeofday(&tv, NULL);
		cout<<"Total .t creation section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
		debug_out<<"Total .t creation section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
	}
	
	
	if(mapping_opt == 'B'){
		
		string infile = output_file_path;
		string outfile = output_file_path+".t";
		processbowtie(infile,outfile);
		
		gettimeofday(&tv, NULL);
		cout<<"Total .t creation section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
		debug_out<<"Total .t creation section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
	}
	
	if(mapping_opt != '#'){
		cout << "Removing mapping results!..." << endl;
		debug_out << "Removing mapping results!..." << endl;
		
		cmd = "rm "  + output_file_path ;
		
		debug_out << "COMMAND:" << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = "rm "  + temp_path + map_filename ;
		
		debug_out << "COMMAND:" << cmd << endl;
		system_status = system(cmd.c_str());
	}
	
	cout << "Creating mapping index!..." << endl;
	debug_out << "Creating mapping index!..." << endl;
	
	gettimeofday(&map_temp_start_tv, NULL);
	

	
	
	
	cout << "Generating mapping file!..." << endl;
	debug_out << "Generating mapping file!..." << endl;
	
	ifstream dot_t_file;
	temp = output_file_path+".t";
	dot_t_file.open(temp.c_str(), ios::in);
	
	if (!dot_t_file.is_open()) {
		cout << "ERROR: Fatal error opening .t file: " << output_file_path+".t" << endl;
		cout << "Have the reads been mapped?" << endl;
		debug_out << "ERROR: Fatal error opening .t file: " << output_file_path+".t" << endl;
		debug_out << "Have the reads been mapped?" << endl;
		exit(2);
	}
	
	string dot_t_line;
	uint_fast32_t mapping_index = 1;
	getline(dot_t_file, dot_t_line);
	trim2(dot_t_line);
	if (dot_t_line.length() == 0) {
		cout << "ERROR: There is something wrong with the mapping results... unexpected length" << endl;
		exit(2);
	}
	
	for(int k = 0; k<num_pair;k++){
		readslist_it = reads_filename[k].begin();
		while (readslist_it != reads_filename[k].end()) {
			ifstream index_file;
			temp = *readslist_it + ".index";
			index_file.open(temp.c_str(), ios::in);
			
			if (!index_file.is_open()) {
				cout << "ERROR: Fatal error opening index file ... "  << temp << endl;
				debug_out << "ERROR: Fatal error opening index file ... "  << temp << endl;
				exit(2);
			}
			
			ofstream out_index_file;
			temp = *readslist_it + ".index.t";
			out_index_file.open(temp.c_str(), ios::out);
			
			while (!index_file.eof()) {
				getline(index_file, line);
				
				trim2(line);
				
				if (line.length() == 0) {
					continue;
				}
				
				
				int curr_read_length = atoi(line.c_str());
				
				
				out_index_file << curr_read_length << "\n";
				
				vector<coord_t> suffix_list = designsuffix(curr_read_length);
				
				for (unsigned int i = 0; i<suffix_list.size(); i++) {
					
					
					vector<dot_t_t> read_segment[2];
					
					for(int j = 0; j<2;j++){
						bool reached_end = false; 
						read_segment[j] = vector<dot_t_t>();
						while (!reached_end && !dot_t_file.eof()) {
						
							
							vector<string> line_list = split(dot_t_line, '\t');
							unsigned int read_idx = (unsigned int)atoi(line_list[0].c_str());
							
							if (read_idx != mapping_index) {
								reached_end = true;
								mapping_index++;
								
							}else {

								if (line_list.size() == 2) {
									
									
								}else if (line_list.size() > 2){
									
									
									
									int curr_num_mismatch = -1;
									switch (line_list[1][1]) {  
										case '0':
											curr_num_mismatch = 1;
											break;
										case '1':
											curr_num_mismatch = 2;
											break;
										case '2':
											curr_num_mismatch = 3;
											break;
										default:
											cout << "ERROR: num mismatch not matching... " << line << endl;
											debug_out << "ERROR: num mismatch not matching... " << line << endl;
											exit(2);
											break;
									}
									
									int chr_pos = atoi(line_list[3].c_str());
									
									int direction = 0;
									switch (line_list[4][0]) {
										case 'F':
											direction = 1;
											break;
										case 'R':
											direction = -1;
											break;
										default:
											cout << "ERROR: direction not matching... " << line << endl;
											debug_out << "ERROR: direction not matching... " << line << endl;
											exit(2);
											break;
									}
									
									dot_t_t contents;
									contents.mismatch_dir = curr_num_mismatch*direction;
									contents.location = chr_pos;
									contents.chr_name = line_list[2];
									
									
									
									
									read_segment[j].push_back(contents);
									
								}else {
									
									cout << "ERROR: Fatal error parsing .t file: " << line << endl;
									debug_out << "ERROR: Fatal error parsing .t file: " << line << endl;
									exit(2);
								}
								
								
								getline(dot_t_file, dot_t_line);
								
								trim2(dot_t_line);
							}
						}
					}
					
					
					
					
					
					
					
					
					
					
					
					for (int j = 0; j < 2; j++) {
						
						
						
						out_index_file << read_segment[j].size() << "\n";
						
						vector<dot_t_t>::iterator dot_t_it = read_segment[j].begin();
						while (dot_t_it != read_segment[j].end()) {
							out_index_file << dot_t_it->chr_name << "\t" << dot_t_it->location << "\t" << (int)dot_t_it->mismatch_dir << "\n";
							
							dot_t_it++;
						}
					}
				}
				
			}
			
			
			index_file.close();
			index_file.clear();
			
			out_index_file.close();
			out_index_file.clear();
			
			readslist_it++;
		}
	}
	
	dot_t_file.close();
	dot_t_file.clear();
	
	
	
	
	
	gettimeofday(&tv, NULL);
	cout<<"Total mapping index creation section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
	debug_out<<"Total mapping index creation section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;

	
	gettimeofday(&tv, NULL);
	cout<<"Total Mapping section time: "<<diffclock(temp_start_tv,tv)<<" s."<<endl;
	debug_out<<"Total Mapping section time: "<<diffclock(temp_start_tv,tv)<<" s."<<endl;

	gettimeofday(&temp_start_tv, NULL);
	

	
	
	
	

	
	
	ofstream list_file;
	temp = temp_path + good_dict_list_filename;
	list_file.open(temp.c_str(), ios::out);
	list_file.close();
	
	
	
	ref_map_it = ref_map.begin();
	
	int running_index = 0;
	int num_running = 0;
	bool running_success = true;
	int num_fork_fail = 0;
	
	pid_t *pid_list = new pid_t[ref_map.size()];
	for(unsigned int i = 0;i<ref_map.size();i++){
		pid_list[i] = -1; 
	}
	
	while (ref_map_it != ref_map.end() && running_success && num_fork_fail < 10) {

		cout << "Running SpliceMap on " << ref_map_it->first << " ..." << endl;
		debug_out << "Running SpliceMap on " << ref_map_it->first << " ..." << endl;
		
		cmd = package_path_filename.first+"SpliceMap "+cfg_filename+" "+ref_map_it->first;
		
		debug_out << "COMMAND: " << cmd << endl;
		
		pid_t pid = fork(); 
		
		if(pid == 0){
			
			string child_cmd = cmd;
			int child_status = system(child_cmd.c_str());
			exit(child_status); 
			
		}else if (pid < 0) { 
			cout << "ERROR: fork failed: " << pid << endl;
			debug_out << "ERROR: fork failed: " << pid << endl;
			
			num_fork_fail++;
			if (num_fork_fail >= 10) {
				cout << "ERROR: Too many failed fork attempts ... exiting" << endl;
				exit(3);
			}
			
			cout << "Trying again... in 2 seconds (cross fingers)" << endl;
			debug_out << "Trying again... in 2 seconds (cross fingers)" <<  endl;
			sleep(2);
			
			continue;
			
			
		}else{
			
			
			pid_list[running_index] = pid;
			
			sleep(1); 
			
		}


		running_index++;
		
		
		do {
			num_running = 0;
			
			map<string,reference_t>::iterator ref_map_it2 = ref_map.begin();
			for(unsigned int i = 0;i<ref_map.size();i++){
				if(pid_list[i] > 0){ 
					int status = 0;
					pid_t curr_pid = waitpid(pid_list[i], &status, WNOHANG);
					
					
					if (curr_pid != 0) {
						if (status == 0) {
							pid_list[i] = 0; 
						}else {
							cout << "ERROR: SpliceMapping of " << ref_map_it2->first << " failed." << endl;
							cout << "Waiting for other processes to finish then exiting..." << endl;
							debug_out << "ERROR: SpliceMapping of " << ref_map_it2->first << " failed." << endl;
							debug_out << "Waiting for other processes to finish then exiting..." << endl;
							running_success = false;
						}
					}
				}

				
				
				if(pid_list[i] > 0){
					num_running++;
				}
				
				ref_map_it2++;
			} 

			
			
			sleep(5);
		}while (num_running >= num_chromosome_together);
		

		ref_map_it++;
	}
	
	num_running = 0;
	
	while (running_index > num_running) {
		num_running = 0;
		
		map<string,reference_t>::iterator ref_map_it2 = ref_map.begin();
		for(unsigned int i = 0;i<ref_map.size();i++){
			if(pid_list[i] > 0){ 
				int status = 0;
				pid_t curr_pid = waitpid(pid_list[i], &status, WNOHANG);
				if (curr_pid != 0) {
					if (status == 0) {
						pid_list[i] = 0; 
					}else {
						cout << "ERROR: SpliceMapping of " << ref_map_it2->first << " failed. status = " << status << endl;
						cout << "Waiting for other processes to finish then exiting..." << endl;
						debug_out << "ERROR: SpliceMapping of " << ref_map_it2->first << " failed. status = " << status << endl;
						debug_out << "Waiting for other processes to finish then exiting..." << endl;
						running_success = false;
					}
				}
			}
			
			
			
			if(pid_list[i] == 0){ 
				num_running++;
			}
			
			
			ref_map_it2++;
		}
		
		
		
		sleep(5);
	}
	
	delete [] pid_list;
	
	if (!running_success) {
		cout << "ERROR: SpliceMap failed. Exiting..." << endl;
		debug_out << "ERROR: SpliceMap failed. Exiting..." << endl;
		exit(2);
	}
	
	gettimeofday(&tv, NULL);
	cout<<"Total SpliceMap section time: "<<diffclock(temp_start_tv,tv)<<" s."<<endl;
	debug_out<<"Total SpliceMap section time: "<<diffclock(temp_start_tv,tv)<<" s."<<endl;
	gettimeofday(&temp_start_tv, NULL);
	

	

	
	cout << "Merging results and computing junctions..." << endl;
	
	
	cmd = package_path_filename.first+"amalgamateSAM "+temp_path+" "+ temp_path+"junction ";
	
	debug_out << "COMMAND: " << cmd << endl;
	system_status = system(cmd.c_str());
	
	if (system_status != 0) {
		cout << "ERROR: Amalgamation failed = " << system_status << endl;
		debug_out << "ERROR: Amalgamation failed = " << system_status << endl;
		exit(2);
	}
	
	
	
	cout << "Removing some temporary files..." << endl;
	
	ref_map_it = ref_map.begin();
	while (ref_map_it != ref_map.end()) {
		
		temp = temp_path+ref_map_it->second.file_name + "_" + LongintToStr(ref_map_it->second.file_index_start) + ".sam";
		cmd = "rm "  +temp ;
		
		debug_out << "COMMAND:" << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = "rm "  +temp +"_uniq" ;
		
		debug_out << "COMMAND:" << cmd << endl;
		system_status = system(cmd.c_str());
		
		ref_map_it++;
	}
	
	
	
	cout << "Sorting SAM file..." << endl;
	
	
	
	cmd = package_path_filename.first+"sortsam -pos "+temp_path+"junction.sam " +temp_path+"junction.sam_sorted";
	
	debug_out << "COMMAND:" << cmd << endl;
	system_status = system(cmd.c_str());
	
	cmd = "mv "  +temp_path+"junction.sam_sorted " +temp_path+"junction.sam" ;
	
	debug_out << "COMMAND:" << cmd << endl;
	system_status = system(cmd.c_str());
	
	
	
	cout << "Writing coverage..." << endl;
	
	
	cmd = package_path_filename.first+"precipitateSAM "+temp_path+"junction.sam "+ temp_path+"junction ";
	
	debug_out << "COMMAND: " << cmd << endl;
	system_status = system(cmd.c_str());
	
	if (system_status != 0) {
		cout << "ERROR: Precipitation failed = " << system_status << endl;
		debug_out << "ERROR: Precipitation failed = " << system_status << endl;
		exit(2);
	}
	
	
	cout << "Running filters and formatting output..." << endl;
	
	cmd = package_path_filename.first+"uniqueJunctionFilter "+temp_path+"junction.bed "+temp_path+"junction_nUM.bed ";
	
	debug_out << "COMMAND: " << cmd << endl;
	system_status = system(cmd.c_str());
	
	if (ref_filename.length()!=0){
		
		cmd=package_path_filename.first+"findNovelJunctions " + ref_filename + " "+temp_path+"junction.bed";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd=package_path_filename.first+"findNovelJunctions " + ref_filename + " "+temp_path+"junction_nUM.bed";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = package_path_filename.first+"colorJunction "+temp_path+"junction.bed "+temp_path+"junction.bed.new.bed";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = package_path_filename.first+"colorJunction "+temp_path+"junction.bed.new.bed "+temp_path+"junction.bed.new.bed";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = package_path_filename.first+"colorJunction "+temp_path+"junction_nUM.bed "+temp_path+"junction_nUM.bed.new.bed";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = package_path_filename.first+"colorJunction "+temp_path+"junction_nUM.bed.new.bed "+temp_path+"junction_nUM.bed.new.bed";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = package_path_filename.first+"statSpliceMap "+temp_path+"junction.bed "+temp_path+"junction.bed.new.bed > "+temp_path+"log2";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		
		cmd = "cat "+temp_path+"log2 > "+out_path+"log";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		
		cmd = "mv "+temp_path+"junction_color.bed "+temp_path+"junction.bed.new_color.bed "+temp_path+"junction_nUM_color.bed "+temp_path+"junction_nUM.bed.new_color.bed "+out_path;
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
	
		cmd = "mv "+out_path+"junction.bed.new_color.bed "+out_path+"junction_color.new.bed ";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = "mv "+out_path+"junction_nUM.bed.new_color.bed "+out_path+"junction_nUM_color.new.bed ";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		
	}else {
		
		cmd = package_path_filename.first+"colorJunction "+temp_path+"junction.bed ";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = package_path_filename.first+"colorJunction "+temp_path+"junction_nUM.bed ";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = package_path_filename.first+"statSpliceMap "+temp_path+"junction.bed > " + temp_path+"log2";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		
		cmd = "cat "+temp_path+"log2 > "+out_path+"log";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
		
		cmd = "mv "+temp_path+"junction_color.bed "+temp_path+"junction_nUM_color.bed "+out_path;
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());

	}
	
	
	
	cmd = "mv "+temp_path+"junction_up.wig "+temp_path+"junction_down.wig "+temp_path+"junction_all.wig "+out_path;
	
	debug_out << "COMMAND: " << cmd << endl;
	system_status = system(cmd.c_str());
	
	

	if (print_sam || cufflinks) {
		
		cmd = "mv "+temp_path+"junction.sam "+out_path+"good_hits.sam";
		
		debug_out << "COMMAND: " << cmd << endl;
		system_status = system(cmd.c_str());
	}
	
	
	cmd = "mv "+out_path+"junction_up.wig "+out_path+"coverage_up.wig ";
	
	debug_out << "COMMAND: " << cmd << endl;
	system_status = system(cmd.c_str());
	
	cmd = "mv "+out_path+"junction_down.wig "+out_path+"coverage_down.wig ";
	
	debug_out << "COMMAND: " << cmd << endl;
	system_status = system(cmd.c_str());
	
	
	cmd = "mv "+out_path+"junction_all.wig "+out_path+"coverage_all.wig ";
	
	debug_out << "COMMAND: " << cmd << endl;
	system_status = system(cmd.c_str());
	
	
	gettimeofday(&tv, NULL);
	cout<<"Total merge/novel/stat/etc. section time: "<<diffclock(temp_start_tv,tv)<<" s."<<endl;
	debug_out<<"Total merge/novel/stat/etc. section time: "<<diffclock(temp_start_tv,tv)<<" s."<<endl;

	
	gettimeofday(&tv, NULL);
	cout << "______________________" << endl;
	cout<<"Total All Chromosone Execution Time: "<<diffclock(start_tv,tv)<<" s."<<endl;
	debug_out << "______________________" << endl;
	debug_out<<"Total All Chromosone Execution Time: "<<diffclock(start_tv,tv)<<" s."<<endl;

	
	
	debug_out.close();
	debug_out.clear();
	
	
    return 0;
}





inline int sam_char2int(char c)
{
	if (c == '0') {
		return 0;
	}else if (c == '1') {
		return 1;
	}else if(c == '2'){
		return 2;
	}
	
	cout << "Warning: Unexpected character in mismatch string, " << c << endl;
	
	return -1;

}

inline bool sam_is_mapped(unsigned int flag)
{
	return !(flag & 4);
}

inline string sam_get_direction(unsigned int flag)
{
	
	if (flag & 16){
		return "R";
	}else {
		return "F";
	}

	
}

inline void processbowtie(string &input_filename, string &output_filename)
{
	ifstream input;
	ofstream output;
	

	string temp;
	string line;
	list<sam_mismatch_t> sam_buf;
	sam_mismatch_t curr_sam;
	int curr_mismatch;
	int curr_read_idx;
	unsigned int curr_flag;
	int tabloc;
	int tabloc2;
	int min_mismatch;
	
	input.open(input_filename.c_str(), ios::in);
	if(!input.is_open()){
		cout << "ERROR: I'm sorry, cannot open the mapping output file, " << input_filename << endl;
		exit(1);
	}
	
	
	while (!input.eof()) {
		getline(input, line);
		trim2(line);
		if (line.length() == 0) {
			continue;
		}
		
		if (line[0] == '@') {
			
		}else {
			break;
		}
		
	}
	
	
	output.open(output_filename.c_str(), ios::out);

	
	tabloc = (int) line.find('\t');
	tabloc2 = (int) line.find('\t',tabloc+1);
	curr_flag = atoi(line.substr(tabloc+1, tabloc2-tabloc-1).c_str());
	
	
	curr_read_idx = atoi(line.substr(0,tabloc).c_str()) + 1;  
	
	if (sam_is_mapped(curr_flag)) {
		
		curr_mismatch = sam_char2int(line[line.length()-1]); 
		
		curr_sam.num_mismatch = curr_mismatch;
		curr_sam.data = line;
		sam_buf.push_back(curr_sam);
		
		
	}else {
		curr_mismatch = -1;
	}
	
	
	while (!input.eof()) {
		getline(input, line);
		trim2(line);
		
		if (line.length() == 0){
			continue;  
		}
		
		tabloc = (int) line.find('\t');
		tabloc2 = (int) line.find('\t',tabloc+1);
		curr_flag = atoi(line.substr(tabloc+1, tabloc2-tabloc-1).c_str());
		int temp_curr_read_idx = atoi(line.substr(0,tabloc).c_str()) + 1;
		
		
		
		if (temp_curr_read_idx != curr_read_idx) {
			
			
			
			
			if (curr_mismatch == -1) {
				
				output << curr_read_idx << "\t" << "NM" << endl;
				
			}else{ 
				list<sam_mismatch_t>::iterator buf_it;

				buf_it = sam_buf.begin();
				
				bool curr_unique = (sam_buf.size()==1);
				
				while (buf_it != sam_buf.end()) {
					temp = buf_it->data;
					string curr_chr_name;
					int curr_pos;
					
					tabloc = (int) temp.find('\t'); 
					tabloc2 = (int) temp.find('\t',tabloc+1); 
					
					int temp_curr_flag = atoi(temp.substr(tabloc+1, tabloc2-tabloc-1).c_str());
					
					tabloc = tabloc2;
					tabloc2 = (int) temp.find('\t',tabloc+1); 
					
					curr_chr_name = temp.substr(tabloc+1, tabloc2-tabloc-1).c_str();
					
					tabloc = tabloc2;
					tabloc2 = (int) temp.find('\t',tabloc+1); 
					
					curr_pos = atoi(temp.substr(tabloc+1, tabloc2-tabloc-1).c_str());
					
					int temp_curr_mismatch = sam_char2int(temp[temp.length()-1]); 
					
					if (curr_unique) {
						output << curr_read_idx << "\t" << "U" << temp_curr_mismatch << "\t" 
						<< curr_chr_name << "\t" << curr_pos << "\t" << sam_get_direction(temp_curr_flag) << endl;
					} else{
						output << curr_read_idx << "\t" << "R" << temp_curr_mismatch << "\t" 
						<< curr_chr_name << "\t" << curr_pos << "\t" << sam_get_direction(temp_curr_flag) << endl;
						
					}
					
					
					buf_it++;
				}
				
				
			}
			
			sam_buf.clear();
			min_mismatch = 2;
			curr_read_idx = temp_curr_read_idx;
		}
		
		tabloc = (int) line.find('\t');
		tabloc2 = (int) line.find('\t',tabloc+1);
		
		curr_read_idx = atoi(line.substr(0,tabloc).c_str()) + 1;
		
		if (sam_is_mapped(curr_flag)) {
			
			curr_mismatch = sam_char2int(line[line.length()-1]); 
			
			curr_sam.num_mismatch = curr_mismatch;
			curr_sam.data = line;
			sam_buf.push_back(curr_sam);
			
		}else {
			curr_mismatch = -1;
		}
		
		
	}
	
	
	if (curr_mismatch == -1) {
		
		output << curr_read_idx << "\t" << "NM" << endl;
		
	}else{ 
		list<sam_mismatch_t>::iterator buf_it;
		
		buf_it = sam_buf.begin();
		
		bool curr_unique = (sam_buf.size()==1);
		
		while (buf_it != sam_buf.end()) {
			temp = buf_it->data;
			string curr_chr_name;
			int curr_pos;
			
			tabloc = (int) temp.find('\t'); 
			tabloc2 = (int) temp.find('\t',tabloc+1); 
			
			int temp_curr_flag = atoi(temp.substr(tabloc+1, tabloc2-tabloc-1).c_str());
			
			tabloc = tabloc2;
			tabloc2 = (int) temp.find('\t',tabloc+1); 
			
			curr_chr_name = temp.substr(tabloc+1, tabloc2-tabloc-1).c_str();
			
			tabloc = tabloc2;
			tabloc2 = (int) temp.find('\t',tabloc+1); 
			
			curr_pos = atoi(temp.substr(tabloc+1, tabloc2-tabloc-1).c_str());
			
			int temp_curr_mismatch = sam_char2int(temp[temp.length()-1]); 
			
			if (curr_unique) {
				output << curr_read_idx << "\t" << "U" << temp_curr_mismatch << "\t" 
				<< curr_chr_name << "\t" << curr_pos << "\t" << sam_get_direction(temp_curr_flag) << endl;
			} else{
				output << curr_read_idx << "\t" << "R" << temp_curr_mismatch << "\t" 
				<< curr_chr_name << "\t" << curr_pos << "\t" << sam_get_direction(temp_curr_flag) << endl;
				
			}
			
			
			buf_it++;
		}
		
		
	}
	
	
	
	input.close();
	input.clear();
	output.close();
	output.clear();
}


inline void processmultimap(string &input_filename, string &output_filename,int max_multi_hit)
{
	ifstream input;
	ofstream output;
	stringstream buf(stringstream::in | stringstream::out); 
	
	
	
	string line;
	
	input.open(input_filename.c_str(), ios::in);
	if(!input.is_open()){
		cout << "ERROR: I'm sorry, cannot open the mapping output file, " << input_filename << endl;
		exit(1);
	}
	
	while (!input.eof()) {
		getline(input, line);
		trim2(line);
		
		if (line.length() == 0){
			continue;  
		}
		
		buf << line << "\n";
	}
	input.close();
	input.clear();

	
	output.open(output_filename.c_str(), ios::out);
	
	while (!buf.eof()) {
		
		getline(buf, line,'\t');
		
		
		
		if (line.length() == 0){
			continue;  
		}
		
		
		
		
		vector<string> line0_tokens = split(line, '-'); 
		
		string index = line0_tokens.back(); 
		
		getline(buf, line,'\t'); 
		
		
		
		
		buf >> line;
		
		
		
		
		if (line.compare("NM")!=0 && line.compare("QC")!=0) {
			int mismatch = 2;
			
			vector<string> line2_tokens = split(line, ':');
			
			int token[3];
			for (int i = 0;i<3;i++){
				token[i] = atoi(line2_tokens[i].c_str());
			}
			
			
			if (token[0] != 0) {
				mismatch = 0;
				
				if (token[0] > max_multi_hit){
					output << index << "\tR" << IntToStr(mismatch) << "\n";
				}else {
					
					
					buf >> line;
					
					vector<string> mapping_list = rearrange_eland_multi_map(line,mismatch,token[mismatch]);   
					
					for (vector<string>::iterator map_it = mapping_list.begin(); map_it < mapping_list.end(); map_it++) {
						output << index << '\t' << *map_it << "\n";
					}
				}
				
				
			}else if(token[1] != 0) {
				mismatch = 1;
				
				if (token[1] > max_multi_hit){
					output << index << "\tR" << IntToStr(mismatch) << "\n";
				}else {
					
					
					buf >> line;
					
					vector<string> mapping_list = rearrange_eland_multi_map(line,mismatch,token[mismatch]);   
					
					for (vector<string>::iterator map_it = mapping_list.begin(); map_it < mapping_list.end(); map_it++) {
						output << index << '\t' << *map_it << "\n";
					}
				}
			}else {
				mismatch = 2;
				
				if (token[2] > max_multi_hit){
					output << index << "\tR" << IntToStr(mismatch) << "\n";
				}else {
					
					
					buf >> line;
					
					vector<string> mapping_list = rearrange_eland_multi_map(line,mismatch,token[mismatch]);   
					
					for (vector<string>::iterator map_it = mapping_list.begin(); map_it < mapping_list.end(); map_it++) {
						output << index << '\t' << *map_it << "\n";
					}
				}
			}
			
			
		}else if(line.compare("NM")==0){
			output << index << "\tNM" << "\n";
		}else if(line.compare("QC")==0){
			output << index << "\tQC" << "\n";
		}
		
		
		
		char c;
		buf.get(c);
		if (c != '\n') {
			getline(buf, line,'\n'); 
			
		}
		
	}
	
	
	output << flush;
	
	output.close();
	output.clear();
	
	
	
	
	
}


inline vector<string> rearrange_eland_multi_map(string& s, int n_miss, int num_map){
	string chr_fa_name;
	string pos;
	string RU;
	string temp;
	
	vector<string> output;
	
	trim2(s);
	vector<string> s_list = split(s,',');
	
	if (num_map == 1){
		RU = "U";  
	}else {
		RU = "R";  
	}
	
	for (vector<string>::iterator it = s_list.begin(); it < s_list.end(); it++) {
		if ((*it).substr(0,3).compare("chr") == 0) {
			vector<string> it_list = split(*it, ':');
			
			chr_fa_name = it_list[0];
			temp = chr_fa_name.substr(chr_fa_name.length()-3);
			if (temp.compare(".fa") == 0) {
				chr_fa_name = chr_fa_name.substr(0, chr_fa_name.length()-3);
			}
			
			pos = it_list[1];
		}else { 
			pos = *it;
		}
		int lenpos = (int)pos.length();
		int mismatch = atoi(pos.substr(lenpos-1).c_str());
		
		if (mismatch == n_miss){
			output.push_back(RU+IntToStr(mismatch)+"\t"+chr_fa_name+"\t"+pos.substr(0, lenpos-2)+"\t"+pos.substr(lenpos-2,1));
		}
	}

	
	return output;
	
}

inline void print_usage_and_exit()
{
	cout << "usage: ./runSpliceMap run.cfg" << endl;
	cout << "  run.cfg  --  Configuration options for this run, see comments in file for details" << endl;
	cout << "See website for further details" << endl;
	exit(0);
}

 
pair<string,string> get_path_and_filename(const string& path)
{
	string::size_type lastslash = path.rfind('/');
	if (lastslash == string::npos) {
		lastslash = path.rfind('\\');
	}
	if (lastslash == string::npos) {
		return pair<string,string>("",path);
	}
	
	return pair<string,string>(path.substr(0,lastslash+1),path.substr(lastslash+1,path.length()-lastslash));
}





void output_index(string reads_filename,ofstream &map_out)
{
	string temp;
	string line;
	
	
	
	ifstream reads_in;
	
	
	reads_in.open(reads_filename.c_str(), ios::in);
	
	if (!reads_in.is_open()) {
		cout << "ERROR: Fatal error opening reads file... " << reads_filename << endl;
		exit(1);
	}
	
	ofstream index_out;
	temp = reads_filename+".index";
	index_out.open(temp.c_str(), ios::out);
	
	
	
	while (!reads_in.eof()) {
		getline(reads_in, line);
		
		trim2(line);
		
		if (line.length() == 0) {
			continue;
		}
		
		transform(line.begin(), line.end(), line.begin(), (int(*)(int))toupper); 
		
		index_out << line.length() << "\n";
		
		vector<coord_t> suffix = designsuffix((int)line.length());
		
		for (vector<coord_t>::iterator it = suffix.begin(); it < suffix.end(); it++) {
			
			
			
			
			
			map_out << line.substr(it->first - 1,25) << '\n';    
			map_out << line.substr(it->first+25 - 1,25) << '\n'; 
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
		}
		
		
	}
	
	index_out.close();
	index_out.clear();
	
	reads_in.close();
	reads_in.clear();

}

