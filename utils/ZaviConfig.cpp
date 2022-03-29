#include "ZaviConfig.hpp"
std::shared_ptr<libconfig::Config> cfg;
void initConfig(const char* config_folder,const char * config_file) {
	cfg = std::make_shared<libconfig::Config>();
	cfg->setIncludeDir(config_folder);
	cfg->setOptions(cfg->OptionAutoConvert);
	std::string name(config_folder);
	name.append(config_file);
	cfg->readFile(name);
}



std::shared_ptr<double> readTrialArray(const char * config_name,size_t size){
	libconfig::Setting & conf=cfg->lookup((const char *)cfg->lookup("trial"));
	auto pointer=(double *)calloc(size,sizeof(double));
	auto array=std::shared_ptr<double>(pointer);
	for(size_t i=0; i < size; i++){
		pointer[i]=conf[config_name][i];
	}
	return array;

}
