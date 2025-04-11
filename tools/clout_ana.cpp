#include <analyse_clfiles.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
    
    if (argc != 2) {
        printf("Usage: %s <clusterfiles_*_particles_type_*.txt>\n", argv[0]);
        return 1;
    }
    
    std::filesystem::path path(argv[1]);
    
    if (!std::filesystem::exists(path)) {
        fprintf(stderr, "File not found: %s! exiting...\n", path.string().c_str());
        return 2;
    }
    
    std::ifstream cls_files(path.string());
    
    if (!cls_files.is_open()) {
        fprintf(stderr, "Error opening file: %s! exiting...\n", path.string().c_str());
        return 3;
    }

    std::string clsdata, xmlfile, ptype, layer;
    
    path.replace_extension().filename();
    
    std::string sumary_filename = path.string() + ".summary";
    
    std::ofstream summary_file(sumary_filename);
    if (!summary_file.is_open()) {
        fprintf(stderr, "Error opening summary file: %s! exiting...\n", sumary_filename.c_str());
        return 4;
    }
    
    summary_file<<"# index NClusters minSize maxSize averageSize cls_filename\n";
       
    for ( size_t indx = 0 ; cls_files >> clsdata >> xmlfile >> ptype >> layer; ++indx)
    {
        std::filesystem::path cls_data_path(clsdata);
        std::ifstream cls_data_in(cls_data_path);
        if (!cls_data_in.is_open()) {
            fprintf(stderr, "Error opening file: %s! skipping...\n", clsdata.c_str());
            continue;
        }
        
        cls_data_path.replace_extension();
        
        int Nclusters, minSize, maxSize;
        double averageSize;
        
        std::cout<<"Analyzing: "<<cls_data_path.string()<<"\n";
        
        analyzeClusters(clsdata, cls_data_path.filename().string() + ".hist", 
                        cls_data_path.filename().string() + ".largest_cluster", 
                        Nclusters, minSize, maxSize, averageSize);
        
        summary_file<<indx<<" "
                    <<Nclusters<<" "<<minSize<<" "<<maxSize<<" "<<averageSize<<" "
                    <<cls_data_path.string()<<"\n";
        
    }
    summary_file.close();
    return 0;
}