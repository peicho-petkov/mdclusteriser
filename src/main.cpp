#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <filesystem>
#include "parser.h"
#include "clustering.h"

int main(int argc, char **argv) {
    if (argc < 6) {
        printf("Usage: %s --xml system.xml --cut <float> --types <str> ... <str> --up_down_layers --up_layer --down_layer\n", argv[0]);
        return 1;
    }

    std::vector<std::string> considered_types;
    std::vector<std::string> input_files;
    bool up_layer = false;
    bool down_layer = false;
    bool all = true;
    
    float cluster_cutoff = 1.0;

    for ( int i = 1; i < argc;) {
        std::cout<<"argv[" << i <<"] "<<argv[i]<< std::endl;
        if (!strcmp(argv[i], "--cut")) {
            for (++i; i < argc && argv[i][0] != '-'; ++i) { // Skip non-option arguments
            std::cout<<"cut argv[" << i <<"] "<<argv[i]<< std::endl;
                cluster_cutoff = atof(argv[i]);
            }
            continue;
        } else if (!strcmp(argv[i], "--types")) {
            for (++i; i < argc && argv[i][0] != '-'; ++i) { // Skip non-option arguments
            std::cout<<"types argv[" << i <<"] "<<argv[i]<< std::endl;
                considered_types.push_back(argv[i]);
            }
            continue;
        } else if (!strcmp(argv[i], "--xml")) {
            for (++i; i < argc && argv[i][0] != '-'; ++i) { // Skip non-option arguments
                std::cout<<"xml argv[" << i <<"] "<<argv[i]<< std::endl;
                input_files.push_back(argv[i]);
            }
            continue;
        } else if (!strcmp(argv[i], "--up_down_layers")) {
            for (++i; i < argc && argv[i][0] != '-'; ++i) { // Skip non-option arguments
                std::cout<<"xml argv[" << i <<"] "<<argv[i]<< std::endl;
            }
            up_layer = true;
            down_layer = true;
            all = false;        
            continue;
        } else if (!strcmp(argv[i], "--up_layer")) {
            for (++i; i < argc && argv[i][0] != '-'; ++i) { // Skip non-option arguments
                std::cout<<"xml argv[" << i <<"] "<<argv[i]<< std::endl;
            }
            up_layer = true;
            all = false;
            continue;
        } else if (!strcmp(argv[i], "--down_layer")) {
            for (++i; i < argc && argv[i][0] != '-'; ++i) { // Skip non-option arguments
                std::cout<<"xml argv[" << i <<"] "<<argv[i]<< std::endl;
            }
            down_layer = true;
            all = false;
            continue;
        } else {
            fprintf(stderr, "Invalid option: %s\n", argv[i]);
            return 1;
        }
    }

    std::map<std::string, std::ofstream> all_files_output, up_files_output, down_files_output;
    if (all)
    for (const auto &t : considered_types)
        all_files_output[t].open("clusterfiles_all_particles_type_"+ t +".txt");
    if (up_layer)
    for (const auto &t : considered_types)
        up_files_output[t].open("clusterfiles_up_layer_particles_type_"+ t +".txt");
    if (down_layer)
    for (const auto &t : considered_types)
        down_files_output[t].open("clusterfiles_down_layer_particles_type_"+ t +".txt");

    for ( std::string &file : input_files) {
        std::filesystem::path path(file);
        if (!std::filesystem::exists(path)) {
            fprintf(stderr, "File not found: %s! skipping...\n", path.string().c_str());
            continue;
        }

        std::cout<<"System: " << path.string() << std::endl;
        std::cout<<"Cut-off: " << cluster_cutoff << std::endl;
        std::cout<<"Types: ";
        for (const auto &t : considered_types)
            std::cout << t << " ";
        std::cout<<std::endl;

        float *x, *y, *z;
        float *vx, *vy, *vz;
        char **types;
        bond *bonds;
        int n_particles, n_bonds;
        float lx = 0, ly = 0, lz = 0, xy = 0, xz = 0, yz = 0;
        std::string xmlfilename = path.string();
        
        if (parse_hoomd_xml(path.string().c_str(), 
                            &x, &y, &z, 
                            &vx, &vy, &vz, 
                            &types, 
                            &n_particles,
                            &bonds, &n_bonds,
                            &lx, &ly, &lz, &xy, &xz, &yz) != 0) {
            fprintf(stderr, "Error during parsing: %s! skipping...\n", path.string().c_str());
            continue;
        }   
        int n_types = considered_types.size();

        std::map<std::string,std::vector<float>> x_map, y_map, z_map, x_up, y_up, z_up, x_down, y_down, z_down;
        std::map<std::string,std::vector<int>> andx_map, andx_up, andx_down;

        for (int i = 0; i < n_particles; i++) {
            std::string type = types[i];
            x_map[types[i]].push_back(x[i]);
            y_map[types[i]].push_back(y[i]);
            z_map[types[i]].push_back(z[i]);
            andx_map[types[i]].push_back(i);
        } 

        if (!all)
        for (int ib=0; ib < n_bonds; ib++) {
            std::string typei = bonds[ib].typei, typej = bonds[ib].typej;
            int ai = bonds[ib].ai, aj = bonds[ib].aj, hndx = 0;
            std::string htype;
            float xhead, yhead, zhead, xtail, ytail, ztail;
            if (std::find(considered_types.begin(), considered_types.end(),typei) != considered_types.end()) {
                xhead = x[ai]; yhead = y[ai]; zhead = z[ai]; hndx = ai; htype = typei;
                xtail = x[aj]; ytail = y[aj]; ztail = z[aj];
            } else if (std::find(considered_types.begin(), considered_types.end(),typej) != considered_types.end()) {
                xhead = x[aj]; yhead = y[aj]; zhead = z[aj]; hndx = aj; htype = typej;
                xtail = x[ai]; ytail = y[ai]; ztail = z[ai];
            } else {
                continue;
            }
            float dx = xhead - xtail, dy = yhead - ytail, dz = zhead - ztail;
            float dx_dot_rhead = dx * xhead + dy * yhead + dz * zhead;
            if (dx_dot_rhead > 0 && up_layer) {
                x_up[htype].push_back(xhead); y_up[htype].push_back(yhead); z_up[htype].push_back(zhead);
                andx_up[htype].push_back(hndx);
            } else if ( down_layer ) {
                x_down[htype].push_back(xhead); y_down[htype].push_back(yhead); z_down[htype].push_back(zhead);
                andx_down[htype].push_back(hndx);
            }

        }
        // You can convert float** -> double** if needed for clustering
        // Example usage of `neighboring()` goes here if positions are passed

        printf("Parsed %d particles.\n", n_particles);
        printf("Parsed %d bonds.\n", n_bonds);

        path.replace_extension().filename();

        for (int i = 0; i < n_types; i++) {
            std::string ptype = considered_types[i];
            std::cout << "Type " << ptype ;
            if ( all )
            std::cout << " all: "<< x_map[ptype].size();
            if ( up_layer )
            std::cout << " upper layer: "<< x_up[ptype].size();
            if ( down_layer )
            std::cout << "  down layer: "<< x_down[ptype].size();
            std::cout << std::endl;
            
            std::string filename;
            if (all) {
                filename = path.filename().string() + "_type_" + ptype + "_neighboring.txt";
                neighboring_particles(x_map[ptype].data(), y_map[ptype].data(), z_map[ptype].data(), cluster_cutoff, x_map[ptype].size(), 32, lx, ly, lz, andx_map[ptype].data(),filename.c_str());
                all_files_output[ptype] << filename <<'\t'<<xmlfilename<<'\t'<<ptype<< "all" <<'\t'<<'\n';
            }
            
            if ( up_layer ) {
                filename = path.filename().string() + "_up_type_" + ptype + "_neighboring.txt";
                neighboring_particles(x_up[ptype].data(), y_up[ptype].data(), z_up[ptype].data(), cluster_cutoff, x_up[ptype].size(), 32, lx, ly, lz, andx_up[ptype].data(),filename.c_str());
                up_files_output[ptype] << filename <<'\t'<<xmlfilename<<'\t'<<ptype<<'\t'<< "up" <<'\n';
            }
            
            if ( down_layer ) {
                filename = path.filename().string() + "_down_type_" + ptype + "_neighboring.txt";
                neighboring_particles(x_down[ptype].data(), y_down[ptype].data(), z_down[ptype].data(), cluster_cutoff, x_down[ptype].size(), 32, lx, ly, lz, andx_down[ptype].data(),filename.c_str());
                down_files_output[ptype] << filename <<'\t'<<xmlfilename<<'\t'<<ptype<<'\t'<< "down" << '\n';
            }
        }
        for (int i = 0; i < n_particles; i++) free(types[i]);
        free(x); free(y); free(z);
        free(vx); free(vy); free(vz);
        free(types);
    }
    
    if (all)
    for (const auto &t : considered_types)
        all_files_output[t].close();
    if (up_layer)
    for (const auto &t : considered_types)
        up_files_output[t].close();
    if (down_layer)
    for (const auto &t : considered_types)
        down_files_output[t].close();

    return 0;
}
