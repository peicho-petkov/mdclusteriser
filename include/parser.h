#ifndef PARSER_H
#define PARSER_H

struct sbond {
    int ai, aj;
    char typei[9], typej[9];
};

typedef sbond bond;

int parse_hoomd_xml(const char *filename,
    float **x, float **y, float **z,
    float **vx, float **vy, float **vz,
    char ***types, int *n_particles,
    bond **bonds, int *n_bonds,
    float *lx, float *ly, float *lz,
    float *xy, float *xz, float *yz);

#endif // PARSER_H
