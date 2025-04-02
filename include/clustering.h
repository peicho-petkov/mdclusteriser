#ifndef CLUSTERING_H
#define CLUSTERING_H

// Modified Clustering Algorithm for Molecular Simulation
// By: Fellipe Carvalho de Oliveira - COPPE/PEQ/UFRJ
void neighboring(float **rx, float **ry, float **rz,
                 float dist_cluster, int n_molecules,
                 int max_contacts,  float Lx, float Ly, float Lz, int n_parti_per_molecule, int* aindex, const char *out_name);

// The same algorithm like neighboring() but for particles (instead of molecules)
void neighboring_particles(float *rx, float *ry, float *rz,
                                float dist_cluster, int n_particles,
                                int max_contacts, float Lx, float Ly, float Lz, int *aindex, const char *out_name);

#endif // CLUSTERING_H
