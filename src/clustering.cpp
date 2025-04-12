#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

void clustering(int **node_next, int *n_contacts_per_molecule,
                int n_links, int n_molecules, int max_contacts,  int *aindex, const char *out_name)
{
    int *nodeL = (int *)calloc(n_molecules, sizeof(int));
    int *labels = (int *)calloc(max_contacts, sizeof(int));
    int **clusters = (int **)malloc(n_molecules * sizeof(int *));
    int *id = (int *)calloc(n_molecules, sizeof(int));
    int *n_mol_per_cluster = (int *)calloc(n_molecules, sizeof(int));
    for (int i = 0; i < n_molecules; i++) {
        clusters[i] = (int *)calloc(n_molecules, sizeof(int));
        for (int j = 0; j < n_molecules; j++)
            clusters[i][j] = -1;
        nodeL[i] = -1;
        id[i] = -1;
    }

    int n_clusters = 0;
    int tol = 1000, N = 1;
    bool condition1 = true;

    // Initial labeling
    for (int i = 0; i < n_molecules; i++)
    {
        if (n_contacts_per_molecule[i] == 0)
        {
            nodeL[i] = n_clusters++;
        }
        else
        {
            for (int j = 0; j < max_contacts; j++)
                labels[j] = -1;
            for (int j = 0; j < n_contacts_per_molecule[i]; j++)
                labels[j] = nodeL[node_next[i][j]];

            int max_val = -1;
            for (int j = 0; j < n_contacts_per_molecule[i]; j++)
                if (labels[j] > max_val)
                    max_val = labels[j];

            if (max_val == -1)
            {
                nodeL[i] = n_clusters++;
            }
            else
            {
                int min_val = n_molecules + 1;
                for (int j = 0; j < n_contacts_per_molecule[i]; j++)
                    if (labels[j] > -1 && labels[j] < min_val)
                        min_val = labels[j];
                nodeL[i] = min_val;
                for (int j = 0; j < n_contacts_per_molecule[i]; j++)
                    nodeL[node_next[i][j]] = nodeL[i];
            }
        }
    }

    while (N < tol && condition1)
    {
        condition1 = false;
        for (int i = 0; i < n_molecules; i++)
        {
            if (n_contacts_per_molecule[i] != 0)
            {
                for (int j = 0; j < max_contacts; j++)
                    labels[j] = -1;
                for (int j = 0; j < n_contacts_per_molecule[i]; j++)
                    labels[j] = nodeL[node_next[i][j]];
                int min_val = n_molecules + 1;
                for (int j = 0; j < n_contacts_per_molecule[i]; j++)
                    if (labels[j] > -1 && labels[j] < min_val)
                        min_val = labels[j];
                nodeL[i] = min_val;
                for (int j = 0; j < n_contacts_per_molecule[i]; j++)
                    nodeL[node_next[i][j]] = nodeL[i];
            }
        }

        // Check if labeling converged
        for (int i = 0; i < n_molecules; i++)
        {
            for (int j = 0; j < n_contacts_per_molecule[i]; j++)
                labels[j] = nodeL[node_next[i][j]];
            for (int j = 0; j < n_contacts_per_molecule[i]; j++)
            {
                if (nodeL[i] != labels[j])
                {
                    condition1 = true;
                    N++;
                    goto exit_check;
                }
            }
        }
    exit_check:;
    }

    // Recount clusters
    n_clusters = 0;
    for (int i = 0; i < n_molecules; i++)
    {
        bool condition2 = false, condition3 = true;
        for (int j = i + 1; j < n_molecules; j++)
            if (nodeL[i] == nodeL[j])
                condition2 = true;
        for (int j = i - 1; j >= 0; j--)
            if (nodeL[i] == nodeL[j])
                condition3 = false;

        if (condition2 && condition3)
        {
            id[n_clusters++] = nodeL[i];
        }
        if (!condition2 && condition3)
        {
            id[n_clusters++] = nodeL[i];
        }
    }

    for (int i = 0; i < n_clusters; i++)
    {
        for (int j = 0; j < n_molecules; j++)
        {
            if (nodeL[j] == id[i])
            {
                clusters[i][j] = j;
                n_mol_per_cluster[i]++;
            }
        }
    }

    FILE *f = fopen(out_name, "w");

    if (!f) {
        printf("ERROR: enable to create file %s\n",out_name);
        exit(2);
    }
  
    fprintf(f, "Numbers of Links %d\n", n_links);
    fprintf(f, "Number of clusters %d\n", n_clusters);
    fprintf(f, "Number of iterations for convergence %d\n\n", N);
    for (int i = 0; i < n_clusters; i++)
    {
        fprintf(f, "Cluster : %d\n", i + 1);
        fprintf(f, "Molecules (%d):\n", n_mol_per_cluster[i]);
        for (int j = 0; j < n_molecules; j++)
            if (clusters[i][j] > -1)
                fprintf(f, "%d\n", aindex[clusters[i][j]]);
    }
    fclose(f);

    // Free memory
    free(nodeL);
    free(labels);
    for (int i = 0; i < n_molecules; i++)
        free(clusters[i]);
    free(clusters);
    free(id);
    free(n_mol_per_cluster);
}

void neighboring(float **rx, float **ry, float **rz,
                 float dist_cluster, int n_molecules,
                 int max_contacts, float Lx, float Ly, float Lz, int n_parti_per_molecule,  int *aindex, const char *out_name)
{
    int **node_next = (int **)malloc(n_molecules * sizeof(int *));
    int *n_contacts_per_molecule = (int *)calloc(n_molecules, sizeof(int));
    for (int i = 0; i < n_molecules; i++)
        node_next[i] = (int *)calloc(max_contacts, sizeof(int));

    int n_links = 0;

#pragma omp parallel for schedule(dynamic)
    for (int m = 0; m < n_molecules; m++)
    {
        for (int n = m + 1; n < n_molecules; n++)
        {
            int found = 0;
            for (int i = 0; i < n_parti_per_molecule && !found; i++)
            {
                for (int j = 0; j < n_parti_per_molecule; j++)
                {
                    float dx = rx[m][i] - rx[n][j];
                    dx -= Lx * rint(dx / Lx);
                    float dy = ry[m][i] - ry[n][j];
                    dy -= Ly * rint(dy / Ly);
                    float dz = rz[m][i] - rz[n][j];
                    dz -= Lz * rint(dz / Lz);
                    float dist2 = dx * dx + dy * dy + dz * dz;

                    if (dist2 < dist_cluster * dist_cluster)
                    {
#pragma omp critical
                        {
                            if (n_contacts_per_molecule[m] >= max_contacts || n_contacts_per_molecule[n] >= max_contacts)
                            {
                                fprintf(stderr, "Too many contacts!\n");
                                exit(EXIT_FAILURE);
                            }
                            node_next[m][n_contacts_per_molecule[m]++] = n;
                            node_next[n][n_contacts_per_molecule[n]++] = m;
                            n_links++;
                        }
                        found = 1;
                        break;
                    }
                }
                if (found)
                    break;
            }
        }
    }

    clustering(node_next, n_contacts_per_molecule, n_links, n_molecules, max_contacts, aindex, out_name);

    for (int i = 0; i < n_molecules; i++)
        free(node_next[i]);
    free(node_next);
    free(n_contacts_per_molecule);
}

void neighboring_particles(float *rx, float *ry, float *rz,
                           float dist_cluster, int n_particles,
                           int max_contacts, float Lx, float Ly, float Lz, 
                           int *aindex, const char *out_name,
                           bool use_pbc)
{
    int **node_next = (int **)malloc(n_particles * sizeof(int *));
    int *n_contacts_per_molecule = (int *)calloc(n_particles, sizeof(int));
    for (int i = 0; i < n_particles; i++)
        node_next[i] = (int *)calloc(max_contacts, sizeof(int));

    int n_links = 0;

    int pbc = 0;
    if (use_pbc) pbc = 1;

#pragma omp parallel for schedule(dynamic)
    for (int m = 0; m < n_particles; m++)
    {
        for (int n = m + 1; n < n_particles; n++)
        {
            int found = 0;
            float dx = rx[m] - rx[n];
            dx -= pbc * Lx * rint(dx / Lx);
            float dy = ry[m] - ry[n];
            dy -= pbc * Ly * rint(dy / Ly);
            float dz = rz[m] - rz[n];
            dz -= pbc * Lz * rint(dz / Lz);
            float dist2 = dx * dx + dy * dy + dz * dz;

            if (dist2 < dist_cluster * dist_cluster)
            {
#pragma omp critical
                {
                    if (n_contacts_per_molecule[m] >= max_contacts || n_contacts_per_molecule[n] >= max_contacts)
                    {
                        fprintf(stderr, "Too many contacts!\n");
                        exit(EXIT_FAILURE);
                    }
                    node_next[m][n_contacts_per_molecule[m]++] = n;
                    node_next[n][n_contacts_per_molecule[n]++] = m;
                    n_links++;
                }
                // found = 1;
                // break;
            }
        }
    }

    clustering(node_next, n_contacts_per_molecule, n_links, n_particles, max_contacts, aindex, out_name);

    for (int i = 0; i < n_particles; i++)
        free(node_next[i]);
    free(node_next);
    free(n_contacts_per_molecule);
}
