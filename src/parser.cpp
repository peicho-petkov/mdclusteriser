#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "parser.h"

int count_lines(const char *text)
{
    int count = 0;
    const char *p = text;
    while (*p)
    {
        if (*p == '\n')
            count++;
        p++;
    }
    return count;
}

void parse_vector_components(xmlNode *node, float *x, float *y, float *z, int n_particles)
{
    xmlChar *content = xmlNodeGetContent(node);
    char *line = strtok((char *)content, "\n");
    int idx = 0;

    while (line && idx < n_particles)
    {
        float px, py, pz;
        if (sscanf(line, "%f %f %f", &px, &py, &pz) == 3)
        {
            x[idx] = px;
            y[idx] = py;
            z[idx] = pz;
            idx++;
        }
        line = strtok(NULL, "\n");
    }

    xmlFree(content);
}

void parse_bond_components(xmlNode *node, bond *bonds, int n_bonds)
{
    xmlChar *content = xmlNodeGetContent(node);
    char *line = strtok((char *)content, "\n");
    int idx = 0;

    while (line && idx < n_bonds)
    {
        int a1, a2;
        char type1[9], type2[9];
        int len = strlen(line);
        for (int i = 0; i < len; i++) {
            if (line[i] == '-') {
                line[i] = ' ';
                break;
            }
        }
        if (sscanf(line, "%s %s %d %d", type1, type2, &a1, &a2) == 4)
        {
            bonds[idx].ai = a1;
            bonds[idx].aj = a2;
            strcpy(bonds[idx].typei, type1);
            strcpy(bonds[idx].typej, type2);
            idx++;
        }
        line = strtok(NULL, "\n");
    }

    xmlFree(content);
}

void parse_type_block(xmlNode *node, char **types, int n_particles)
{
    xmlChar *content = xmlNodeGetContent(node);
    char *line = strtok((char *)content, "\n");
    int idx = 0;

    while (line && idx < n_particles)
    {
        types[idx] = strdup(line);
        idx++;
        line = strtok(NULL, "\n");
    }

    xmlFree(content);
}

// int parse_hoomd_xml(const char *filename,
//                     float **x, float **y, float **z,
//                     float **vx, float **vy, float **vz,
//                     char ***types, int *n_particles)
int parse_hoomd_xml(const char *filename,
                    float **x, float **y, float **z,
                    float **vx, float **vy, float **vz,
                    char ***types, int *n_particles,
                    bond **bonds, int *n_bonds,
                    float *lx, float *ly, float *lz,
                    float *xy, float *xz, float *yz)
{
    xmlDoc *doc = xmlReadFile(filename, NULL, 0);
    if (!doc)
    {
        fprintf(stderr, "Could not parse file %s\n", filename);
        return 1;
    }

    xmlNode *root = xmlDocGetRootElement(doc);
    xmlNode *conf = root->children;
    xmlNode *position_node = NULL, *velocity_node = NULL, *type_node = NULL;
    xmlNode *box_node = NULL, *bond_node = NULL;
    
    for (; conf; conf = conf->next)
    {
        if (conf->type == XML_ELEMENT_NODE &&
            strcmp((char *)conf->name, "configuration") == 0)
        {
            for (xmlNode *child = conf->children; child; child = child->next)
            {
                if (child->type != XML_ELEMENT_NODE)
                    continue;
                if (!strcmp((char *)child->name, "position"))
                    position_node = child;
                else if (!strcmp((char *)child->name, "velocity"))
                    velocity_node = child;
                else if (!strcmp((char *)child->name, "type"))
                    type_node = child;
                else if (!strcmp((char *)child->name, "box"))
                    box_node = child;
                else if (!strcmp((char *)child->name, "bond"))
                    bond_node = child;
            }
        }
    }

    if (box_node)
    {
        xmlChar *val;

        val = xmlGetProp(box_node, (const xmlChar *)"lx");
        if (val)
        {
            *lx = atof((const char *)val);
            xmlFree(val);
        }

        val = xmlGetProp(box_node, (const xmlChar *)"ly");
        if (val)
        {
            *ly = atof((const char *)val);
            xmlFree(val);
        }

        val = xmlGetProp(box_node, (const xmlChar *)"lz");
        if (val)
        {
            *lz = atof((const char *)val);
            xmlFree(val);
        }

        val = xmlGetProp(box_node, (const xmlChar *)"xy");
        if (val)
        {
            *xy = atof((const char *)val);
            xmlFree(val);
        }

        val = xmlGetProp(box_node, (const xmlChar *)"xz");
        if (val)
        {
            *xz = atof((const char *)val);
            xmlFree(val);
        }

        val = xmlGetProp(box_node, (const xmlChar *)"yz");
        if (val)
        {
            *yz = atof((const char *)val);
            xmlFree(val);
        }
    }

    if (!position_node || !type_node)
    {
        fprintf(stderr, "Missing required <position> or <type> blocks\n");
        xmlFreeDoc(doc);
        return 2;
    }

    xmlChar *raw = xmlNodeGetContent(position_node);
    *n_particles = count_lines((char *)raw)-1;
    xmlFree(raw);
    
    raw = xmlNodeGetContent(bond_node);
    *n_bonds = count_lines((char *)raw)-1;
    xmlFree(raw);

    *x = (float*) malloc(sizeof(float) * (*n_particles));
    *y = (float*) malloc(sizeof(float) * (*n_particles));
    *z = (float*) malloc(sizeof(float) * (*n_particles));
    *vx = (float*) malloc(sizeof(float) * (*n_particles));
    *vy = (float*) malloc(sizeof(float) * (*n_particles));
    *vz = (float*) malloc(sizeof(float) * (*n_particles));
    *types = (char**) malloc(sizeof(char *) * (*n_particles));

    parse_vector_components(position_node, *x, *y, *z, *n_particles);

    if (velocity_node)
        parse_vector_components(velocity_node, *vx, *vy, *vz, *n_particles);
    else
        memset(*vx, 0, sizeof(float) * (*n_particles)),
            memset(*vy, 0, sizeof(float) * (*n_particles)),
            memset(*vz, 0, sizeof(float) * (*n_particles));

    parse_type_block(type_node, *types, *n_particles);

    if (!bond_node)
    {
        *n_bonds = 0;
        fprintf(stderr, "Missing required <bond> block\n");
        xmlFreeDoc(doc);
    } else {
        *bonds = (bond*) malloc(sizeof(bond) * (*n_bonds));
        parse_bond_components(bond_node, *bonds, *n_bonds);
    }

    xmlFreeDoc(doc);
    xmlCleanupParser();
    return 0;
}
