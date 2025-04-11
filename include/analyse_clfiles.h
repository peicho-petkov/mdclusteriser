#ifndef _analyse_clusters_h
#define _analyse_clusters_h

#include <string>


void analyzeClusters(const std::string &inputFile,
    const std::string &histogramFile,
    const std::string &largestClusterFile,
    int &Nclusters,
    int &minSize,
    int &maxSize,
    double &averageSize);

#endif // _analyse_clusters_h