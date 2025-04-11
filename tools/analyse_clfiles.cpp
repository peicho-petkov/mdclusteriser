#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept> // for runtime_error
#include <utility>   // for move

// This function:
//  1) Reads clusters from inputFile
//  2) Builds histogram of cluster sizes
//  3) Writes histogram to histogramFile
//  4) Finds and writes the largest cluster’s particle IDs to largestClusterFile
//  5) Returns minSize, maxSize, and averageSize via reference parameters
void analyzeClusters(const std::string &inputFile,
                     const std::string &histogramFile,
                     const std::string &largestClusterFile,
                     int &Nclusters,
                     int &minSize,
                     int &maxSize,
                     double &averageSize)
{
    // Open the input file
    std::ifstream infile(inputFile);
    if (!infile.is_open()) {
        throw std::runtime_error("Failed to open input file: " + inputFile);
    }

    // Parse cluster data into a vector of vectors: clusters[clusterIndex] = list of particle IDs
    std::vector<std::vector<int>> clusters;
    std::string line;
    while (std::getline(infile, line)) {
        const std::string marker = "Molecules (";
        std::size_t pos = line.find(marker);
        if (pos != std::string::npos) {
            // Parse the integer count in "Molecules (X):"
            pos += marker.size(); 
            std::size_t endPos = line.find(')', pos);
            if (endPos != std::string::npos) {
                std::string numberStr = line.substr(pos, endPos - pos);
                int count = 0;
                try {
                    count = std::stoi(numberStr);
                } catch (...) {
                    // Ignore malformed lines
                    continue;
                }

                // Read 'count' lines for particle IDs
                std::vector<int> particleIDs;
                particleIDs.reserve(count);
                for (int i = 0; i < count; ++i) {
                    if (!std::getline(infile, line)) {
                        // If file ends unexpectedly, break out
                        break;
                    }
                    try {
                        int id = std::stoi(line);
                        particleIDs.push_back(id);
                    } catch (...) {
                        // Ignore parse errors
                    }
                }
                clusters.push_back(std::move(particleIDs));
            }
        }
    }
    infile.close();

    // Check if we actually found any clusters
    // if (clusters.empty()) {
    //     throw std::runtime_error("No clusters found in the input file.");
    // }
    Nclusters = static_cast<int>(clusters.size());
    if (Nclusters == 0) {
        return;
    }
    // Build histogram, track min, max, sum
    std::map<int,int> histogram; // clusterSize -> count
    minSize = static_cast<int>(clusters[0].size());
    maxSize = static_cast<int>(clusters[0].size());
    long long totalSize = 0; // to handle potentially large sums
    size_t minIndex = 0;
    size_t maxIndex = 0;

    // Compute stats
    for (size_t i = 0; i < clusters.size(); ++i) {
        int sz = static_cast<int>(clusters[i].size());
        histogram[sz]++;
        
        if (sz < minSize) {
            minSize = sz;
            minIndex = i;
        }
        if (sz > maxSize) {
            maxSize = sz;
            maxIndex = i;
        }
        totalSize += sz;
    }

    averageSize = static_cast<double>(totalSize) / static_cast<double>(clusters.size());

    // Write the histogram to histogramFile
    {
        std::ofstream hfile(histogramFile);
        if (!hfile.is_open()) {
            throw std::runtime_error("Failed to open histogram file: " + histogramFile);
        }

        hfile << "Cluster size histogram (size count)\n";
        for (const auto &kv : histogram) {
            hfile << kv.first << " " << kv.second << "\n";
        }
    }

    // Write the largest cluster’s particle IDs to largestClusterFile
    {
        std::ofstream lcfile(largestClusterFile);
        if (!lcfile.is_open()) {
            throw std::runtime_error("Failed to open largest-cluster file: " + largestClusterFile);
        }

        lcfile << "Largest cluster size: " << maxSize << "\n";
        lcfile << "Particle IDs in the largest cluster:\n";
        for (int id : clusters[maxIndex]) {
            lcfile << id << "\n";
        }
    }
}

// // Example main() to demonstrate usage of analyzeClusters()
// int main() {
//     // Input file
//     const std::string inputFile = "snapshot.0326000000_down_type_A_neighboring.txt";
//     // Output files
//     const std::string histogramFile = "cluster_histogram.txt";
//     const std::string largestClusterFile = "largest_cluster_ids.txt";

//     int Nclusters = 0;  // Number of clusters found in the input file
//     int minSize = 0;
//     int maxSize = 0;
//     double averageSize = 0.0;

//     try {
//         analyzeClusters(inputFile, histogramFile, largestClusterFile,
//                         Nclusters, minSize, maxSize, averageSize);

//         // Now minSize, maxSize, and averageSize are filled in
//         std::cout << "Analysis complete.\n";
//         std::cout << "Number of clusters = " << Nclusters << "\n";
//         std::cout << "Min cluster size = " << minSize << "\n";
//         std::cout << "Max cluster size = " << maxSize << "\n";
//         std::cout << "Average cluster size = " << averageSize << "\n";
//         std::cout << "Histogram saved to: " << histogramFile << "\n";
//         std::cout << "Largest cluster IDs saved to: " << largestClusterFile << "\n";
//     }
//     catch (const std::exception &ex) {
//         std::cerr << "ERROR: " << ex.what() << "\n";
//         return 1;
//     }

//     return 0;
// }
