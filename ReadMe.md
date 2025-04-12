# HOOMD-blue Cluster Analyzer

This project is a C-based toolset for parsing HOOMD-blue 2.9 XML files and performing molecular cluster analysis using a modified distance-based clustering algorithm [Fellipe C. de Oliveira , Shaghayegh Khani , João M. Maia & Frederico W.
Tavares (2020) Modified clustering algorithm for molecular simulation, Molecular Simulation, 46:18, 1453-1466, DOI: 10.1080/08927022.2020.1839661].

## 🔍 Features

- Parses HOOMD-blue 2.9 `.xml` system configuration files
- Extracts:
  - Particle positions
  - Particle velocities
  - Particle types
  - Bonds
- Stores position and velocity components in separate arrays
- Performs clustering based on inter-particle distances
- Outputs cluster information and composition to a file (`clustering.out`)
- Geometric center can be sutracted from the particles position (the term COM (center of mass) in the code)

---

## ⚙️ Build Instructions

### Requirements

- C compiler (GCC or Clang)
- libxml2 (`sudo apt install libxml2-dev` on Ubuntu)
- CMake ≥ 3.10

### Build

```bash
mkdir build
cd build
cmake ..
make
```

### Run

```bash
hoomd_cluster <path/to/hoomd_xml> --cut <float> --types <str> ... <str>
```

Results will be saved in `clustering.out`.

---

## 📦 Output Format

The output file `clustering.out` includes:

* Total number of links (connections)
* Number of clusters
* Number of iterations for convergence
* Molecule indices grouped by cluster

---

## 📚 Example XML Format

```xml
<hoomd_xml version="1.0">
  <configuration time_step="0" dimensions="3">
    <box lx="20.0" ly="20.0" lz="20.0"/>
    <position>
      0.0 0.0 0.0
      1.0 0.0 0.0
    </position>
    <velocity>
      0.0 0.0 0.0
      0.1 0.0 0.0
    </velocity>
    <type>
      A
      B
    </type>
  </configuration>
</hoomd_xml>
```

---

## 🤝 Acknowledgments

This project was created with assistance from ChatGPT ([https://openai.com/chatgpt](https://openai.com/chatgpt)), which helped with:

* Translating clustering algorithms from Fortran to C/C++
* Designing data structures for particle buffers
* Implementing XML parsing using `libxml2`
* Writing the CMake build system
* Generating parts of this README 😊
