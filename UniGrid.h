#ifndef UNIGRID_H
#define UNIGRID_H

#include "Cube.h"
#include "HashMap.h"
#include "LinkedList.h"
#include "PLYWriter.h"
#include "TriangulationTree.h"
#include "curvature.h"
#include "geocommon.h"
#include "volume.h"
#include <math.h>
#include <stdio.h>
#include <time.h>

#define G_PI 3.1415926535897932385

struct UniTreeNode {
    float min, max;
    UniTreeNode** chd;
};

class UniGrid {
public:
    // Data
    int dim[3];
    float dis[3];
    float* offset[3];
    float*** data;
    int d1, d2;

    float* verts;
    int totalverts;
    int* tris;
    int totaltris;

    // Octree structure for fast contouring
    UniTreeNode* root;

    // Contouring tables
    TriangulationTree** ttrees;
    Cube* cube;

    // Timings
    clock_t meshTime, searchTime;

    /// Constructor: from a CT file
    UniGrid(const char* fname);

    /// Constructor: from an existing Volume object
    UniGrid(Volume* vol);
    UniGrid(Volume* vol, float spacing);
    UniGrid(Volume* vol, float* xoff, float* yoff, float* zoff);

    /// Destructor
    ~UniGrid();

    /// Initialization
    void init();

    /// Contouring given a iso-value
    void genContour(float iso);
    void genContour(float iso, int isConvex);
    inline void genContourCell(float iso, int isConvex, int i, int j, int k,
                               LinkedList<Point3f>* vertlist,
                               LinkedList<Triangle3i>* trianglist,
                               HashMap* myhash, Status* status);
    inline float getData(int i, int j, int k);

    /// Contouring acceleration
    void buildTree();
    void clearTree(UniTreeNode* node);
    UniTreeNode* buildTree(int x, int y, int z, int len);
    inline void genContourTree(UniTreeNode* node, int x, int y, int z, int len,
                               float iso, int isConvex,
                               LinkedList<Point3f>* vertlist,
                               LinkedList<Triangle3i>* trianglist,
                               HashMap* myhash, Status* status);

    /// Laplacian smoothing of the mesh
    void smooth(int iters);

    /// Retrieve the contour
    void getMesh(int& numverts, int& numtris, float*& vertices,
                 int*& triangles);

    /// Writing procedures
    void writePLY(char* fname);
    void writeMultiPLY(char* fname, int minsize);
    void writeOBJ(char* filename);
    // void writeCurvatures( char* fname );

    struct GridShift {
        double x;
        double y;
        double z;
    };

    GridShift m_GridShift;
};

#endif
