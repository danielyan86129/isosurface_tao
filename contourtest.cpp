/************************************************************************
* Just a simple example how to use the contouring code to extract
dose/structure meshes
* Note: files cycle8.txt and pretriang_?.txt need to be in the same folder
*
* Tao Ju (2006)
************************************************************************/

#include "UniGrid.h"
#include "reader.h"
#include <stdio.h>
#include <stdlib.h>

void testDose(char* fname) {
    // Here shows how to extract meshes for a sequence of iso-values from a dose
    // volume

    /////////////////////// Initialize coordinates array ///////////////////////
    float *xoff, *yoff, *zoff;

    // First, read a binary dose file and build a Volume object
    printf("Reading Dose file...\n");
    DoseReader* reader = new DoseReader(fname);
    Volume* vol = reader->getVolume();
    /////////////////////// Read coordinates ///////////////////////
    reader->getCoordinates(xoff, yoff, zoff);

    // Next, setting things up for contouring
    /////////////////////// Use correct coordinates ///////////////////////
    UniGrid* grid = new UniGrid(vol, xoff, yoff, zoff);
    // Here is the optimization for extracting multi-contours from a single
    // volume
    grid->buildTree();

    // Given some iso-values (here is an example of 20 iso-values)
    int numiso = 13;
    float mindose = 0.5f;
    float maxdose = 70;

    // Now, contouring
    for (int i = 0; i < numiso; i++) {
        float isovalue = mindose + (maxdose - mindose) * i / (numiso - 1);

        // Build contour
        grid->genContour(isovalue);

        // Apply several rounds of smoothing, if you need it
        grid->smooth(5);

        // Retrieve the mesh
        int numvertices, numtriangles;
        float* vertices; // a flat list, e.g., y coordinate of ith vertex is
                         // stored at vertices[3*i+1]
        int* triangles;  // again a flat list
        grid->getMesh(numvertices, numtriangles, vertices, triangles);

        // remember to dispose off stuff when you are done with it
        delete vertices;
        delete triangles;
    }

    // if you want to write out to a mesh file
    grid->writePLY("test.ply");

    // clean up
    delete grid;
    delete vol;
    delete reader;
}

void testStructure(char* fname) {
    // Here shows how to extract meshes for all structures

    ///////////////////////Initialize coordinates array/////////////////////
    float *xoff, *yoff, *zoff;

    // First, read a binary structure file (containging all structures)
    printf("Reading Structure file...\n");
    StructureReader* reader = new StructureReader(fname);

    // Now, contouring
    for (int i = 0; i < reader->getNumStructures(); i++) {

        // First, get the volume for ith structure
        Volume* vol = reader->getVolume(i);
        /////////////////////// Read coordinates ///////////////////////
        reader->getCoordinates(xoff, yoff, zoff);

        // Next, setting things up for contouring
        /////////////////////// Use correct coordinates ///////////////////////
        UniGrid* grid = new UniGrid(vol, xoff, yoff, zoff);

        // Build contour (always use 0.5 as isovalue for structures)
        grid->genContour(0.5);

        // Apply several rounds of smoothing, if you need it
        grid->smooth(5);

        // Retrieve the mesh
        int numvertices, numtriangles;
        float* vertices; // a flat list, e.g., y coordinate of ith vertex is
                         // stored at vertices[3*i+1]
        int* triangles;  // again a flat list
        grid->getMesh(numvertices, numtriangles, vertices, triangles);

        // remember to dispose off stuff when you are done with it
        delete vertices;
        delete triangles;

        delete grid;
        delete vol;
    }

    // clean up
    delete reader;
}

int main(int args, char* argv[]) {
    if (strstr(argv[1], ".dos") != NULL) {
        testDose(argv[1]);
    } else if (strstr(argv[1], ".str") != NULL) {
        testStructure(argv[1]);
    }
    return 0;
}
