#ifndef PLYWRITER_H
#define PLYWRITER_H

/**
 * Static class for writing PLY files.
 *
 * @author Tao Ju
 */

class PLYWriter {
public:
    /// Constructor
    PLYWriter(){};

    /// Write ply header
    static void writeHeader(FILE* fout, int numVert, int numFace) {
        // Ply
        fprintf(fout, "ply\n");

        // Always big endian
        fprintf(fout, "format binary_big_endian 1.0\n");

        // vertex properties
        fprintf(fout, "element vertex %d\n", numVert);
        fprintf(fout, "property float x\n");
        fprintf(fout, "property float y\n");
        fprintf(fout, "property float z\n");

        // face properties
        fprintf(fout, "element face %d\n", numFace);
        fprintf(fout, "property list uchar int vertex_indices\n");

        // End
        fprintf(fout, "end_header\n");
    };

    /// Write vertex
    static void writeVertex(FILE* fout, float vt[3]) {
        float nvt[3];
        for (int i = 0; i < 3; i++) {
            nvt[i] = vt[i];
            flipBits32(&(nvt[i]));
        }
        fwrite(nvt, sizeof(float), 3, fout);
    };

    /// Write face
    static void writeFace(FILE* fout, int num, int fc[]) {
        unsigned char cnum = num;
        fwrite(&cnum, sizeof(unsigned char), 1, fout);
        for (int i = 0; i < num; i++) {
            flipBits32(&(fc[i]));
        }
        fwrite(fc, sizeof(int), num, fout);
    };

    static void flipBits32(void* x) {
        unsigned char* temp = (unsigned char*)x;
        unsigned char swap;

        swap = temp[0];
        temp[0] = temp[3];
        temp[3] = swap;

        swap = temp[1];
        temp[1] = temp[2];
        temp[2] = swap;
    };
};

#endif