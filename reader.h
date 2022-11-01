/**
 * Reading routines for volumetric data
 *
 * Author: Tao Ju
 * Date: 02/16/2005
 */

#ifndef READER_H
#define READER_H

//#define CT_STRUCTURE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "volume.h"

class VolumeReader {
public:
    VolumeReader(){};

    /* Read volume from input */
    virtual Volume* getVolume() = 0;

    /* Get resolution */
    virtual void getSpacing(float& ax, float& ay, float& az) = 0;

    void flipBits32(void* x) {
        unsigned char* temp = (unsigned char*)x;
        unsigned char swap;

        swap = temp[0];
        temp[0] = temp[3];
        temp[3] = swap;

        swap = temp[1];
        temp[1] = temp[2];
        temp[2] = swap;
    }
};

class SOFReader : public VolumeReader {
private:
    char soffile[1024];

    /* Recursive reader */
    void readSOF(FILE* fin, Volume* vol, int off[3], int len) {
        // printf("%d %d %d: %d\n", off[0], off[1], off[2], len) ;
        char type;
        int noff[3];
        int nlen = len / 2;

        // Get type
        fread(&type, sizeof(char), 1, fin);

        if (type == 0) {
            // Internal node
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    for (int k = 0; k < 2; k++) {
                        noff[0] = off[0] + i * nlen;
                        noff[1] = off[1] + j * nlen;
                        noff[2] = off[2] + k * nlen;
                        readSOF(fin, vol, noff, nlen);
                    }
        } else if (type == 1) {
            // Empty node
            char sg;
            fread(&sg, sizeof(char), 1, fin);

            for (int i = 0; i <= len; i++)
                for (int j = 0; j <= len; j++)
                    for (int k = 0; k <= len; k++) {
                        noff[0] = off[0] + i;
                        noff[1] = off[1] + j;
                        noff[2] = off[2] + k;
                        vol->setDataAt(noff[0], noff[1], noff[2], -sg);
                    }
        } else if (type == 2) {
            // Leaf node
            char sg;
            fread(&sg, sizeof(char), 1, fin);

            int t = 0;
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    for (int k = 0; k < 2; k++) {
                        noff[0] = off[0] + i;
                        noff[1] = off[1] + j;
                        noff[2] = off[2] + k;
                        vol->setDataAt(noff[0], noff[1], noff[2],
                                       -((sg >> t) & 1));
                        t++;
                    }
        } else {
            printf("Wrong!\n");
        }
    }

public:
    /* Initializer */
    SOFReader(char* fname) {
        sprintf(soffile, "%s", fname);
        FILE* fin = fopen(fname, "rb");

        if (fin == NULL) {
            printf("Can not open file %s.\n", fname);
        }

        fclose(fin);
    }

    /* Read volume */
    Volume* getVolume() {
        int sx, sy, sz;

        FILE* fin = fopen(soffile, "rb");

        // Process header
        fread(&sx, sizeof(int), 1, fin);
        sy = sx;
        sz = sx;
        printf("Dimensions: %d %d %d\n", sx, sy, sz);

        Volume* rvalue = new Volume(sx + 1, sy + 1, sz + 1);

        // Recursive reader
        int off[3] = {0, 0, 0};
        readSOF(fin, rvalue, off, sx);

        printf("Done reading.\n");
        fclose(fin);
        return rvalue;
    }

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = 1;
        ay = 1;
        az = 1;
    }
};

class MRCReader : public VolumeReader {
public:
    /* Initializer */
    MRCReader(char* fname) {
        sprintf(mrcfile, "%s", fname);

        FILE* fin = fopen(fname, "rb");

        // Parse header
        fread(&totx, sizeof(int), 1, fin);
        fread(&toty, sizeof(int), 1, fin);
        fread(&totz, sizeof(int), 1, fin);

        fread(&mode, sizeof(int), 1, fin);

        fread(&offx, sizeof(int), 1, fin);
        fread(&offy, sizeof(int), 1, fin);
        fread(&offz, sizeof(int), 1, fin);

        fread(&dimx, sizeof(int), 1, fin);
        fread(&dimy, sizeof(int), 1, fin);
        fread(&dimz, sizeof(int), 1, fin);
        dimx++;
        dimy++;
        dimz++;

        fread(&angsx, sizeof(float), 1, fin);
        fread(&angsy, sizeof(float), 1, fin);
        fread(&angsz, sizeof(float), 1, fin);

        fread(&anglex, sizeof(float), 1, fin);
        fread(&angley, sizeof(float), 1, fin);
        fread(&anglez, sizeof(float), 1, fin);

        fseek(fin, 12, SEEK_CUR);

        fread(&dmin, sizeof(float), 1, fin);
        fread(&dmax, sizeof(float), 1, fin);
        fread(&dmean, sizeof(float), 1, fin);

        fseek(fin, 4 * 32, SEEK_CUR);

        fread(&drms, sizeof(float), 1, fin);
        fclose(fin);

        dimx = totx;
        dimy = toty;
        dimz = totz;

        printf("Dimension: %d %d %d\n", dimx, dimy, dimz);
        printf("Mode: %d\n", mode);
        printf("Density: from %f to %f, mean at %f, rms at %f\n", dmin, dmax,
               dmean, drms);
        printf("Cell size: %f %f %f\n", angsx / (dimx - 1), angsy / (dimy - 1),
               angsz / (dimz - 1));
        printf("Cell angles: %f %f %f\n", anglex, angley, anglez);

        if (mode > 2) {
            printf("Complex mode not supported.\n");
            exit(0);
        }

        /* Hacking code
        fseek( fin, 0, SEEK_END ) ;
        long len = ftell( fin ) ;
        len -= 1024 ;

        dimen = 1 ;
        while ( dimen * dimen * dimen < len / 4 )
        {
                dimen ++ ;
        }
        printf("Size: %d\n", dimen) ;
        */
    }

    /* Read volume */
    Volume* getVolume() {
        FILE* fin = fopen(mrcfile, "rb");
        fseek(fin, 1024, SEEK_SET);

        char chard;
        short shortd;
        float floatd;
        double d;

        Volume* vol = new Volume(dimx, dimy, dimz);
        for (int i = 0; i < dimz; i++)
            for (int j = 0; j < dimy; j++)
                for (int k = 0; k < dimx; k++) {
                    switch (mode) {
                        case 0:
                            fread(&chard, sizeof(char), 1, fin);
                            d = (double)chard;
                            break;
                        case 1:
                            fread(&shortd, sizeof(short), 1, fin);
                            d = (double)shortd;
                            break;
                        case 2:
                            fread(&floatd, sizeof(float), 1, fin);
                            d = (double)floatd;
                            break;
                    }

                    vol->setDataAt(i, j, k, d);
                }
        fclose(fin);

        return vol;
    }

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = angsx / (dimx - 1);
        ay = angsy / (dimy - 1);
        az = angsz / (dimz - 1);
    }

private:
    int totx, toty, totz;
    int offx, offy, offz;
    int dimx, dimy, dimz;

    float angsx, angsy, angsz;
    float anglex, angley, anglez;
    float dmin, dmax, dmean, drms;

    int mode;

    char mrcfile[1024];
};

class InvMRCReader : public VolumeReader {
public:
    /* Initializer */
    InvMRCReader(char* fname) {
        sprintf(mrcfile, "%s", fname);

        FILE* fin = fopen(fname, "rb");

        // Parse header
        ifread(&totx, sizeof(int), 1, fin);
        ifread(&toty, sizeof(int), 1, fin);
        ifread(&totz, sizeof(int), 1, fin);

        ifread(&mode, sizeof(int), 1, fin);

        ifread(&offx, sizeof(int), 1, fin);
        ifread(&offy, sizeof(int), 1, fin);
        ifread(&offz, sizeof(int), 1, fin);
        printf("%d %d %d\n", offx, offy, offz);

        ifread(&dimx, sizeof(int), 1, fin);
        ifread(&dimy, sizeof(int), 1, fin);
        ifread(&dimz, sizeof(int), 1, fin);

        ifread(&angsx, sizeof(float), 1, fin);
        ifread(&angsy, sizeof(float), 1, fin);
        ifread(&angsz, sizeof(float), 1, fin);

        ifread(&anglex, sizeof(float), 1, fin);
        ifread(&angley, sizeof(float), 1, fin);
        ifread(&anglez, sizeof(float), 1, fin);

        fseek(fin, 12, SEEK_CUR);

        ifread(&dmin, sizeof(float), 1, fin);
        ifread(&dmax, sizeof(float), 1, fin);
        ifread(&dmean, sizeof(float), 1, fin);

        fseek(fin, 4 * 32, SEEK_CUR);

        ifread(&drms, sizeof(float), 1, fin);
        fclose(fin);

        dimx = totx;
        dimy = toty;
        dimz = totz;
        printf("Dimension: %d %d %d\n", dimx, dimy, dimz);
        printf("Mode: %d\n", mode);
        printf("Density: from %f to %f, mean at %f, rms at %f\n", dmin, dmax,
               dmean, drms);
        printf("Cell size: %f %f %f\n", angsx / (dimx - 1), angsy / (dimy - 1),
               angsz / (dimz - 1));
        printf("Cell angles: %f %f %f\n", anglex, angley, anglez);

        if (mode > 2) {
            printf("Complex mode not supported.\n");
            // exit(0) ;
        }

        /* Hacking code
        fseek( fin, 0, SEEK_END ) ;
        long len = ftell( fin ) ;
        len -= 1024 ;

        dimen = 1 ;
        while ( dimen * dimen * dimen < len / 4 )
        {
                dimen ++ ;
        }
        printf("Size: %d\n", dimen) ;
        */
    }

    /* Read volume */
    Volume* getVolume() {
        FILE* fin = fopen(mrcfile, "rb");
        fseek(fin, 1024, SEEK_SET);

        char chard;
        short shortd;
        float floatd;
        double d;

        Volume* vol = new Volume(dimx, dimy, dimz);
        for (int i = 0; i < dimz; i++)
            for (int j = 0; j < dimy; j++)
                for (int k = 0; k < dimx; k++) {
                    switch (mode) {
                        case 0:
                            fread(&chard, sizeof(char), 1, fin);
                            d = (double)chard;
                            break;
                        case 1:
                            fread(&shortd, sizeof(short), 1, fin);
                            d = (double)shortd;
                            break;
                        default:
                            ifread(&floatd, sizeof(float), 1, fin);
                            d = (double)floatd;
                            break;
                    }

                    //					printf("%g\n", d)
                    //;exit(0)
                    //;

                    vol->setDataAt(i, j, k, d);
                }
        fclose(fin);

        return vol;
    }

    void ifread(void* d, size_t s1, size_t s2, FILE* f) {
        fread(d, s1, s2, f);
        flipBits32(d);
    }

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = angsx / (dimx - 1);
        ay = angsy / (dimy - 1);
        az = angsz / (dimz - 1);
    }

private:
    int totx, toty, totz;
    int offx, offy, offz;
    int dimx, dimy, dimz;

    float angsx, angsy, angsz;
    float anglex, angley, anglez;
    float dmin, dmax, dmean, drms;

    int mode;

    char mrcfile[1024];
};

class MRCReaderPicker {
public:
    MRCReaderPicker(){};

    static VolumeReader* pick(char* fname) {
        FILE* fin = fopen(fname, "rb");
        if (fin == NULL) {
            printf("Error reading MRC file %s.\n", fname);
            exit(0);
        }
        int totx, toty, totz;

        // Parse header
        fread(&totx, sizeof(int), 1, fin);
        fread(&toty, sizeof(int), 1, fin);
        fread(&totz, sizeof(int), 1, fin);

        fclose(fin);

        if (totx <= 0 || totx > 1024) {
            printf("Calling inverse MRCreader.\n");
            return new InvMRCReader(fname);
        } else {
            printf("Calling MRCreader.\n");
            return new MRCReader(fname);
        }
    }
};

class HackMRCReader : public VolumeReader {
public:
    /* Initializer */
    HackMRCReader(char* fname) {
        sprintf(mrcfile, "%s", fname);

        FILE* fin = fopen(fname, "rb");

        fseek(fin, 0, SEEK_END);
        long len = ftell(fin);
        len -= 1024;

        dimen = 1;
        while (dimen * dimen * dimen < len / 4) {
            dimen++;
        }
        printf("Volume Size: %d\n", dimen);

        fclose(fin);
    }

    /* Read volume */
    Volume* getVolume() {
        FILE* fin = fopen(mrcfile, "rb");
        fseek(fin, 1024, SEEK_SET);

        Volume* vol = new Volume(dimen, dimen, dimen);
        float d;
        for (int i = 0; i < dimen; i++)
            for (int j = 0; j < dimen; j++)
                for (int k = 0; k < dimen; k++) {
                    fread(&d, sizeof(float), 1, fin);
                    flipBits32(&d);
                    // printf("%g\n", d) ;exit(0) ;
                    vol->setDataAt(k, j, i, d);
                }
        fclose(fin);

        return vol;
    }

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = 1;
        ay = 1;
        az = 1;
    }

private:
    int dimen;
    char mrcfile[1024];
};

class AnalyzeReader : public VolumeReader {
private:
    char imgfile[1024];
    int dimx, dimy, dimz;
    int datatype;
    float spx, spy, spz;

public:
    /* Initializer */
    AnalyzeReader(char* fname) {
        struct dsr hdr;
        // int size;
        // double cmax, cmin;
        FILE* fp;
        if ((fp = fopen(fname, "r")) == NULL) {
            fprintf(stderr, "Can't open:<%s>\n", fname);
            exit(0);
        }
        fread(&hdr, 1, sizeof(struct dsr), fp);
        fclose(fp);
        if (hdr.dime.dim[0] < 0 || hdr.dime.dim[0] > 15)
            swap_hdr(&hdr);
        struct dsr* hdrp = &hdr;

        dimx = (int)hdrp->dime.dim[1];
        dimy = (int)hdrp->dime.dim[2];
        dimz = (int)hdrp->dime.dim[3];

        spx = hdrp->dime.pixdim[1];
        spy = hdrp->dime.pixdim[2];
        spz = hdrp->dime.pixdim[3];

        datatype = hdrp->dime.datatype;

        strcpy(imgfile, fname);
        strcpy(imgfile + strlen(fname) - 4, ".img");
        printf("Dimensions: %d %d %d, Datatype: %d, ImgFile: %s Spacing: %f %f "
               "%f\n",
               dimx, dimy, dimz, datatype, imgfile, spx, spy, spz);
    }

    /* Read volume */
    Volume* getVolume() {
        FILE* fin = fopen(imgfile, "rb");

        Volume* vol = new Volume(dimx, dimy, dimz);
        float d;
        unsigned char d2;
        short d4;
        int d8;
        for (int k = 0; k < dimz; k++)
            for (int j = 0; j < dimy; j++)
                for (int i = 0; i < dimx; i++) {
                    switch (datatype) {
                        case 2:
                            fread(&d2, sizeof(unsigned char), 1, fin);
                            d = (float)d2;
                            break;
                        case 4:
                            fread(&d4, sizeof(short), 1, fin);
                            d = (float)d4;
                            break;
                        case 8:
                            fread(&d8, sizeof(int), 1, fin);
                            d = (float)d8;
                            break;
                        case 16:
                            fread(&d, sizeof(float), 1, fin);
                            break;
                    }
                    // printf("%g\n", d) ;exit(0) ;
                    vol->setDataAt(i, j, k, d);
                }
        fclose(fin);

        vol->setSpacing(spx, spy, spz);

        return vol;
    }

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = spx;
        ay = spy;
        az = spz;
    }

private:
    int dimen;
    char mrcfile[1024];
};

class TimDoseReader : public VolumeReader {
private:
    char imgfile[1024];
    int dimx, dimy, dimz;
    float spx, spy, spz;

public:
    /* Initializer */
    TimDoseReader(char* fname) {
        strcpy(imgfile, fname);
        FILE* fin;
        if ((fin = fopen(fname, "r")) == NULL) {
            fprintf(stderr, "Can't open:<%s>\n", fname);
            exit(0);
        }

        // Scan
#ifdef CT_STRUCTURE
        {
            // Scan Header
            fscanf(fin, "%d %d %d \n", &dimx, &dimy, &dimz);

            readSpacing(fin);
        }
#else
        {
            // Scan Header
            int xSize, ySize, zSize;
            fscanf(fin, "%d %d %d \n", &xSize, &ySize, &zSize);

            printf("skipping...");
            skipFloats(fin, xSize);
            skipFloats(fin, ySize);
            skipFloats(fin, zSize);
            skipFloats(fin, (xSize * ySize * zSize));
            printf("skipped.");
        }
#endif
        // Structures
        {
            // Structures Header
            int xSize, ySize, zSize, numStructures;
            fscanf(fin, "%d %d %d %d \n", &xSize, &ySize, &zSize,
                   &numStructures);

            skipFloats(fin, xSize);
            skipFloats(fin, ySize);
            skipFloats(fin, zSize);
            skipFloats(fin, (xSize * ySize * zSize));
        }

#ifndef CT_STRUCTURE
        // Dose
        {
            // Dose Header
            fscanf(fin, "%d %d %d \n", &dimx, &dimy, &dimz);
            readSpacing(fin);
        }
#endif
        fclose(fin);

        printf("Dimensions: %d %d %d, ImgFile: %s Spacing: %f %f %f\n", dimx,
               dimy, dimz, imgfile, spx, spy, spz);
    }

    void readSpacing(FILE* fid) {
        // Read dose coords
        spx = 0;
        spy = 0;
        spz = 0;
        float lx = 0, hx = 0, ly = 0, hy = 0, lz = 0, hz = 0;
        int i;
        for (int i = 0; i < dimx; i++) {
            lx = hx;
            fscanf(fid, "%f ", &hx);
            if (i) {
                spx = hx - lx;
            }
        }
        for (int i = 0; i < dimy; i++) {
            ly = hy;
            fscanf(fid, "%f ", &hy);
            if (i) {
                spy = hy - ly;
            }
        }
        for (int i = 0; i < dimz; i++) {
            lz = hz;
            fscanf(fid, "%f ", &hz);
            if (i) {
                spz = hz - lz;
            }
        }
    }

    void skipInts(FILE* fin, int n) {
        skipFloats(fin, n);
        /*
        int temp;
        for(int i = 0; i < n; i++)
        {
                fscanf(fin, "%d ", &temp);
        }
        */
    }

    void skipFloats(FILE* fin, int n) {
        char* line = new char[n * 10];
        fgets(line, n * 10, fin);
        delete line;
        /*
        float temp;
        for(int i = 0; i < n; i++)
        {
                fscanf(fin, "%f ", &temp);
        }
        */
    }

    void readVolume(Volume* vol, FILE* fid) {
        float maxd = -100000, mind = 100000;
        for (int x = 0; x < dimx; x++)
            for (int y = 0; y < dimy; y++)
                for (int z = 0; z < dimz; z++) {
                    float data;
                    fscanf(fid, "%f ", &data);
                    vol->setDataAt(x, y, z, data);

                    if (data < mind) {
                        mind = data;
                    }
                    if (data > maxd) {
                        maxd = data;
                    }
                }

        fclose(fid);
    }

    /* Read volume */
    Volume* getVolume() {
        Volume* vol = new Volume(dimx, dimy, dimz);
        vol->setSpacing(spx, spy, spz);

        FILE* fin;
        if ((fin = fopen(imgfile, "r")) == NULL) {
            fprintf(stderr, "Can't open:<%s>\n", imgfile);
            exit(0);
        }

        // Scan Header
        int xSize, ySize, zSize;
        fscanf(fin, "%d %d %d \n", &xSize, &ySize, &zSize);

        skipFloats(fin, xSize);
        skipFloats(fin, ySize);
        skipFloats(fin, zSize);
        skipFloats(fin, (xSize * ySize * zSize));

        // Structures
        {
            // Structures Header
            int xSize, ySize, zSize, numStructures;
            fscanf(fin, "%d %d %d %d \n", &xSize, &ySize, &zSize,
                   &numStructures);

            skipFloats(fin, xSize);
            skipFloats(fin, ySize);
            skipFloats(fin, zSize);
#ifdef CT_STRUCTURE
            readVolume(vol, fin);
#else
            skipFloats(fin, (xSize * ySize * zSize));
#endif
        }

#ifndef CT_STRUCTURE
        // Dose
        {
            // Dose Header
            int xSize, ySize, zSize;
            fscanf(fin, "%d %d %d \n", &xSize, &ySize, &zSize);
            skipFloats(fin, xSize);
            skipFloats(fin, ySize);
            skipFloats(fin, zSize);
            readVolume(vol, fin);
        }
#endif
        fclose(fin);

        return vol;
    }

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = spx;
        ay = spy;
        az = spz;
    }

private:
    int dimen;
    char mrcfile[1024];
};

class StructureReader : public VolumeReader {
private:
    char imgfile[1024];
    int dimx, dimy, dimz;
    float spx, spy, spz;
    float* xoff;
    float* yoff;
    float* zoff;
    int numstructs;

public:
    /* Initializer */
    StructureReader(const char* fname) {
        strcpy(imgfile, fname);

        FILE* fid;
        if ((fid = fopen(fname, "rb")) == NULL) {
            fprintf(stderr, "Can't open:<%s>\n", fname);
            exit(0);
        }

        // Read structure Header
        fread(&dimx, sizeof(int), 1, fid);
        fread(&dimy, sizeof(int), 1, fid);
        fread(&dimz, sizeof(int), 1, fid);
        fread(&numstructs, sizeof(int), 1, fid);
        // Read structure coords
        xoff = new float[dimx];
        fread(xoff, sizeof(float), dimx, fid);
        spx = xoff[dimx - 1] - xoff[dimx - 2];

        yoff = new float[dimy];
        fread(yoff, sizeof(float), dimy, fid);
        spy = yoff[dimy - 1] - yoff[dimy - 2];

        zoff = new float[dimz];
        fread(zoff, sizeof(float), dimz, fid);
        spz = zoff[dimz - 1] - zoff[dimz - 2];

        fclose(fid);

        printf("Dimensions: %d %d %d, Spacing: %f %f %f, structures: %d\n",
               dimx, dimy, dimz, spx, spy, spz, numstructs);
        printf("Origin: %f %f %f\n", xoff[0], yoff[0], zoff[0]);
    }

    void skipInts(FILE* fin, int n) {
        int* data = new int[n];
        fread(data, sizeof(int), n, fin);
        delete data;
    }

    void skipFloats(FILE* fin, int n) {
        float* data = new float[n];
        fread(data, sizeof(float), n, fin);
        delete data;
    }

    void readVolume(Volume* vol, FILE* fid, int ind) {
        int mind = 100000, maxd = -100000;
        int* data = new int[dimx * dimy * dimz];
        fread(data, sizeof(int), dimx * dimy * dimz, fid);
        int ct = 0;
        for (int x = 0; x < dimx; x++)
            for (int y = 0; y < dimy; y++)
                for (int z = 0; z < dimz; z++) {
                    int d = ((data[ct] >> ind) & 1);
                    vol->setDataAt(x, y, z, (double)d);
                    ct++;

                    if (d < mind) {
                        mind = d;
                    }
                    if (d > maxd) {
                        maxd = d;
                    }
                }
        printf("Min: %d, max: %d\n", mind, maxd);
        delete data;
    }

    int getNumStructures() { return this->numstructs; }

    /* Read volume */
    Volume* getVolume() { return NULL; }

    Volume* getVolume(int ind) {
        Volume* vol = new Volume(dimx, dimy, dimz);
        vol->setSpacing(spx, spy, spz);

        FILE* fid;
        if ((fid = fopen(imgfile, "rb")) == NULL) {
            fprintf(stderr, "Can't open:<%s>\n", imgfile);
            exit(0);
        }

        // By pass Header
        skipFloats(fid, 4 + dimx + dimy + dimz);

        readVolume(vol, fid, ind + 1);

        fclose(fid);

        return vol;
    }

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = spx;
        ay = spy;
        az = spz;
    }

    void getCoordinates(float*& xc, float*& yc, float*& zc) {
        xc = xoff;
        yc = yoff;
        zc = zoff;
    }

private:
    int dimen;
    char mrcfile[1024];
};

class DoseReader : public VolumeReader {
private:
    char imgfile[1024];
    int dimx, dimy, dimz;
    float spx, spy, spz;
    float* xoff;
    float* yoff;
    float* zoff;

public:
    /* Initializer */
    DoseReader(const char* fname) {
        strcpy(imgfile, fname);

        FILE* fid;
        if ((fid = fopen(fname, "rb")) == NULL) {
            fprintf(stderr, "Can't open:<%s>\n", fname);
            exit(0);
        }

        // Read dose Header
        fread(&dimx, sizeof(int), 1, fid);
        fread(&dimy, sizeof(int), 1, fid);
        fread(&dimz, sizeof(int), 1, fid);

        // Read dose coords
        xoff = new float[dimx];
        fread(xoff, sizeof(float), dimx, fid);
        spx = xoff[dimx - 1] - xoff[dimx - 2];

        yoff = new float[dimy];
        fread(yoff, sizeof(float), dimy, fid);
        spy = yoff[dimy - 1] - yoff[dimy - 2];

        zoff = new float[dimz];
        fread(zoff, sizeof(float), dimz, fid);
        spz = zoff[dimz - 1] - zoff[dimz - 2];

        fclose(fid);

        printf("Dimensions: %d %d %d, ImgFile: %s Spacing: %f %f %f\n", dimx,
               dimy, dimz, imgfile, spx, spy, spz);
        printf("Origin: %f %f %f\n", xoff[0], yoff[0], zoff[0]);
    }

    void skipInts(FILE* fin, int n) {
        int* data = new int[n];
        fread(data, sizeof(int), n, fin);
        delete data;
    }

    void skipFloats(FILE* fin, int n) {
        float* data = new float[n];
        fread(data, sizeof(float), n, fin);
        delete data;
    }

    void readVolume(Volume* vol, FILE* fid) {
        float maxd = -100000, mind = 100000;
        float* data = new float[dimx * dimy * dimz];
        fread(data, sizeof(float), dimx * dimy * dimz, fid);
        int ct = 0;
        for (int x = 0; x < dimx; x++)
            for (int y = 0; y < dimy; y++)
                for (int z = 0; z < dimz; z++) {
                    float d = data[ct];
                    vol->setDataAt(x, y, z, d);
                    ct++;

                    if (d < mind) {
                        mind = d;
                    }
                    if (d > maxd) {
                        maxd = d;
                    }
                }
        delete data;
    }

    /* Read volume */
    Volume* getVolume() {
        Volume* vol = new Volume(dimx, dimy, dimz);
        vol->setSpacing(spx, spy, spz);

        FILE* fid;
        if ((fid = fopen(imgfile, "rb")) == NULL) {
            fprintf(stderr, "Can't open:<%s>\n", imgfile);
            exit(0);
        }

        // By pass Header
        skipFloats(fid, 3 + dimx + dimy + dimz);

        readVolume(vol, fid);

        fclose(fid);

        return vol;
    }

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = spx;
        ay = spy;
        az = spz;
    }

    void getCoordinates(float*& xc, float*& yc, float*& zc) {
        xc = xoff;
        yc = yoff;
        zc = zoff;
    }

private:
    int dimen;
    char mrcfile[1024];
};

class CylinderVolumeGenerator : public VolumeReader {
public:
    /* Initializer */
    CylinderVolumeGenerator(int dim, double rad) {
        size = dim;
        radius = rad;
    };

    Volume* getVolume() {
        Volume* vol = new Volume(size, size, size);
        int thick = size / 6;

        double cent = (size - 1) / 2.0;
        // double r2 = radius * radius ;
        for (int x = 0; x < size; x++)
            for (int y = 0; y < size; y++)
                for (int z = 0; z < size; z++) {
                    double dis =
                        sqrt((z - cent) * (z - cent) + (x - cent) * (x - cent));

                    if (fabs(y - cent) > thick) {
                        vol->setDataAt(x, y, z, -(fabs(y - cent) - thick));
                    } else {
                        vol->setDataAt(x, y, z, radius - dis);
                    }

                    // double dis = ( z - cent ) * ( z - cent ) + ( x - cent ) *
                    // ( x - cent ); vol->setDataAt( x, y, z, 10 - 10 * dis / r2
                    // ) ;
                };

        return vol;
    };

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = 1;
        ay = 1;
        az = 1;
    }

private:
    int size;
    double radius;
};

class SheetVolumeGenerator : public VolumeReader {
public:
    /* Initializer */
    SheetVolumeGenerator(int dim, double hit) {
        size = dim;
        height = hit;
    };

    Volume* getVolume() {
        Volume* vol = new Volume(size, size, size);

        double cent = (size - 1) / 2.0;
        for (int x = 0; x < size; x++)
            for (int y = 0; y < size; y++)
                for (int z = 0; z < size; z++) {
                    double dis = fabs((double)z - cent);
                    vol->setDataAt(x, y, z, 10 - 10 * dis / height);

                    // double dis = ( z - cent ) * ( z - cent ) + ( x - cent ) *
                    // ( x - cent ); vol->setDataAt( x, y, z, 10 - 10 * dis / r2
                    // ) ;
                };

        return vol;
    };

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = 1;
        ay = 1;
        az = 1;
    }

private:
    int size;
    double height;
};

class SphereVolumeGenerator : public VolumeReader {
public:
    /* Initializer */
    SphereVolumeGenerator(int dim, double ratiox, double ratioy,
                          double ratioz) {
        size = dim;
        rx = (int)ratiox;
        ry = (int)ratioy;
        rz = (int)ratioz;
    };

    Volume* getVolume() {
        Volume* vol = new Volume(size, size, size);

        double dis = (size - 1) * (size - 1) / 4;

        for (int x = 0; x < size; x++)
            for (int y = 0; y < size; y++)
                for (int z = 0; z < size; z++) {
                    double d =
                        rx * (x - (size - 1) / 2) * (x - (size - 1) / 2) +
                        ry * (y - (size - 1) / 2) * (y - (size - 1) / 2) +
                        rz * (z - (size - 1) / 2) * (z - (size - 1) / 2);
                    vol->setDataAt(x, y, z, 10 - 10 * sqrt(d) / sqrt(dis));
                    // vol->setDataAt( x, y, z, 10 - 10 * (x + y) / size ) ;
                };

        return vol;
    };

    /* Get resolution */
    void getSpacing(float& ax, float& ay, float& az) {
        ax = 1;
        ay = 1;
        az = 1;
    }

private:
    int size, rx, ry, rz;
};

#endif
