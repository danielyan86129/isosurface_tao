#ifndef POLYTOPE_H
#define POLYTOPE_H

#include "Cycle.h"
#include "geocommon.h"
#include <stdio.h>

class Cube {
public:
    /// the cycle table
    LinkedList<Cycle*>** cycleTable;

public:
    /**
     * Constructor.
     */
    Cube(const char* filename);

    /**
     * Destructor
     */
    ~Cube(void);

    /**
     * Returns cycles for a particular sign configuration
     */
    LinkedList<Cycle*>* getCycle(unsigned char mask);
};

#endif