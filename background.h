
#pragma once

#include "background.h"
#include "main.h"

struct Background {
    //pass cosmological parameters 
    Background() { 
        integrate();
    } 

private:
    void integrate() {}
};

#include <vector>
