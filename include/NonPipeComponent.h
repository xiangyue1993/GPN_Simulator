#pragma once

#include <iostream>
#include <cstdlib>
#include <vector>

#include "SelfDefinedVariables.h"

using namespace std;

class NonPipeComponent
{
public:
	int nonpipe_comp_id;
	string comp_name;
	NonPipeComponentType comp_type;
};