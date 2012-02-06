/*
 * RBDL - Rigid Body Library
 * Copyright (c) 2011 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "rbdlconfig.h"

#include <iostream>
namespace RigidBodyDynamics {

void rbdl_print_version() {
	std::cout << "RigidBodyDynamicsLibrary version:" << std::endl
		<< "  revision     : " << RBDL_BUILD_REVISION
		<< " (branch: " << RBDL_BUILD_BRANCH << ")" << std::endl
		<< "  build type   : " << RBDL_BUILD_TYPE << std::endl
#ifdef RBDL_ENABLE_LOGGING
		<< "  logging      : on (warning: reduces performance!)" << std::endl
#else
		<< "  logging      : off" << std::endl
#endif
#ifdef RBDL_USE_SIMPLE_MATH
		<< "  simplemath   : on (warning: reduces performance!)" << std::endl
#else
		<< "  simplemath   : off" << std::endl
#endif
		;
}

}
