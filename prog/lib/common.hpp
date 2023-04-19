#ifndef COMMON_HPP
#define COMMON_HPP

#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <array>
#include <functional>
#include <string>
#include <cmath>

#if !defined(__PRETTY_FUNCTION__) && !defined(__GNUC__)
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

// ================ not implemented macro
#define _THROW_NOT_IMP_ \
	{ \
		printf("NOT IMPLEMENTED ERROR:\n"); \
		printf("function: %s\nat:       %s:%i\n", __PRETTY_FUNCTION__, __FILE__, __LINE__); \
		throw std::runtime_error("not implemented"); \
	}

#define _THROW_UNREACHABLE_ \
	{ \
		printf("UNREACHABLE CODE ERROR:\n"); \
		printf("function: %s\nat:       %s:%i\n", __PRETTY_FUNCTION__, __FILE__, __LINE__); \
		throw std::runtime_error("unreachable code error"); \
	}

#endif

