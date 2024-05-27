#ifndef __DG_COMMON_HPP__
#define __DG_COMMON_HPP__

#include <iostream>
#include <vector>
#include <memory>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include "catch2/catch.hpp"

#define _THROW_NOT_IMP_ \
{\
	std::ostringstream oss; \
	oss << "NOT IMPLEMENTED method called." << std::endl; \
	oss << __PRETTY_FUNCTION__ << std::endl; \
	oss << "At" << std::endl; \
	oss << __FILE__ << ": " << __LINE__ << std::endl; \
	throw std::runtime_error(oss.str()); \
}

constexpr size_t INVALID_INDEX = (size_t)-1;
#if !defined(__PRETTY_FUNCTION__) && !defined(__GNUC__)
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif


#endif
