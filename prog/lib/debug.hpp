#ifndef BF_DEBUG_HPP
#define BF_DEBUG_HPP

#include <vector>
#include "grid.hpp"
#include "slae/csrmat.hpp"

namespace dbg{

void print(const bf::CsrStencil& mat, const std::vector<double>& val);
void print(const bf::CsrMatrix& mat);

}

#endif
