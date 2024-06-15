#include "bflow/fem_grid.hpp"

using namespace bflow;

FemGrid::FemGrid(const GraphGrid& grid) : _power(grid.n_midnodes)
{
    if (_power > 6)
    {
        throw std::runtime_error("Only power=1,2,3,4,5,6 is allowed");
    }
    if (grid.n_edges() < 1)
    {
        throw std::runtime_error("Only single vessel grids are allowed");
    }
    std::vector<double> len_by_edge;
    _h.resize(grid.n_elem());
    double len = 0.0;
    int elem = 0;
    // green - points from grid
    // red   - points from this
    // blue  - nodes
    std::set<int> used_greens;
    std::vector<int> red_to_green;
    for (size_t i = 0; i < grid.n_edges(); i++)
    {
        auto a = grid.points_by_edge(i);
        if (grid.tab_point_nodes(a[0]).size() < 3)
        {
            _points.push_back(len);
            if (used_greens.find(a[0]) == used_greens.end())
            {
                red_to_green.push_back(a[0]);
                used_greens.insert(a[0]);
            }
            int left = _nodes.size();
            _nodes.push_back(len);
            for (int j = 1; j < a.size() - 1; j++)
            {
                len += grid.find_cell_length(elem);
                _h[elem] = len;
                elem++;
                _points.push_back(len);
                if (used_greens.find(a[j]) == used_greens.end())
                {
                    red_to_green.push_back(a[j]);
                    used_greens.insert(a[j]);
                }
                _nodes.push_back(len);
                _nodes.push_back(len);
            }
            len += grid.find_cell_length(elem);
            _h[elem] = len;
            elem++;
            _points.push_back(len);
            if (used_greens.find(a.back()) == used_greens.end())
            {
                red_to_green.push_back(a.back());
                used_greens.insert(a.back());
            }
            //_nodes_by_point.push_back(x[j]);
            _nodes.push_back(len);
            int right = _nodes.size() - 1;
            _nodes_by_edge.push_back({left, right});
        }
        else
        {
            int left = _nodes.size();
            double l = len_by_edge[0];
            _nodes.push_back(l);
            if (used_greens.find(a[0]) == used_greens.end())
            {
                red_to_green.push_back(a[0]);
                used_greens.insert(a[0]);
            }
            for (int j = 1; j < a.size() - 1; j++)
            {
                l += grid.find_cell_length(elem);
                _h[elem] = len;
                elem++;
                _points.push_back(l);
                if (used_greens.find(a[j]) == used_greens.end())
                {
                    red_to_green.push_back(a[j]);
                    used_greens.insert(a[j]);
                }
                _nodes.push_back(l);
                _nodes.push_back(l);
            }
            l += grid.find_cell_length(elem);
            _h[elem] = len;
            elem++;
            _points.push_back(l);
            if (used_greens.find(a.back()) == used_greens.end())
            {
                red_to_green.push_back(a.back());
                used_greens.insert(a.back());
            }
            _nodes.push_back(l);
            int right = _nodes.size() - 1;
            _nodes_by_edge.push_back({left, right});
        }
        len_by_edge.push_back(len);
    }

    for (int i = 0; i < grid.n_elem(); i++)
    {
        _h[i] = grid.find_cell_length(i);
        _len += grid.find_cell_length(i);
    }
    for (int i = 0; i < red_to_green.size(); i++)
    {
        _nodes_by_point.push_back(grid.tab_point_nodes(red_to_green[i]));
    }

    double l = _h[0] / _power;
    for (size_t i = 0; i < _points.size() - 1; i++)
    {
        double p = _points[i] + l;
        for (size_t j = 1; j < _power; j++)
        {
            _nodes.push_back(p);
            p += l;
        }
    }
    fill_f_vec();
}

double FemGrid::h(int ielem) const
{
    return _h[ielem];
}

size_t FemGrid::n_nodes() const
{
    return _nodes.size();
}

size_t FemGrid::n_elements() const
{
    return _nodes.size() / 2;
}

size_t FemGrid::n_points() const
{
    return _points.size();
}

size_t FemGrid::n_local_bases() const
{
    return _power + 1;
}

double FemGrid::node(size_t i) const
{
    return _nodes[i];
}

std::vector<double> FemGrid::load_vector() const
{
    return _f_vec;
}

double FemGrid::full_length() const
{
    return _len;
}

CsrMatrix FemGrid::mass_matrix() const
{
    CsrMatrix ret(stencil());
    std::vector<double> local = local_mass_matrix();

    for (size_t ielem = 0; ielem < n_elements(); ++ielem)
    {
        std::vector<int> lg = tab_elem_nodes(ielem);

        for (size_t irow = 0; irow < n_local_bases(); ++irow)
            for (size_t icol = 0; icol < n_local_bases(); ++icol)
            {
                double v = h(ielem) / 2 * local[irow * n_local_bases() + icol];
                size_t iaddr = ret.find_index(lg[irow], lg[icol]);
                ret.vals()[iaddr] += v;
            }
    }

    return ret;
}

CsrMatrix FemGrid::block_transport_matrix() const
{
    CsrMatrix ret(stencil());
    std::vector<double> local = local_transport_matrix();

    for (int ielem = 0; ielem < n_elements(); ++ielem)
    {
        std::vector<int> lg = tab_elem_nodes(ielem);

        // block diagonal
        for (int irow = 0; irow < n_local_bases(); ++irow)
            for (int icol = 0; icol < n_local_bases(); ++icol)
            {
                double v = local[irow * n_local_bases() + icol];
                size_t iaddr = ret.find_index(lg[irow], lg[icol]);
                ret.vals()[iaddr] += v;
            }
    }

    return ret;
}

CsrMatrix FemGrid::block_u_transport_matrix(const std::vector<double>& u) const
{
    size_t dim = n_local_bases();
    CsrMatrix ret(stencil());
    std::vector<double> local_u = local_u_transport_matrix();

    for (int ielem = 0; ielem < n_elements(); ++ielem)
    {
        std::vector<int> lg = tab_elem_nodes(ielem);

        std::vector<double> local(dim * dim, 0);
        for (size_t i = 0; i < dim; ++i)
        {
            for (size_t j = 0; j < dim; ++j)
            {
                for (size_t k = 0; k < dim; ++k)
                {
                    local[dim * i + j] += local_u[dim * dim * i + dim * j + k] * u[lg[k]];
                }
            }
        }

        // block diagonal
        for (int irow = 0; irow < n_local_bases(); ++irow)
            for (int icol = 0; icol < n_local_bases(); ++icol)
            {
                double v = local[irow * n_local_bases() + icol];
                size_t iaddr = ret.find_index(lg[irow], lg[icol]);
                ret.vals()[iaddr] += v;
            }
    }

    return ret;
}

CsrMatrix FemGrid::coupling_transport_matrix() const
{
    CsrMatrix ret(stencil());
    for (size_t ielem = 0; ielem < n_elements(); ++ielem)
    {
        std::vector<int> lg = tab_elem_nodes(ielem);
        // upwind coupling
        if (ielem > 0)
        {
            // left
            size_t iaddr = ret.find_index(lg[0], lg[0] - 1);
            ret.vals()[iaddr] -= 1;
        }
        {
            // right
            size_t iaddr = ret.find_index(lg[1], lg[1]);
            ret.vals()[iaddr] += 1;
        }
    }

    return ret;
}

std::vector<int> FemGrid::tab_elem_nodes(int ielem) const
{
    std::vector<int> ret;
    // clang-format off
    if (_power == 1)
    {
        ret = {2 * ielem, 2 * ielem + 1};
    }
    else if (_power == 2)
    {
        ret = {2 * ielem, 2 * ielem + 1, 2 * int(n_elements()) + ielem};
    }
    else if (_power == 3)
    {
        ret = {2 * ielem, 2 * ielem + 1,
               2 * int(n_elements()) + 2 * ielem,
               2 * int(n_elements()) + 2 * ielem + 1};
    }
    else if (_power == 4)
    {
        ret = {2 * ielem, 2 * ielem + 1,
               2 * int(n_elements()) + 3 * ielem,
               2 * int(n_elements()) + 3 * ielem + 1,
               2 * int(n_elements()) + 3 * ielem + 2};
    }
    else if (_power == 5)
    {
        ret = {2 * ielem,
               2 * ielem + 1,
               2 * int(n_elements()) + 4 * ielem,
               2 * int(n_elements()) + 4 * ielem + 1,
               2 * int(n_elements()) + 4 * ielem + 2,
               2 * int(n_elements()) + 4 * ielem + 3};
    }
    else if (_power == 6)
    {
        ret = {2 * ielem,
               2 * ielem + 1,
               2 * int(n_elements()) + 5 * ielem,
               2 * int(n_elements()) + 5 * ielem + 1,
               2 * int(n_elements()) + 5 * ielem + 2,
               2 * int(n_elements()) + 5 * ielem + 3,
               2 * int(n_elements()) + 5 * ielem + 4};
    }
    // clang-format on
    else
    {
        _THROW_NOT_IMP_;
    }
    return ret;
}

CsrMatrix FemGrid::stencil() const
{
    if (_stencil.n_rows() < 1)
    {
        std::vector<std::set<int>> lod(n_nodes());
        // block diagonal
        for (size_t ielem = 0; ielem < n_elements(); ++ielem)
        {
            std::vector<int> lg = tab_elem_nodes(ielem);

            for (size_t irow = 0; irow < n_local_bases(); ++irow)
                for (size_t icol = 0; icol < n_local_bases(); ++icol)
                {
                    lod[lg[irow]].insert(lg[icol]);
                }
        }
        // coupling
        for (int ipoint = 0; ipoint < n_points() - 1; ++ipoint)
        {
            std::vector<int> nodes = tab_point_nodes(ipoint);
            if (nodes.size() == 2)
            {
                lod[nodes[0]].insert(nodes[1]);
                lod[nodes[1]].insert(nodes[0]);
            }
            else if (nodes.size() == 3)
            {
                lod[nodes[0]].insert(nodes[1]);
                lod[nodes[0]].insert(nodes[2]);
                lod[nodes[1]].insert(nodes[0]);
                lod[nodes[1]].insert(nodes[2]);
                lod[nodes[2]].insert(nodes[0]);
                lod[nodes[2]].insert(nodes[1]);
            }
        }
        _stencil.set_stencil(lod);
    }
    return _stencil;
}

std::vector<int> FemGrid::tab_point_nodes(int ipoint) const
{
    return _nodes_by_point[ipoint];
}

std::vector<double> FemGrid::local_mass_matrix() const
{
    // clang-format off
    if (_power == 1)
    {
        return {2.0 / 3.0, 1.0 / 3.0,
                1.0 / 3.0, 2.0 / 3.0};
    }
    else if (_power == 2)
    {
        return {4.0 / 15.0, -1.0 / 15.0,  2.0 / 15.0,
               -1.0 / 15.0,  4.0 / 15.0,  2.0 / 15.0,
                2.0 / 15.0,  2.0 / 15.0, 16.0 / 15.0};
    }
    else if (_power == 3)
    {
        return {16.0 / 105.0,  19.0 / 840.0,  33.0 / 280.0, -3.0 / 70,
                19.0 / 840.0,  16.0 / 105.0,  -3.0 / 70.0,  33.0 / 280.0,
                33.0 / 280.0 , -3.0 / 70.0,   27.0 / 35.0, -27.0 / 280.0,
                -3.0 / 70.0,   33.0 / 280.0, -27.0 / 280.0, 27.0 / 35.0};
    }
    else if (_power == 4)
    {
        return {292.0/2835, -29.0/2835, 296.0/2835,  -58.0/945,    8.0/405,
                -29.0/2835, 292.0/2835,   8.0/405,   -58.0/945,  296.0/2835,
                296.0/2835,   8.0/405,  256.0/405,  -128.0/945,  256.0/2835,
                -58.0/945,  -58.0/945, -128.0/945,   208.0/315, -128.0/945,
                  8.0/405,  296.0/2835, 256.0/2835, -128.0/945,  256.0/405};
    }
    else if (_power == 5)
    {
        return { 1907.0/24948,    493.0/88704,   24775.0/266112,  -9925.0/133056,  -1525.0/133056,  17125.0/399168,
                  493.0/88704,   1907.0/24948,   -1525.0/133056,  17125.0/399168,  24775.0/266112,  -9925.0/133056,
                24775.0/266112, -1525.0/133056, 111625.0/199584, -24625.0/133056, -62875.0/798336,   2125.0/14784,
                -9925.0/133056, 17125.0/399168, -24625.0/133056,  62375.0/99792,    2125.0/14784,  -13625.0/66528,
                -1525.0/133056, 24775.0/266112, -62875.0/798336,   2125.0/14784,  111625.0/199584, -24625.0/133056,
                17125.0/399168, -9925.0/133056,   2125.0/14784,  -13625.0/66528,  -24625.0/133056,  62375.0/99792};
    }
    else if (_power == 6)
    {
        return 
        {90269.0/1501500, -10237.0/3003000, 42087.0/500500, -16971.0/200200, 10237.0/150150,   -687.0/20020,   3867.0/500500,
        -10237.0/3003000,  90269.0/1501500,  3867.0/500500,   -687.0/20020,  10237.0/150150, -16971.0/200200, 42087.0/500500,
         42087.0/500500,    3867.0/500500,  64692.0/125125,  -4887.0/20020,   5688.0/25025,  -14607.0/100100,  8532.0/125125,
        -16971.0/200200,    -687.0/20020,   -4887.0/20020,    2619.0/4004,   -3231.0/10010,    9693.0/40040, -14607.0/100100,
         10237.0/150150,   10237.0/150150,   5688.0/25025,   -3231.0/10010,  10544.0/15015,   -3231.0/10010,   5688.0/25025,
          -687.0/20020,   -16971.0/200200, -14607.0/100100,   9693.0/40040,  -3231.0/10010,    2619.0/4004,   -4887.0/20020,
          3867.0/500500,   42087.0/500500,   8532.0/125125, -14607.0/100100,  5688.0/25025,   -4887.0/20020,  64692.0/125125};
    }
    // clang-format on
    else
    {
        _THROW_NOT_IMP_;
    }
}

std::vector<double> FemGrid::local_transport_matrix() const
{
    // clang-format off
    if (_power == 1)
    {
        return {-0.5, -0.5, 0.5, 0.5};
    }
    else if (_power == 2)
    {
        return {-1.0 / 2.0, 1.0 / 6.0, -2.0 / 3.0,
                -1.0 / 6.0, 1.0 / 2.0, 2.0 / 3.0, 
                2.0 / 3.0, -2.0 / 3.0, 0.0};
    }
    else if (_power == 3)
    {
        return {-1.0 / 2.0,   -7.0 / 80.0, -57.0 / 80.0, 3.0 / 10.0,
                 7.0 / 80.0,   1.0 / 2.0,   -3.0 / 10.0, 57.0 / 80.0,
                57.0 / 80.0,   3.0 / 10.0,     0.0,     -81.0 / 80.0,
                -3.0 / 10.0, -57.0 / 80.0,  81.0 / 80.0,    0.0};
    }
    else if (_power == 4)
    {
        return { -1.0/2,    107.0/1890, -736.0/945, 134.0/315, -64.0/315,
               -107.0/1890,   1.0/2,      64.0/315,-134.0/315, 736.0/945,
               736.0/945,  -64.0/315,      0.0,   -352.0/315, 512.0/945,
               -134.0/315,  134.0/315,   352.0/315,    0.0,   -352.0/315,
                 64.0/315, -736.0/945,  -512.0/945, 352.0/315,    0.0 };        
    }
    else if (_power == 5)
    {
        return {-1.0/2,       -5951.0/145152, -123425.0/145152,  6325.0/10368,  1675.0/10368,   -575.0/1512,
              5951.0/145152,      1.0/2,        -1675.0/10368,    575.0/1512, 123425.0/145152, -6325.0/10368,
            123425.0/145152,   1675.0/10368,         0.0,       -6875.0/5184, -19375.0/48384,  51875.0/72576,
             -6325.0/10368,    -575.0/1512,      6875.0/5184,        0.0,      51875.0/72576, -38125.0/36288,
             -1675.0/10368, -123425.0/145152,   19375.0/48384, -51875.0/72576,      0.0,        6875.0/5184,
               575.0/1512,     6325.0/10368,   -51875.0/72576,  38125.0/36288, -6875.0/5184,        0.0};
    }
    else if (_power == 6)
    {
        return {
           -1.0/2,      587.0/18480, -1776.0/1925, 5151.0/6160, -3967.0/5775,   732.0/1925, -267.0/1925,
         -587.0/18480,    1.0/2,       267.0/1925, -732.0/1925,  3967.0/5775, -5151.0/6160, 1776.0/1925,
         1776.0/1925,  -267.0/1925,       0.0,    -3078.0/1925,  2136.0/1925,  -243.0/385,   648.0/1925,
        -5151.0/6160,   732.0/1925,   3078.0/1925,     0.0,       -87.0/77,    3807.0/6160, -243.0/385, 
         3967.0/5775, -3967.0/5775,  -2136.0/1925,   87.0/77,        0.0,       -87.0/77,   2136.0/1925,
         -732.0/1925,  5151.0/6160,    243.0/385, -3807.0/6160,    87.0/77,        0.0,    -3078.0/1925,
          267.0/1925, -1776.0/1925,   -648.0/1925,  243.0/385,  -2136.0/1925,  3078.0/1925,     0.0 };
    }
    // clang-format on
    else
    {
        _THROW_NOT_IMP_;
    }
}

std::vector<double> FemGrid::local_u_transport_matrix() const
{
    // clang-format off
    if (_power == 1)
    {
        return {-1.0/3, -1.0/6, -1.0/6, -1.0/3, 
                 1.0/3,  1.0/6,  1.0/6,  1.0/3};
    }
    else if (_power == 2){
        return {-1.0/3, 1.0/30, -1.0/5, 1.0/30, 1.0/15, 1.0/15, -1.0/5, 1.0/15, -8.0/15,
                -1.0/15, -1.0/30, -1.0/15, -1.0/30, 1.0/3, 1.0/5.0, -1.0/15, 1.0/5, 8.0/15,                
                2.0/5, 0, 4.0/15, 0, -2.0/5, -4.0/15, 4.0/15, -4.0/15, 0};
               
    }
    // clang-format on
    else
    {
        _THROW_NOT_IMP_;
    }
}

void FemGrid::fill_f_vec()
{
    if (_power == 1)
    {
        for (size_t i = 0; i < _nodes.size(); i++)
            _f_vec.push_back(_h[0] / 2);
    }
    if (_power == 2)
    {
        for (size_t i = 0; i < _points.size(); i++)
        {
            _f_vec.push_back(_h[0] / 3 / 2);
            _f_vec.push_back(_h[0] / 3 / 2);
        }
        for (int i = 0; i < _points.size(); i++)
            _f_vec.push_back(4.0 * _h[0] / 3 / 2);
    }
    if (_power == 3)
    {
        for (size_t i = 0; i < _points.size(); i++)
        {
            _f_vec.push_back(_h[0] / 4 / 2);
            _f_vec.push_back(_h[0] / 4 / 2);
        }
        for (int i = 0; i < _points.size(); i++)
        {
            _f_vec.push_back(3.0 * _h[0] / 4 / 2);
            _f_vec.push_back(3.0 * _h[0] / 4 / 2);
        }
    }
    if (_power == 4)
    {
        for (size_t i = 0; i < _points.size(); i++)
        {
            _f_vec.push_back(7.0 * _h[0] / 45 / 2);
            _f_vec.push_back(7.0 * _h[0] / 45 / 2);
        }
        for (int i = 0; i < _points.size(); i++)
        {
            _f_vec.push_back(32.0 * _h[0] / 45 / 2);
            _f_vec.push_back(4.0 * _h[0] / 15 / 2);
            _f_vec.push_back(32.0 * _h[0] / 45 / 2);
        }
    }
    if (_power == 5)
    {
        for (size_t i = 0; i < _points.size(); i++)
        {
            _f_vec.push_back(19.0 * _h[0] / 144 / 2);
            _f_vec.push_back(19.0 * _h[0] / 144 / 2);
        }
        for (int i = 0; i < _points.size(); i++)
        {
            _f_vec.push_back(25.0 * _h[0] / 48 / 2);
            _f_vec.push_back(25.0 * _h[0] / 72 / 2);
            _f_vec.push_back(25.0 * _h[0] / 48 / 2);
            _f_vec.push_back(25.0 * _h[0] / 72 / 2);
        }
    }
    if (_power == 6)
    {
        for (size_t i = 0; i < _points.size(); i++)
        {
            _f_vec.push_back(41.0 * _h[0] / 420 / 2);
            _f_vec.push_back(41.0 * _h[0] / 420 / 2);
        }
        for (int i = 0; i < _points.size(); i++)
        {
            _f_vec.push_back(18.0 * _h[0] / 35 / 2);
            _f_vec.push_back(9.0 * _h[0] / 140 / 2);
            _f_vec.push_back(68.0 * _h[0] / 105 / 2);
            _f_vec.push_back(9.0 * _h[0] / 140 / 2);
            _f_vec.push_back(18.0 * _h[0] / 35 / 2);
        }
    }
}

size_t FemGrid::closest_node(double x) const
{
    double t = 1e16;
    size_t ret = 0;
    for (size_t i = 0; i < n_nodes(); ++i)
    {
        double t1 = std::abs(_nodes[i] - x);
        if (t1 < t)
        {
            t = t1;
            ret = i;
        }
    }
    return ret;
};

size_t FemGrid::closest_node(size_t section, double p) const
{
    auto nbe = _nodes_by_edge[section];
    double min = std::abs(_nodes[nbe[0]] - p);
    size_t res = nbe[0];
    for (int i = nbe[0] + 1; i <= nbe[1]; i++)
    {
        double x = std::abs(_nodes[i] - p);
        if (std::abs(_nodes[i] - p) < min)
        {
            min = std::abs(_nodes[i] - p);
            res = i;
        }
    }
    if (p > _nodes[res] && res != nbe[0])
        res++;

    return res;
}
size_t FemGrid::tab_node_elem(size_t inode) const
{
    if (inode < (n_points() - 1) * 2)
    {
        return inode / 2;
    }
    else
    {
        return (inode - n_points() * 2 + 2) / _power;
    }
};
