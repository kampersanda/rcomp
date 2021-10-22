#pragma once

#include "FData_Naive.hpp"
#include "LData_Naive.hpp"
#include "LFIntervalGraph.hpp"
#include "Rindex_LFIG.hpp"

#include "GroupedFData_Naive.hpp"
#include "GroupedFData_Serialized.hpp"
#include "GroupedLData_Naive.hpp"
#include "GroupedLData_Serialized.hpp"
#include "GroupedLFIntervalGraph.hpp"
#include "Rindex_GLFIG.hpp"

namespace rcomp::rindex_types {

template <size_type t_DivBound = 7>
struct lfig_naive {
    // Components
    using ldata_type = LData_Naive<true>;
    using fdata_type = FData_Naive<ldata_type>;
    using graph_type = LFIntervalGraph<ldata_type, fdata_type, t_DivBound>;
    // Compressor
    using type = Rindex_LFIG<graph_type>;
};

template <size_type t_GroupedBound, size_type t_DivBound = 7>
struct glfig_naive {
    // Components
    using ldata_type = GroupedLData_Naive<t_GroupedBound, true>;
    using fdata_type = GroupedFData_Naive<ldata_type>;
    using graph_type = GroupedLFIntervalGraph<ldata_type, fdata_type, t_DivBound>;
    // Compressor
    using type = Rindex_GLFIG<graph_type>;
};

template <size_type t_GroupedBound, size_type t_DivBound = 7>
struct glfig_serialized {
    // Components
    using ldata_type = GroupedLData_Serialized<t_GroupedBound, true>;
    using fdata_type = GroupedFData_Serialized<ldata_type>;
    using graph_type = GroupedLFIntervalGraph<ldata_type, fdata_type, t_DivBound>;
    // Compressor
    using type = Rindex_GLFIG<graph_type>;
};

}  // namespace rcomp::rindex_types