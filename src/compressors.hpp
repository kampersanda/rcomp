/**
 * @file compressors.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include "Compressor_LFIG.hpp"
#include "FData_Naive.hpp"
#include "FIndex.hpp"
#include "LData_Naive.hpp"
#include "LIndex.hpp"

#include "Compressor_GroupLFIG.hpp"
#include "GroupFData_Naive.hpp"
#include "GroupFData_Serialized.hpp"
#include "GroupFIndex.hpp"
#include "GroupLData_Naive.hpp"
#include "GroupLData_Serialized.hpp"
#include "GroupLIndex.hpp"

namespace rcomp::compressors {

template <size_type t_DivBound = 7>
struct lfig_naive {
    using ldata_type      = LData_Naive;
    using fdata_type      = FData_Naive<ldata_type>;
    using lindex_type     = LIndex<ldata_type>;
    using findex_type     = FIndex<fdata_type>;
    using compressor_type = Compressor_LFIG<lindex_type, findex_type, t_DivBound>;
};

template <size_type t_GroupBound, size_type t_DivBound = 7, bool t_SumExpHint = true, bool t_LookupHint = true>
struct glfig_naive {
    using ldata_type      = GroupLData_Naive<t_GroupBound, t_SumExpHint, t_LookupHint>;
    using fdata_type      = GroupFData_Naive<ldata_type>;
    using lindex_type     = GroupLIndex<ldata_type>;
    using findex_type     = GroupFIndex<fdata_type>;
    using compressor_type = Compressor_GroupLFIG<lindex_type, findex_type, t_DivBound>;
};

template <size_type t_GroupBound, size_type t_DivBound = 7, bool t_SumExpHint = true, bool t_LookupHint = true>
struct glfig_serialized {
    using ldata_type      = GroupLData_Serialized<t_GroupBound, t_SumExpHint, t_LookupHint>;
    using fdata_type      = GroupFData_Serialized<ldata_type>;
    using lindex_type     = GroupLIndex<ldata_type>;
    using findex_type     = GroupFIndex<fdata_type>;
    using compressor_type = Compressor_GroupLFIG<lindex_type, findex_type, t_DivBound>;
};

}  // namespace rcomp::compressors
