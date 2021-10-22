/**
 * @file GroupedLData_UnitTraits.hpp
 */
#pragma once

#include "basics.hpp"

namespace rcomp {

template <bool t_WithSAE, class t_Fnode>
struct GroupedLData_UnitTraits;

template <class t_Fnode>
struct GroupedLData_UnitTraits<false, t_Fnode> {
    struct unit_type {
        uchar_type chr;  // of run
        size_type exp;  // of run
        t_Fnode* fnode;
        size_type frnid;  // run-ID in the F-node
        size_type lrnid;  // run-ID in the L-node (of this)
    };
};

template <class t_Fnode>
struct GroupedLData_UnitTraits<true, t_Fnode> {
    struct unit_type {
        uchar_type chr;  // of run
        size_type exp;  // of run
        t_Fnode* fnode;
        size_type frnid;  // run-ID in the F-node
        size_type lrnid;  // run-ID in the L-node (of this)
        size_type sae;  // SA-entry
    };
};

}  // namespace rcomp
