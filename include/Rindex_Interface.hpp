/**
 * @file Rindex_Interface.hpp
 */
#pragma once

#include "Rlbwt_Interface.hpp"

namespace rcomp {

//! An interface of RLBWT compressor.
class Rindex_Interface : public Rlbwt_Interface {
  public:
    //! Default destructor
    virtual ~Rindex_Interface() = default;

    //! Count the number of occurrences for the query pattern.
    virtual size_type count(range_type<const uchar_type*> pat) const = 0;

    //! Report the positions of occurrence for the query pattern via the callback function.
    virtual void locate(range_type<const uchar_type*> pat, const std::function<void(size_type)>& fn) const = 0;

    //! Extract SA-entries from the head of the BWT-text.
    virtual void extract_sa_entries(const std::function<void(size_type)>& fn) const = 0;
};

}  // namespace rcomp
