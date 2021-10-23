/**
 * @file Rlbwt_Interface.hpp
 */
#pragma once

#include "basics.hpp"

namespace rcomp {

//! An interface of RLBWT compressor.
class Rlbwt_Interface {
  public:
    //! Default destructor
    virtual ~Rlbwt_Interface() = default;

    //! Check if the data structure is empty.
    virtual bool is_empty() const = 0;

    //! Get the number of stored characters.
    virtual size_type get_num_chars() const = 0;

    //! Get the allocated memory in bytes.
    virtual size_type get_memory_in_bytes(bool include_this = true) const = 0;

    //! Print the statistics related to memory (via stdout).
    virtual void show_memory_statistics() const = 0;

    //! Print the detailed statistics (via stdout).
    virtual void show_detailed_statistics() const = 0;

    //! Print the statistics measured with internal monitors (via stdout).
    virtual void show_monitored_statistics() const = 0;

    //! Extend the RLBWT text by appending the given character (i.e., T := c + T).
    virtual void extend(const uchar_type new_chr) = 0;

    //! Output the original text (in the input order) via the callback function.
    virtual void decode_text(const std::function<void(uchar_type)>& fn) const = 0;

    //! Output the RLBWT text via the callback function.
    virtual void output_runs(const std::function<void(const run_type&)>& fn) const = 0;

    //! Test the data structure.
    virtual void test() const = 0;
};

}  // namespace rcomp
