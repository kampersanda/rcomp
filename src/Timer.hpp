/**
 * @file Timer.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include <chrono>

namespace rcomp {

/**
 * A class for timer.
 */
class Timer {
  private:
    std::chrono::high_resolution_clock::time_point m_start;
    std::chrono::high_resolution_clock::time_point m_stop;

  public:
    Timer() {
        start();
    }

    void start() {
        m_start = std::chrono::high_resolution_clock::now();
    }
    void stop() {
        m_stop = std::chrono::high_resolution_clock::now();
    }

    double stop_and_get_sec() {
        m_stop = std::chrono::high_resolution_clock::now();
        return get_elapsed_sec();
    }
    double stop_and_get_ms() {
        m_stop = std::chrono::high_resolution_clock::now();
        return get_elapsed_ms();
    }
    double stop_and_get_us() {
        m_stop = std::chrono::high_resolution_clock::now();
        return get_elapsed_us();
    }
    double stop_and_get_ns() {
        m_stop = std::chrono::high_resolution_clock::now();
        return get_elapsed_ns();
    }

    double get_elapsed_sec() const {
        auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(m_stop - m_start);
        return dur.count() / 1000.0;
    }
    double get_elapsed_ms() const {
        auto dur = std::chrono::duration_cast<std::chrono::microseconds>(m_stop - m_start);
        return dur.count() / 1000.0;
    }
    double get_elapsed_us() const {
        auto dur = std::chrono::duration_cast<std::chrono::nanoseconds>(m_stop - m_start);
        return dur.count() / 1000.0;
    }
    double get_elapsed_ns() const {
        auto dur = std::chrono::duration_cast<std::chrono::nanoseconds>(m_stop - m_start);
        return static_cast<double>(dur.count());
    }
};

}  // namespace rcomp
