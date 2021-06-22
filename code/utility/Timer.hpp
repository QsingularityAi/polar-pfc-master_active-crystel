#pragma once

#include <chrono>

namespace Shtns
{
  /// time measurement methods
  class Timer
  {
    using Clock     = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    using fsec      = std::chrono::duration<double>;

  public:
    /// initializes the timer with current time
    Timer();

    /// resets the timer to current time
    void reset();

    /// returns the elapsed time (from construction or last reset) to now in seconds
    double elapsed() const;

  private:
    /// start time
    TimePoint t0;
  };

}  // end namespace Shtns
