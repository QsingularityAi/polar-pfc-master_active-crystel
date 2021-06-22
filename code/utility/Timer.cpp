#include "Timer.hpp"

namespace Shtns {

Timer::Timer()
  : t0(Clock::now())
{}

void Timer::reset()
{
  t0 = Clock::now();
}

double Timer::elapsed() const
{
  auto t1 = Clock::now();
  fsec fs = t1 - t0;
  return fs.count();
}

} // end namespace Shtns
