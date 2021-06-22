#pragma once

#include <iosfwd>
#include <vector>

#include "DOFVectorBase.hpp"

namespace Shtns
{
  class FileWriter
  {
  public:
    virtual ~FileWriter() {}

    virtual void write(Mesh const& mesh,
                       std::vector<DOFVector const*> const& values = {}) const = 0;
  };

} // end namespace Shtns
