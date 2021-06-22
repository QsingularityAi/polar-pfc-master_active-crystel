#pragma once

#include <type_traits>

namespace Shtns
{

  template <class E, class Integer,
    class = typename std::enable_if<std::is_enum<E>::value>::type >
  constexpr bool is_a(E a, Integer b)
  {
    return (a & b) != 0;
  }

  /// The basic element data structur
  struct Element
  {
    enum Types : std::size_t {
      NO_TYPE   = 0,
      VERTEX    = 1<<0,
      LINE2     = 1<<1,
      LINE3     = 1<<2,
      LINE4     = 1<<3,
      LINE5     = 1<<4,
      TRIANGLE3 = 1<<5,
      TRIANGLE6 = 1<<6,
      TRIANGLE9 = 1<<7,
      TRIANGLE10= 1<<8,
      TRIANGLE12= 1<<9,
      TRIANGLE15= 1<<10,
      QUAD4     = 1<<11,
      QUAD8     = 1<<12,
      QUAD9     = 1<<13,

      // combined categories
      LINE     = LINE2 | LINE3 | LINE4 | LINE5,
      TRIANGLE = TRIANGLE3 | TRIANGLE6 | TRIANGLE9 | TRIANGLE10 | TRIANGLE12 | TRIANGLE15,
      QUAD     = QUAD4 | QUAD8 | QUAD9,

      ELEMENT_0D  = VERTEX,
      ELEMENT_1D  = LINE,
      ELEMENT_2D  = TRIANGLE | QUAD
    };

    static constexpr std::size_t nodes(Types T)
    {
      return T == VERTEX     ? 1 :
             T == LINE2      ? 2 :
             T == LINE3      ? 3 :
             T == TRIANGLE3  ? 3 :
             T == TRIANGLE6  ? 6 :
             T == QUAD4      ? 4 :
             T == QUAD8      ? 8 :
             T == QUAD9      ? 9 : 0;
    }

    static constexpr std::size_t corners(Types T)
    {
      return is_a(T, VERTEX)   ? 1 :
             is_a(T, LINE)     ? 2 :
             is_a(T, TRIANGLE) ? 3 :
             is_a(T, QUAD)     ? 4 : 0;
    }
  };

} // end namespace Shtns;
