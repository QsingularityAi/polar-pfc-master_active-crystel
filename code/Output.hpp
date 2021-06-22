#pragma once

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

/**
 * \def SHTNS_NO_THROW
 * \brief The preprocessor constant sets whether to use c-asserts (if defined) or
 * to throw an exception in case of an error (if not defined).
 **/
#ifdef SHTNS_NO_THROW
  #include <cassert>
#else
  #include <stdexcept>
#endif

/**
 * \def SHTNS_ENABLE_MSG_DBG
 * \brief The preprocessor constant enables the functions \ref Shtns::msg_dbg
 * and \ref Shtns::assert_msg_dbg.
 *
 * If the value is set to 1 the functions \ref Shtns::msg_dbg and \ref Shtns::assert_msg_dbg
 * are implemented, otherwise empty. Default is value 0 if \ref NDEBUG is not
 * defined, otherwise value 1.
 **/
#ifndef SHTNS_ENABLE_MSG_DBG
  #ifndef NDEBUG
    #define SHTNS_ENABLE_MSG_DBG 1
  #else
    #define SHTNS_ENABLE_MSG_DBG 0
  #endif
#endif


namespace Shtns
{
  namespace aux
  {
    template <class OStream>
    OStream& concat(OStream& out) { return out; }

    template <class OStream, class Arg0, class... Args>
    OStream& concat(OStream& out, Arg0&& arg0, Args&&... args)
    {
      out << arg0; concat(out, std::forward<Args>(args)...);
      return out;
    }

    template <class... Args>
    std::string to_string(Args&&... args)
    {
      std::stringstream ss; concat(ss, std::forward<Args>(args)...);
      return ss.str();
    }
  }

  /// \brief print a message
  /**
   * Example:
   * ```
   * msg("Hello ", "World: ", 123); // prints "Hello World: 123\n"
   * ```
   **/
  template <class... Args>
  void msg(Args&&... args)
  {
    // std::cout << std::setprecision(6) << std::scientific;
    aux::concat(std::cout, std::forward<Args>(args)..., "\n");
  }


  /// \brief print a message (without appended newline)
  /**
   * Example:
   * ```
   * msg("Hello ", "World: ", 123); // prints "Hello World: 123"
   * ```
   **/
  template <class... Args>
  void msg_(Args&&... args)
  {
    // std::cout << std::setprecision(6) << std::scientific;
    aux::concat(std::cout, std::forward<Args>(args)...);
  }


  /// \brief print a message and exit
  /**
   * If the preprocessor constant \ref SHTNS_NO_THROW is defined,
   * the c-assert macro is called, otherwise an exception of
   * type \ref std::runtime_Error is thrown.
   **/
  template <class... Args>
  void error_exit(Args&&... args)
  {
#ifdef SHTNS_NO_THROW
    aux::concat(std::cerr, "ERROR: ", std::forward<Args>(args)..., "\n");
    #ifndef NDEBUG
      assert(false);
    #else
      std::exit(EXIT_FAILURE);
    #endif
#else
    throw std::runtime_error( aux::to_string("ERROR: ", std::forward<Args>(args)...) );
#endif
  }


  /// \brief test for condition and in case of failure print message and exit
  /**
   * This function is equivalent to
   * ```
   * if (condition == false) error_exit(text);
   * ```
   * where `text` correspond to the arguments passed after the
   * \p condition argument.
   **/
  template <class... Args>
  void assert_msg(bool condition, Args&&... args)
  {
    if (!condition) { error_exit(std::forward<Args>(args)...); }
  }


  /// \brief test for condition and in case of failure print message
  /**
   * Same as \ref assert_msg but does not throw an exception, or call assert.
   * It just tests for the condition and prints a message with prepended
   * string "WARNING".
   **/
  template <class... Args>
  void warn_msg(bool condition, Args&&... args)
  {
    if (!condition) { msg("WARNING: ", std::forward<Args>(args)...); }
  }


#if SHTNS_ENABLE_MSG_DBG
  /// \brief print message, in debug mode only
  /**
   * Same as \ref msg, but is available only if preprocessor constant
   * \ref SHTNS_ENABLE_MSG_DBG is set to 1, otherwise the function is empty.
   **/
  template <class... Args>
  void msg_dbg(Args&&... args) { msg(std::forward<Args>(args)...); }


  /// \brief call assert_msg, in debug mode only
  /**
   * Same as \ref assert_msg, but is available only if preprocessor constant
   * \ref SHTNS_ENABLE_MSG_DBG is set to 1, otherwise the function is empty.
   **/
  template <class... Args>
  void assert_msg_dbg(bool condition, Args&&... args)
  {
    assert_msg(condition, std::forward<Args>(args)...);
  }
#else
  template <class... Args>
  void msg_dbg(Args&&...) {}

  template <class... Args>
  void assert_msg_dbg(bool, Args&&...) {}
#endif

} // end namespace Shtns
