/** \file flags.hpp */

#pragma once

  // NOTE: imported from project AMDiS
  /// The Flag class encapsulates flags which represents simple information.
  class Flag
  {	
  public:
    /// Constructs a unset Flag
    inline Flag() : flags(0) {}

    /// Constructs a Flag initialized by f
    inline Flag(const unsigned long f) : flags(f) {}

    /// Copy constructor
    inline Flag(const Flag& f) : flags(f.flags) {}

    /// Destructor
    inline ~Flag() {}

    /// Compares two Flags
    inline bool operator==(const Flag& f) const 
    {
      return (flags == f.flags);
    }

    /// Compares two Flags
    bool operator!=(const Flag& f) const 
    {
      return !(f == *this);
    }

    /// Assignment operator
    inline Flag& operator=(const Flag& f) 
    {
      if (this != &f) 
	flags = f.flags; 
      return *this;
    }

    /// Typecast
    inline operator bool() const 
    { 
      return isAnySet(); 
    }

    /// Set \ref flags 
    inline void setFlags(const unsigned long f) 
    { 
      flags = f; 
    }

    /// Set \ref flags
    inline void setFlags(const Flag& f) 
    { 
      flags = f.flags; 
    }

    /// Sets \ref flags to \ref flags | f
    inline void setFlag(const unsigned long f) 
    { 
      flags |= f; 
    }

    /// Sets \ref flags to \ref flags | f.flags
    inline void setFlag(const Flag& f) 
    { 
      flags |= f.flags; 
    }

    /// Sets \ref flags to \ref flags & ~f
    inline void unsetFlag(const unsigned long f) 
    { 
      flags &= ~f; 
    }

    /// Sets \ref flags to \ref flags & ~f.flags
    inline void unsetFlag(const Flag& f) 
    { 
      flags &= ~f.flags; 
    }

    inline unsigned long getFlags() const 
    { 
      return flags; 
    }

    /// Returns \ref flags | f.flags
    inline Flag operator+(const Flag& f) const 
    {
      Flag r(flags); 
      r.setFlag(f); 
      return r;
    }

    /// Returns \ref flags & ~f.flags
    inline Flag operator-(const Flag& f) const 
    {
      Flag r(flags); 
      r.unsetFlag(f); 
      return r;
    }

    /// Returns \ref flags | f.flags
    inline Flag operator|(const Flag& f) const 
    {
      Flag r(flags); 
      r.setFlag(f); 
      return r;
    }

    /// Returns \ref flags & f.flags
    inline Flag operator&(const Flag& f) const 
    {
      Flag r(flags);
      r.flags &= f.flags; 
      return r;
    }

    /// Sets \ref flags to \ref flags &= f.flags
    inline Flag operator&=(const Flag& f) 
    {
      flags &= f.flags;
      return *this;
    }

    /// Returns \ref flags ^ f.flags
    inline Flag operator^(const Flag& f) const 
    {
      Flag r(flags);
      r.flags ^= f.flags;
      return r;
    }

    /// Sets \ref flags to \ref flags & f.flags
    inline Flag& operator|=(const Flag& f) 
    {
      if (this != &f)
	flags |= f.flags;
      return *this;
    }

    /// Returns ~\ref flags
    inline Flag operator~() const 
    { 
      Flag r;
      r.flags = ~flags; 
      return r;
    }

    /// Checks whether all set bits of f.flags are set in \ref flags too.
    inline bool isSet(const Flag& f) const 
    {
      return ((flags&f.flags) == f.flags);
    }

    /// Returns !\ref isSet(f)
    inline bool isUnset(const Flag& f) const 
    { 
      return !isSet(f); 
    }

    /// Returns true if \ref flags != 0
    inline bool isAnySet() const 
    { 
      return (flags != 0); 
    }

  protected:	
    /// Internal flag representation
    unsigned long  flags;
  };


