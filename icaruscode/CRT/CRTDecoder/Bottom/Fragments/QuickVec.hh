 // This file (QuickVec.hh) was created by Ron Rechenmacher <ron@fnal.gov> on
 // Sep  3, 2014. "TERMS AND CONDITIONS" governing this file are in the README
 // or COPYING file. If you do not have such a file, one can be obtained by
 // contacting Ron or Fermi Lab in Batavia IL, 60510, phone: 630-840-3000.
 // $RCSfile: QuickVec.hh,v $
 // rev="$Revision: 1.8 $$Date: 2014/09/05 19:21:11 $";
#ifndef QuickVec_hh
#define QuickVec_hh

extern "C" {
#include <stdint.h>
}

#include <cassert>
#include <cstddef>		// ptrdiff_t
#include <string.h>		// memcpy
#include <strings.h>		// bzero
#include <utility>		// std::swap 
#include <memory>		// unique_ptr
#include <vector>

//#include "trace.h"		// TRACE
#ifndef TRACE
# define TRACE( lvl, ... )
# define UNDEF_TRACE_AT_END
#endif

#define USE_UNIQUE_PTR 0

#ifndef QUICKVEC_DO_TEMPLATE
# define QUICKVEC_DO_TEMPLATE    1
#endif

#undef NOT_OLD_CXXSTD
#if !defined(__GCCXML__) && defined(__GXX_EXPERIMENTAL_CXX0X__)
# define NOT_OLD_CXXSTD 1
#endif

#if QUICKVEC_DO_TEMPLATE == 0
# ifndef QUICKVEC_TT
#  define QUICKVEC_TT unsigned long long
# endif
# define TT_		 QUICKVEC_TT
# define QUICKVEC_TEMPLATE
# define QUICKVEC          QuickVec
# define QUICKVEC_TN       QuickVec
# define QUICKVEC_VERSION
#else
# define QUICKVEC_TEMPLATE template <typename TT_>
# define QUICKVEC          QuickVec<TT_>
# define QUICKVEC_TN       typename QuickVec<TT_>
# define QUICKVEC_VERSION  static short Class_Version() { return 4; } // proper version for templates
#endif

QUICKVEC_TEMPLATE
struct QuickVec
{
    typedef       TT_*        iterator;
    typedef const TT_*  const_iterator;
    typedef       TT_&       reference;
    typedef const TT_& const_reference;
    typedef       TT_       value_type;
    typedef ptrdiff_t  difference_type;
    typedef    size_t        size_type;

    QuickVec( size_t sz );
    QuickVec( size_t sz, TT_ val );
#  if USE_UNIQUE_PTR == 0
    ~QuickVec();
#   define PTR_(xx) xx
#  else
#   define PTR_(xx) xx.get()
#  endif
    QuickVec( std::vector<TT_> & other )
	: size_(other.size()), data_(new TT_[other.capacity()]), capacity_(other.capacity())
    {   memcpy( PTR_(data_), (void*)&other[0], size_*sizeof(TT_) );
    }
    void clear() { size_=0; }

    QuickVec( const QuickVec & other ) //= delete; // non construction-copyable
	: size_(other.size_), data_(new TT_[other.capacity_]), capacity_(other.capacity_)
    {	TRACE( 10, "QuickVec copy ctor this=%p data_=%p other.data_=%p size_=%d other.size_=%d"
	      , (void*)this, (void*)PTR_(data_), (void*)PTR_(other.data_), size_, other.size_ );
	memcpy( PTR_(data_), PTR_(other.data_), size_*sizeof(TT_) );
    }
    QUICKVEC & operator=( const QuickVec & other ) //= delete; // non copyable
    {	TRACE( 10, "QuickVec copy assign this=%p data_=%p other.data_=%p size_=%d other.size_=%d"
	      , (void*)this, (void*)PTR_(data_), (void*)PTR_(other.data_), size_, other.size_ );
	resize( other.size_ );
	memcpy( PTR_(data_), PTR_(other.data_), size_*sizeof(TT_) );
	return *this;
    }
#  if NOT_OLD_CXXSTD
    QuickVec( QuickVec && other )	 // construction-movable
	: size_(other.size_), data_(std::move(other.data_)), capacity_(other.capacity_)
    {   TRACE( 10, "QuickVec move ctor this=%p data_=%p other.data_=%p"
	      , (void*)this, (void*)PTR_(data_), (void*)PTR_(other.data_) );
#      if USE_UNIQUE_PTR == 0
	other.data_ = nullptr;
#      endif
    }
    QUICKVEC & operator=( QuickVec && other ) // assign movable
    {   TRACE( 10, "QuickVec move assign this=%p data_=%p other.data_=%p"
	      , (void*)this, (void*)PTR_(data_), (void*)PTR_(other.data_) );
	size_ = other.size_;
	delete [] data_;
	data_ = std::move(other.data_);
	capacity_ = other.capacity_;
#      if USE_UNIQUE_PTR == 0
	other.data_ = nullptr;
#      endif
	return *this;
    }
#  endif

    TT_& operator[](int idx);
    const TT_& operator[](int idx) const;
    size_t size() const;
    size_t capacity() const;
    iterator       begin();
    const_iterator begin() const;
    iterator       end();
    const_iterator end() const;
    void reserve( size_t size );
    void resize( size_t size );
    void resize( size_t size, TT_ val );
    iterator insert( const_iterator position, size_t nn, const TT_& val );
    iterator insert(  const_iterator position, const_iterator first
		    , const_iterator last);
    iterator erase( const_iterator first, const_iterator last );
    void     swap( QuickVec& x );
    void push_back( const value_type& val );

    QUICKVEC_VERSION

private:
    // Root needs the size_ member first. It must be of type int.
    // Root then needs the [size_] comment after data_.
    // Note: NO SPACE between "//" and "[size_]"
    unsigned size_;
#  if USE_UNIQUE_PTR == 0
    TT_ * data_; //[size_]
#  else
    std::unique_ptr<TT_[]> data_;
#  endif
    unsigned capacity_;
};

QUICKVEC_TEMPLATE
inline QUICKVEC::QuickVec( size_t sz )
    : size_(sz), data_(new TT_[sz]), capacity_(sz)
{   TRACE( 15, "QuickVec %p ctor sz=%d data_=%p", (void*)this, size_, (void*)PTR_(data_) );
}
QUICKVEC_TEMPLATE
inline QUICKVEC::QuickVec( size_t sz, TT_ val )
    : size_(sz), data_(new TT_[sz]), capacity_(sz)
{   TRACE( 15, "QuickVec %p ctor sz=%d/v data_=%p", (void*)this, size_, (void*)PTR_(data_) );
    for (iterator ii=begin(); ii!=end(); ++ii) *ii=val;
    //bzero( &data_[0], (sz<4)?(sz*sizeof(TT_)):(4*sizeof(TT_)) );
}

#if USE_UNIQUE_PTR == 0
QUICKVEC_TEMPLATE
inline QUICKVEC::~QuickVec()
{   TRACE( 15, "QuickVec %p dtor start data_=%p size_=%d"
	  , (void*)this, (void*)PTR_(data_), size_ );
    delete [] data_;
    TRACE( 15, "QuickVec %p dtor return", (void*)this );
}
#endif

QUICKVEC_TEMPLATE
inline TT_& QUICKVEC::operator[](int idx)
{   assert(idx<(int)size_); return data_[idx];
}

QUICKVEC_TEMPLATE
inline const TT_& QUICKVEC::operator[](int idx) const
{   assert(idx < (int)size_);
    return data_[idx];
}

QUICKVEC_TEMPLATE
inline size_t QUICKVEC::size()     const { return size_; }
QUICKVEC_TEMPLATE
inline size_t QUICKVEC::capacity() const { return capacity_; }

QUICKVEC_TEMPLATE
inline QUICKVEC_TN::iterator       QUICKVEC::begin()       { return iterator(PTR_(data_)); }
QUICKVEC_TEMPLATE
inline QUICKVEC_TN::const_iterator QUICKVEC::begin() const { return iterator(PTR_(data_)); }
QUICKVEC_TEMPLATE
inline QUICKVEC_TN::iterator       QUICKVEC::end()       { return iterator(PTR_(data_)+size_); }
QUICKVEC_TEMPLATE
inline QUICKVEC_TN::const_iterator QUICKVEC::end() const { return const_iterator(PTR_(data_)+size_); }

QUICKVEC_TEMPLATE
inline void QUICKVEC::reserve( size_t size )
{   if (size > capacity_)  // reallocation if true
    {
#      if USE_UNIQUE_PTR == 0
	TT_ * old=data_;
	data_ = new TT_[size];
	memcpy( data_, old, size_*sizeof(TT_) );
	TRACE( 13, "QUICKVEC::reserve this=%p old=%p data_=%p"
	      , (void*)this, (void*)old, (void*)data_ );
	delete [] old;
#      else
	std::unique_ptr<TT_[]> old=std::move(data_);
	data_ = std::unique_ptr<TT_[]>(new TT_[size]);
	memcpy( data_.get(), old.get(), size_*sizeof(TT_) );
#      endif
	capacity_ = size;
	// bye to old(unique_ptr)
    }
}

QUICKVEC_TEMPLATE
inline void QUICKVEC::resize( size_t size )
{   if      (size <  size_)     size_ = size; // decrease
    else if (size <= capacity_) size_ = size;
    else // increase/reallocate 
    {
#      if USE_UNIQUE_PTR == 0
	TT_ * old=data_;
	data_ = new TT_[size];
	memcpy( data_, old, size_*sizeof(TT_) );
	TRACE( 13, "QUICKVEC::resize this=%p old=%p data_=%p"
	      , (void*)this, (void*)old, (void*)data_ );
	delete [] old;
#      else
	std::unique_ptr<TT_[]> old=std::move(data_);
	data_ = std::unique_ptr<TT_[]>(new TT_[size]);
	memcpy( data_.get(), old.get(), size_*sizeof(TT_) );
#      endif
	size_ = capacity_ = size;
	// bye to old(unique_ptr)
    }
}
QUICKVEC_TEMPLATE
inline void QUICKVEC::resize( size_type size, TT_ val )
{   size_type old_size=size;
    resize( size );
    if (size > old_size)
	for (iterator ii=begin()+old_size; ii!=end(); ++ii) *ii=val;
}

QUICKVEC_TEMPLATE
inline QUICKVEC_TN::iterator QUICKVEC::insert(  const_iterator position
							, size_t nn
							, const TT_& val)
{   assert(position<=end());  // the current end
    size_t offset=position-begin();
    reserve( size_+nn );      // may reallocate and invalidate "position"

    iterator dst=end()+nn;    // for shifting existing data after
    iterator src=end();       // insertion point
    size_t cnt=end()-(begin()+offset);
    while (cnt--) *--dst = *--src;

    dst=begin()+offset;
    size_+=nn;
    while (nn--) *dst++ = val;
    return begin()+offset;
}
QUICKVEC_TEMPLATE
inline QUICKVEC_TN::iterator QUICKVEC::insert(  const_iterator position
				       , const_iterator first
				       , const_iterator last)
{   assert(position<=end());  // the current end
    size_t nn=(last-first);
    size_t offset=position-begin();
    reserve( size_+nn );      // may reallocate and invalidate "position"

    iterator dst=end()+nn;    // for shifting existing data after
    iterator src=end();       // insertion point
    size_t cnt=end()-(begin()+offset);
    while (cnt--) *--dst = *--src;

    dst=begin()+offset;
    size_+=nn;
    while (nn--) *dst++ = *first++;
    return begin()+offset;
}

QUICKVEC_TEMPLATE
inline QUICKVEC_TN::iterator QUICKVEC::erase( const_iterator first
				      ,const_iterator last )
{   assert(last<=end());  // the current end
    size_t nn=(last-first);
    size_t offset=first-begin();

    iterator dst=begin()+offset;    // for shifting existing data from last
    iterator src=dst+nn;            // to first
    size_t cnt=end()-src;
    while (cnt--) *dst++ = *src++;

    size_-=nn;
    return begin()+offset;
}

QUICKVEC_TEMPLATE
inline void QUICKVEC::swap( QuickVec& x )
{   TRACE( 12, "QUICKVEC::swap this=%p enter data_=%p x.data_=%p"
	  , (void*)this, (void*)PTR_(data_), (void*)PTR_(x.data_) );
    std::swap( data_, x.data_ );
    std::swap( size_, x.size_ );
    std::swap( capacity_, x.capacity_ );
    TRACE( 12, "QUICKVEC::swap return data_=%p x.data_=%p"
	  , (void*)PTR_(data_), (void*)PTR_(x.data_) );
}

QUICKVEC_TEMPLATE
inline void QUICKVEC::push_back( const value_type& val )
{   if (size_ == capacity_)
    {   reserve( size_ + size_/10 + 1 );
    }
    *end() = val;
    ++size_;
}

#ifdef UNDEF_TRACE_AT_END
# undef TRACE
#endif
#endif /* QuickVec_hh */
