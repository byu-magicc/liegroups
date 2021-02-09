#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_UTILITIES_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_UTILITIES_

#include "lie_groups/state.h"

namespace lie_groups { namespace utilities
{


template< template<typename , int > class tG>
struct GroupIsSEN {
    static constexpr bool value = false;
};

template<> 
inline struct GroupIsSEN<lie_groups::SE2> {
    static constexpr bool value = true;
};

template<> 
inline struct GroupIsSEN<lie_groups::SE3> {
    static constexpr bool value = true;
};

//---------------------------------------------------------------------------

template< template<typename , int > class tG>
struct GroupIsSON {
    static constexpr bool value = false;
};

template<> 
inline struct GroupIsSON<lie_groups::SO2> {
    static constexpr bool value = true;
};

template<> 
inline struct GroupIsSON<lie_groups::SO3> {
    static constexpr bool value = true;
};


//---------------------------------------------------------------------------


template< template<typename , int > class tG>
struct GroupIsRN {
    static constexpr bool value = false;
};

template<> 
inline struct GroupIsRN<lie_groups::Rn> {
    static constexpr bool value = true;
};


//---------------------------------------------------------------------------


template <typename State>
struct StateIsSEN_seN{
    static constexpr bool value = GroupIsSEN<State::template GroupTemplate>::value;
};

//---------------------------------------------------------------------------

template <typename State>
struct StateIsSON_soN{
    static constexpr bool value = GroupIsSON<State::template GroupTemplate>::value;
};

//---------------------------------------------------------------------------

template <typename State>
struct StateIsRN_rN{
    static constexpr bool value = GroupIsRN<State::template GroupTemplate>::value;
};

} // namespace utilities    
} // namespace lie_groups
#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_UTILITIES_