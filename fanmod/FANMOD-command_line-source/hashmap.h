#ifndef HASHMAP64_HPP
#define HASHMAP64_HPP

#include<ext/hash_map>
//#include <hash_map>


using __gnu_cxx::hash_map;
namespace __gnu_cxx {
    template <> struct hash <unsigned long long > {
	size_t operator  () (unsigned long long __x) const {
	    return __x;
	}
    };
    template <> struct hash <const unsigned long long > {
	size_t operator  () (const unsigned long long __x) const {
	    return __x;
	}
    };
}

#endif
