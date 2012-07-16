#ifndef EDGECONTAINER_HPP
#define EDGECONTAINER_HPP

#include "maingraph.hpp"
#include "hashmap.h"
#include <vector>
using std::vector;

typedef int colored_etype;
typedef colored_etype bagID;
typedef long long bag_index;
typedef vector<edge> bag;
static const colored_etype ILLEGAL_BAG_ID = -1;
static const bag_index ILLEGAL_BAG_INDEX = -1;
static const int MAX_NUM_BAGS = 16384; //bagIDs use only lowest 14 bits, and 2^14=16384

typedef unsigned long long ECposCode;

inline ECposCode IDtoPosCode(const bagID bid, const bag_index idx) {
	return ((0ULL | bid) | (idx << 14));
}
inline bagID EPCtoBagID(const ECposCode epc) {
	return (epc & 0x00003FFFULL);
}
inline bag_index EPCtoBagIndex(const ECposCode epc) {
	return (epc >> 14);
}

inline bagID getBagID(const bagID color_uv, const bagID color_vu,
                      const bagID color_u, const bagID color_v) {
    if (color_u != color_v) {                       
        if (color_u > color_v) {
           return color_u << 10 | color_v << 6 | color_uv << 3 | color_vu;                
        } else { //if (color_v > color_u) {
           return color_v << 10 | color_u << 6 | color_vu << 3 | color_uv;
        }
    }
    if (color_uv > color_vu) {
       return color_u << 10 | color_v << 6 | color_uv << 3 | color_vu;
    } else {
       return color_v << 10 | color_u << 6 | color_vu << 3 | color_uv;           
    }
    return 0;                           
}

struct EdgeContainer
{
	bag bags[MAX_NUM_BAGS]; 
	vector<bagID> full_bags;
	int currentFBIndex;
	bagID currentBag;
	bag_index currentIndex;
	hash_map<edge, ECposCode> edgePositions;
	inline void put(const edge e, const bagID et);  
	inline void loggedPut(const edge e, const bagID et);	
	//inline void remove(bagID id, bag_index bi);  
	bag_index getBagSize(bagID et) {return bags[et].size();}  

	EdgeContainer() : currentFBIndex(ILLEGAL_BAG_ID), 
		              currentIndex(ILLEGAL_BAG_INDEX)
		              {};
	inline void resetIteration() { currentFBIndex = ILLEGAL_BAG_ID;
                                   currentBag = ILLEGAL_BAG_ID;  
	                               currentIndex = ILLEGAL_BAG_INDEX; }
	inline bool hasNextBag(); 
	inline void goNextBag()   { ++currentFBIndex; 
	                            currentBag = full_bags[currentFBIndex]; 
								currentIndex = ILLEGAL_BAG_INDEX;}
	inline bagID getCurrentBag() { return currentBag; }  
	inline bag_index getCurrentIndex() { return currentIndex; }  
	
	inline bool hasNextElement()   { return currentIndex + 1 != bags[currentBag].size(); }
	inline void goNextElement()    { ++currentIndex; }
	inline edge getCurrentElement() { return (bags[currentBag])[currentIndex]; }   
	inline edge getIndexedElement(const bag_index idx) { return (bags[currentBag])[idx]; } 
	inline edge getElement(const bag_index idx, const bagID bid) { return (bags[bid])[idx]; }
	inline void setIndexedElement(const bag_index idx, const edge e) { (bags[currentBag])[idx] = e; }
	inline void loggedSetIndexedElement(const bag_index idx, const edge e) { 
                                                 edgePositions.erase((bags[currentBag])[idx]);
                                                 (bags[currentBag])[idx] = e; 
                                                 edgePositions[e] = IDtoPosCode(currentBag,idx); }
	inline void loggedReplaceElement(const edge old, const edge neu) { 
                                                 const ECposCode epc = edgePositions[old];
                                                 edgePositions.erase(old);
                                                 const bag_index idx = EPCtoBagIndex(epc);
                                                 const bagID bid = EPCtoBagID(epc); 
                                                 (bags[bid])[idx] = neu; 
                                                 edgePositions[neu] = epc; }
	inline void loggedSwapElements(const edge old, const edge neu) { 
                                                 const ECposCode epc_o = edgePositions[old],
                                                                 epc_n = edgePositions[neu];
                                                 const bag_index idx_o = EPCtoBagIndex(epc_o),
                                                                 idx_n = EPCtoBagIndex(epc_n);
                                                 const bagID bid_o = EPCtoBagID(epc_o),
                                                             bid_n = EPCtoBagID(epc_n);                                                 
                                                 (bags[bid_o])[idx_o] = neu; 
                                                 (bags[bid_n])[idx_n] = old; 
                                                 edgePositions[old] = epc_n; 
                                                 edgePositions[neu] = epc_o; }
                                                                                                  
    inline void eraseEdgePosition(const edge e) { edgePositions.erase(e); } 
};

inline bool EdgeContainer::hasNextBag()  { 
	bag_index bipp =  currentFBIndex + 1;
	// Check for Dirtyness)
	while (bipp != full_bags.size()) {
		bag_index bi = full_bags[bipp];
		if (bags[bi].size() > 0)
			return true;
		full_bags[bipp] = full_bags.back();
		full_bags.pop_back();
	}
	return false;  
} 

/*inline void EdgeContainer::remove(bagID id, bag_index bi) {
	(bags[id])[bi] = bags[id].back();
	bags[id].pop_back();
	//Adjust iterator if necessary
	if (bi == currentBag)
	  --currentIndex;
}*/

inline void EdgeContainer::put(const edge e, const bagID et) {
	if (bags[et].size() == 0) { //First element to be put in there
		full_bags.push_back(et);
		bags[et].push_back(e);
	} else {
		bags[et].push_back(e);	
	}
}

inline void EdgeContainer::loggedPut(const edge e, const bagID et) {
    const bag::size_type bagsize = bags[et].size();
	if (bagsize == 0) { //First element to be put in there
		full_bags.push_back(et);
		bags[et].push_back(e);
		edgePositions[e] = IDtoPosCode(et,bagsize);
	} else {
		bags[et].push_back(e);	
		edgePositions[e] = IDtoPosCode(et,bagsize);
	}
}

#endif
