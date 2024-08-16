#ifndef __NEIGHBORMAP_H__
#define __NEIGHBORMAP_H__

#include <map>
#include <string> // for std::string

#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h"
#include "hoomd/Compute.h"
#include "hoomd/HOOMDMPI.h"

#include <iostream> // for std::cout
#include <stdexcept> // for std::out_of_range
#include <algorithm> // for std::min and std::max
#include <vector>    // for std::vector


class NeighborManager {
private:
    std::map<unsigned int, std::vector<unsigned int>> neighborMap; // map of neighbors 
    //std::map<std::string, int> dataMap;  // Using std::map to store key-value pairs
    //std::map<std::string, int> dataMap = {{"key1", 100}, {"key2", 200}};  // Using std::map to store key-value pairs with dummy values


public:
    // Deleted copy constructor and assignment operator to prevent copying
    NeighborManager(const NeighborManager&) = delete;
    NeighborManager& operator=(const NeighborManager&) = delete;

    // Static method to get the single instance of NeighborManager
    static NeighborManager& getInstance() {
        static NeighborManager instance;
        return instance;
    }
    
    NeighborManager() {}  // Default constructor

    // Method to get data corresponding to a key
    bool getNeighbors(unsigned int tag_i, std::vector<unsigned int>& value) const {
        auto it = neighborMap.find(tag_i);
        if (it != neighborMap.end()) {
            value = it->second;  // Return value associated with the key
            return true;
        } else {
            return false;
        }
    }

    // Method to set or update data for a key
    void setNeighbors(unsigned int tag_i, const std::vector<unsigned int>& value) {
        neighborMap[tag_i] = value;  // Insert or update the value associated with the key
    }

    // Method to add a single unsigned int to the vector for a key
    void addNeighbor(unsigned int tag_i, unsigned int tag_j) {
        neighborMap[tag_i].push_back(tag_j);  // Add tag_j to the vector associated with tag_i
    }


    // Method to clear all data from the map
    void clearNeighbors() {
        neighborMap.clear();
    }

    // Method to print all data for debugging purposes
    void printNeighbors() const {
        for (const auto& pair : neighborMap) {
            std::cout << "Tag " << pair.first << ": ";
            for (const auto& tag : pair.second) {
                std::cout << tag << " ";
            }
            std::cout << std::endl;
        }
    }


    // Method to count the number of keys in the map
    //int countKeys() const {
    //    return neighborMap.size();
    //}


    // Method to returnt he vector of neighbors for a specific tag_i
    std::vector<unsigned int> getNeighborsForTag(unsigned int tag_i) const {
        std::vector<unsigned int> result;
        auto it = neighborMap.find(tag_i);
        if (it != neighborMap.end()) {
            result = it->second;
        }
        return result;
    }

    // Method to retrieve all keys currently in the BondMap
    std::vector<unsigned int> getAllKeys() const {
        std::vector<unsigned int> keys;
        for (const auto& pair : neighborMap) {
            keys.push_back(pair.first);
        }
        return keys;
    }


    // Method to remove data from the map for a specific key
    bool removeKey(unsigned int tag_i) {
        auto it = neighborMap.find(tag_i);
        if (it != neighborMap.end()) {
            neighborMap.erase(it);
            return true;
        } else {
            return false;
        }
    }


};

extern NeighborManager neighbor_manager;  // Declaration of a global instance of NeighborManager

#endif // __NEIGHBORMAP_H__

