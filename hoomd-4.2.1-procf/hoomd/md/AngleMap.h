#ifndef __ANGLEMAP_H__
#define __ANGLEMAP_H__

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


class AngleManager {
private:
    std::map<unsigned int, double> angleMap; // map of a-i-b bond angles
    //std::map<std::string, int> dataMap;  // Using std::map to store key-value pairs
    //std::map<std::string, int> dataMap = {{"key1", 100}, {"key2", 200}};  // Using std::map to store key-value pairs with dummy values


public:
    // Deleted copy constructor and assignment operator to prevent copying
    AngleManager(const AngleManager&) = delete;
    AngleManager& operator=(const AngleManager&) = delete;

    // Static method to get the single instance of AngleManager
    static AngleManager& getInstance() {
        static AngleManager instance;
        return instance;
    }
    
    AngleManager() {}  // Default constructor

    // Method to get data corresponding to a key
    bool getAngle(unsigned int angle_index, double value) const {
        auto it = angleMap.find(angle_index);
        if (it != angleMap.end()) {
            value = it->second;  // Return value associated with the key
            return true;
        } else {
            return false;
        }
    }

    // Method to set or update data for a key
    void setAngle(unsigned int angle_index, double value) {
        angleMap[angle_index] = value;  // Insert or update the value associated with the key
    }

    // Method to clear all data from the map
    void clearAngles() {
        angleMap.clear();
    }

    // Method to print all data for debugging purposes
    //void printAngles() const {
    //    for (const auto& pair : angleMap) {
    //        std::cout << "Angle " << pair.first << ": ";
    //        for (const auto& tag : pair.second) {
    //            std::cout << tag << " ";
    //        }
    //        std::cout << std::endl;
    //    }
    //}


    // Method to retrieve all keys currently in the BondMap
    std::vector<unsigned int> getAllKeys() const {
        std::vector<unsigned int> keys;
        for (const auto& pair : angleMap) {
            keys.push_back(pair.first);
        }
        return keys;
    }


    // Method to remove data from the map for a specific key
    bool removeKey(unsigned int angle_index) {
        auto it = angleMap.find(angle_index);
        if (it != angleMap.end()) {
            angleMap.erase(it);
            return true;
        } else {
            return false;
        }
    }


};

extern AngleManager angle_manager;  // Declaration of a global instance of AngleManager

#endif // __ANGLEMAP_H__

