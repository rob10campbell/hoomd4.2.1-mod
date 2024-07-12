#ifndef __BONDMAP_H__
#define __BONDMAP_H__

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


struct BondTrackData {
    uint64_t formedTime;
    uint64_t brokeTime = 0;
    double pos1_x;
    double pos1_y;
    double pos1_z;
    double pos2_x;
    double pos2_y;
    double pos2_z;
    double d1;
    double d2;
    int type1;
    int type2;

    //BondTrackData() : formedTime(0), pos1_x(0.0f), pos1_y(0.0f), pos1_z(0.0f), pos2_x(0.0f), pos2_y(0.0f), pos2_z(0.0f), d1(0), d2(0), type1(0), type2(0) {}

};

class BondDataManager {
private:
    std::map<std::pair<unsigned int, unsigned int>, BondTrackData> dataMap;  // Using std::map to store key-value pairs
    //std::map<unsigned int, int> dataMap;  // Using std::map to store key-value pairs
    //std::map<std::string, int> dataMap;  // Using std::map to store key-value pairs
    //std::map<std::string, int> dataMap = {{"key1", 100}, {"key2", 200}};  // Using std::map to store key-value pairs with dummy values
    //int importantData;  // Example of the data you want to abstract

    // Helper function to create a canonical key (ensures tags are stored in order for easier parsing)
    std::pair<unsigned int, unsigned int> makeCanonicalKey(unsigned int tag_i, unsigned int tag_j) const {
        return std::make_pair(std::min(tag_i, tag_j), std::max(tag_i, tag_j));
    }



public:
    // Deleted copy constructor and assignment operator to prevent copying
    BondDataManager(const BondDataManager&) = delete;
    BondDataManager& operator=(const BondDataManager&) = delete;

    // Static method to get the single instance of BondDataManager
    static BondDataManager& getInstance() {
        static BondDataManager instance;
        return instance;
    }

    //BondDataManager(int data = 0) : importantData(data) {}  // Constructor with default value

    //int getData() const { return importantData; }   // Inline method to get the data
    //void setData(int newData) { importantData = newData; }  // Inline method to set the data
    
    BondDataManager() {}  // Default constructor

    // Method to get data corresponding to a key
    //bool getData(const std::string& key) const {
    //bool getData(const std::string& key, int& value) const {
    //bool getData(unsigned int key, int& value) const {
    bool getData(unsigned int tag_i, unsigned int tag_j, BondTrackData& value) const {
        //auto it = dataMap.find(key);
        //auto it = dataMap.find(std::make_pair(tag_i, tag_j));
        auto it = dataMap.find(makeCanonicalKey(tag_i, tag_j));
        if (it != dataMap.end()) {
            value = it->second;  // Return value associated with the key
            return true;
        } else {
            // Handle error case, such as key not found
            // For example, you can throw an exception or return a default value
            //throw std::out_of_range("Key not found in data map");
            return false;
        }
    }


    // Method to set or update data for a key
    //void setData(const std::string& key, int value) {
    //void setData(unsigned int key, int value) {
    void setData(unsigned int tag_i, unsigned int tag_j, const BondTrackData& value) {
        //dataMap[key] = value;  // Insert or update the value associated with the key
        //dataMap[std::make_pair(tag_i, tag_j)] = value;  // Insert or update the value associated with the key
        dataMap[makeCanonicalKey(tag_i, tag_j)] = value;  // Insert or update the value associated with the key
    }


    // Method to clear all data from the map
    void clearData() {
        dataMap.clear();
        //std::cout << "All data cleared from DataManager." << std::endl;
    }


    // Method to print all data for debugging purposes
    void printData() const {
        for (const auto& pair : dataMap) {
            //std::cout << pair.first << ": " << pair.second << std::endl;
            std::cout << "(" << pair.first.first << ", " << pair.first.second << "): formedTime=" << pair.second.formedTime 
                      << ", brokeTime(" << pair.second.brokeTime 
                      << ", pos1=(" << pair.second.pos1_x 
                      << "," << pair.second.pos1_y
                      << "," << pair.second.pos1_z 
                      << "), pos2=(" << pair.second.pos2_x 
                      << "," << pair.second.pos2_y
                      << "," << pair.second.pos2_z 
                      << "), d1=" << pair.second.d1
                      << ", d2=" << pair.second.d2
                      << ", type1=" << pair.second.type1
                      << ", type2=" << pair.second.type2 << std::endl;
        }
    }


    // Method to increment nbonds for a pair of tags
    //void incrementNbonds(unsigned int tag_i, unsigned int tag_j) {
    //    //auto it = dataMap.find(std::make_pair(tag_i, tag_j));
    //    auto it = dataMap.find(makeCanonicalKey(tag_i, tag_j));
    //    if (it != dataMap.end()) {
    //        it->second.nbonds++;  // Increment nbonds
    //    } else {
    //        // Handle error case, such as key not found
    //        std::cerr << "Key (" << tag_i << ", " << tag_j << ") not found in BondMap." << std::endl;
    //        }
    //    }

    // count keys where tag_i = tag_val 
    int countTagIKeys(unsigned int tag_val) const {
        int count = 0;
        for (const auto& pair : dataMap) {
            if (pair.first.first == tag_val) {
                ++count;
            }
            if (pair.first.second == tag_val) {
                ++count;
            }
        }
        return count;
    }

    // Method to search for all keys containing a specific tag_i and collect BondTrackData into a vector
    //std::vector<BondTrackData> getDataForTag(unsigned int tag_i) const {
    //    std::vector<BondTrackData> result;
    std::vector<std::pair<std::pair<unsigned int, unsigned int>, BondTrackData>> getDataForTag(unsigned int tag_i) const {
        std::vector<std::pair<std::pair<unsigned int, unsigned int>, BondTrackData>> result;
        for (const auto& pair : dataMap) {
            if (pair.first.first == tag_i || pair.first.second == tag_i) {
                result.push_back(pair);
                //result.push_back(pair.second);
            }
        }
        return result;
    }


    // Method to retrieve all keys currently in the BondMap
    std::vector<std::tuple<unsigned int, unsigned int>> getAllKeys() const {
        std::vector<std::tuple<unsigned int, unsigned int>> keys;
        for (const auto& pair : dataMap) {
            keys.push_back(pair.first);
        }
        return keys;
    }


    // Method to remove data from the map for a specific key
    bool removeKey(unsigned int tag_i, unsigned int tag_j) {
        auto it = dataMap.find(makeCanonicalKey(tag_i, tag_j));
        if (it != dataMap.end()) {
            dataMap.erase(it);
            //std::cout << "Key (" << std::get<0>(key) << ", " << std::get<1>(key) << ") removed from BondDataManager." << std::endl;
            //std::cout << "Removed (" << tag_i << "," << tag_j << ")" << std::endl;
            return true;
        } else {
            //std::cout << "Key (" << std::get<0>(key) << ", " << std::get<1>(key) << ") not found in BondDataManager." << std::endl;
            return false;
        }
    }


};

extern BondDataManager manager;  // Declaration of a global instance of BondDataManager

#endif // __BONDMAP_H__

