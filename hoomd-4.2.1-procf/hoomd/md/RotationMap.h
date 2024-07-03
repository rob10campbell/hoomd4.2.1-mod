#ifndef __ROTATIONMAP_H__
#define __ROTATIONMAP_H__

#include "hoomd/GlobalArray.h"
//#include "hoomd/SystemDefinition.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/Compute.h"
#include "hoomd/HOOMDMPI.h"

#include <memory>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>

using namespace std;

namespace hoomd
{
namespace md
{

typedef double Scalar;

// Define a struct to hold bond information
struct BondInfo
{
    uint64_t timestep;
    unsigned int id1;
    unsigned int id2;
    unsigned int typei;
    unsigned int typej;
    uint64_t formedTimestep;
    Scalar theta_i_0;
    Scalar theta_j_0;
    Scalar gamma_ij_0;
    bool broken;
    uint64_t brokeTimestep;
};

class RotationMap : public Compute
{
public:
    // Default Constructor
    //RotationMap() : Compute(nullptr), num_colloid(0)
    RotationMap() : Compute(nullptr)
    {
    // Do not perform any initialization that requires SystemDefinition
    BondFile.open("BondHistory.csv", ios::out | ios::trunc);
    BondFile << "timestep,tag_i,tag_j,typei,typej,formedTime,theta_i_0,theta_j_0,gamma_ij_0,brokeTime" << "\n";
    BondFile.close();
    }

    // Constructor w/ SystemDefinition
    //RotationMap(shared_ptr<SystemDefinition> sysdef) : Compute(sysdef)
    //{
    //    ArrayHandle<Scalar4> h_pos2(m_pdata->getPositions(), access_location::host, access_mode::read);
    //    // start the counter at zero
    //    num_colloid = 0;
    //    // count the number of colloid particles
    //    for (unsigned int i = 0; i < m_pdata->getN(); i++)
    //    {
    //        unsigned int typei_temp = __scalar_as_int(h_pos2.data[i].w);
    //        if (typei_temp != 0)
    //        {
    //            num_colloid += 1;
    //        }
    //    }
    //
    //    BondFile.open("BondHistory.csv", ios::out | ios::trunc);
    //    BondFile << "timestep,tag_i,tag_j,typei,typej,formedTime,theta_i_0,theta_j_0,gamma_ij_0,brokeTime" << "\n";
    //    BondFile.close();
    //}

    // Destructor
    ~RotationMap(){}

    // Methods to manage bond information
    void addBondFormed(uint64_t timestep, unsigned int id1, unsigned int id2, unsigned int typei, unsigned int typej, Scalar theta_i_0, Scalar theta_j_0, Scalar gamma_ij_0)
    {
        auto it = bonds.find({id1, id2});
        if (it == bonds.end())
        {
            BondInfo bond{timestep, id1, id2, typei, typej, 0, theta_i_0, theta_j_0, gamma_ij_0, false, uint64_t(0.0)};
            bonds[{id1, id2}] = bond;
            bonds[{id1, id2}].formedTimestep = timestep;
        }
    }

    void addBondBroke(uint64_t timestep, unsigned int id1, unsigned int id2)
    {
        auto it = bonds.find({id1, id2});
        if (it != bonds.end())
        {
            bonds[{id1, id2}].broken = true;
            bonds[{id1, id2}].brokeTimestep = timestep;
        }
    }

    // Get the BondInfo for a given pair of IDs
    bool getBondInfo(unsigned int id1, unsigned int id2, BondInfo & bond_info)
    {
        auto it = bonds.find({id1, id2});
        if (it != bonds.end())
        {
            bond_info = it->second;
            return true; // Return true if bond info found
        }
        return false; // Return false if bond info not found
    }

    void writeBondHistory(uint64_t timestep, unsigned int id1, unsigned int id2)
    {
        auto it = bonds.find({id1, id2});
        if (it != bonds.end())
        {
            BondInfo bond_info = it->second;
            BondFile.open("BondHistory.csv", ios::out | ios::app);
            BondFile << timestep << "," << id1 << "," << id2 << "," << bond_info.typei << "," << bond_info.typej << "," << bond_info.formedTimestep << "," << bond_info.theta_i_0 << "," << bond_info.theta_j_0 << "," << bond_info.gamma_ij_0 << "," << bond_info.brokeTimestep << "\n";
            BondFile.close();
        }
    }

    bool isBondFormed(unsigned int id1, unsigned int id2)
    {
        auto it = bonds.find({id1, id2});
        return it != bonds.end() && !it->second.broken;
    }

    bool isBondBroke(unsigned int id1, unsigned int id2)
    {
        auto it = bonds.find({id1, id2});
        return it != bonds.end() && it->second.broken;
    }

    // Get the timestep when the bond was broken
    uint64_t getBondBrokeTimestep(uint64_t timestep, unsigned int id1, unsigned int id2)
    {
        return bonds[{id1, id2}].brokeTimestep;
    }

protected:
    unsigned int num_colloid; // the number of colloid particles
    std::map<std::pair<unsigned int, unsigned int>, BondInfo> bonds; // Map to store bonds
    ofstream BondFile;
};

} // end namespace md
} // end namespace hoomd

#endif // __ROTATIONMAP_H__

