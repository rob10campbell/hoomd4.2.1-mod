#ifndef __BONDMAP_H__
#define __BONDMAP_H__

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
    bool broken;
    uint64_t brokeTimestep;

    // Default constructor to initialize the members
    BondInfo()
        : timestep(0), id1(0), id2(0), typei(0), typej(0), formedTimestep(0), broken(false), brokeTimestep(0)
        {}
 };

class BondMap : public Compute
{
public:
    // Default Constructor
    //BondMap() : Compute(nullptr), num_colloid(0)
    BondMap() : Compute(nullptr), fileNameBondHist("BondHistory.csv")
    {
    // Do not perform any initialization that requires SystemDefinition

    unsigned int num_colloid = 1;

    // fill the map with dummy data (see function in protected section)
    populateDummyData();

    // Determine the file names dynamically
    stringstream bondHistFileName;
    bondHistFileName << "BondHistory";

    int file_counter = 0;
    while (true) {
        fileNameBondHist = bondHistFileName.str();
        if (file_counter > 0) {
            stringstream counter_broke;
            counter_broke << "_" << file_counter;
            fileNameBondHist += counter_broke.str();
        }
        fileNameBondHist += ".csv";

        // Check if the file exists
        ifstream file(fileNameBondHist.c_str());
        if (!file) {
            break; // Unique file name found, exit the loop
        }

        file_counter++;
    }

//#ifdef ENABLE_MPI
//    if(m_sysdef->isDomainDecomposed())
//        {
//        MPI_Allreduce(MPI_IN_PLACE,
//                      //&num_colloid,
//                      1,
//                      MPI_INT,
//                      MPI_SUM,
//                      m_exec_conf->getMPICommunicator());
//        }
//#endif

    // calculate the number of solvent particles
    //num_solvent = m_pdata->getNGlobal() - num_colloid;

    // find the number of available j colloids to pair with i
    //num_slot = num_colloid*(num_colloid-1)/2;

    // split tracking evenly across MPI processes (requires rank/elements/index calc)
 //   num_rank = m_exec_conf->getNRanks();
 //   my_rank = m_exec_conf->getRank();
    //num_ele_per_process = num_slot / num_rank;
    //index_start = my_rank * num_ele_per_process;
    //index_end = (my_rank+1) * num_ele_per_process;
    //if(my_rank == (num_rank-1))
    //    {
    //    index_end = num_slot;
    //    }
    //track_size = index_end - index_start;
    //Bond_track.resize(track_size);

    // create output files for recording bond data
//    if(my_rank == 0)
//        {
        BondFile.open(fileNameBondHist.c_str(), ios::out | ios::app); //"BondHistory.csv",ios::out | ios::trunc);
        BondFile << "timestep,tag_i,tag_j,typei,typej,formedTime,brokeTime" << "\n";
        BondFile.close();
//        }
    }

    // Constructor w/ SystemDefinition
    //BondMap(shared_ptr<SystemDefinition> sysdef) : Compute(sysdef)
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
    //    BondFile << "timestep,tag_i,tag_j,typei,typej,formedTime,brokeTime" << "\n";
    //    BondFile.close();
    //}

    // Destructor
    ~BondMap(){}

    // Methods to manage bond information
    void addBondFormed(uint64_t timestep, unsigned int id1, unsigned int id2, unsigned int typei, unsigned int typej)
    {
        auto it = bonds.find({id1, id2});
        if (it == bonds.end())
        {
            BondInfo bond;
            bond.timestep = timestep;
            bond.id1 = id1;
            bond.id2 = id2;
            bond.typei = typei;
            bond.typej = typej;
            bond.formedTimestep = timestep;
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
            BondFile << timestep << "," << id1 << "," << id2 << "," << bond_info.typei << "," << bond_info.typej << "," << bond_info.formedTimestep << "," << bond_info.brokeTimestep << "\n";
            BondFile.close();
        }
    }

    // record the lifetime at which a bond breaks/forms
//    virtual void writeBondtime()
//        {
//    #ifdef ENABLE_MPI
//        if (m_sysdef->isDomainDecomposed())
//            {
//            vector<vector<uint3> > BondHist_temp;
//            gather_v(BondHist, BondHist_temp, 0, m_exec_conf->getMPICommunicator());
//            BondHist.clear();
//
//            if(my_rank == 0)
//                {
//                BondFile.open(fileNameBondHist.c_str(), ios::out | ios::app); //"BondHistory.csv",ios::app);
//                for (unsigned int i=0; i < num_rank; i++)
//                    {
//                    for (unsigned int j=0; j < BondHist_temp[i].size(); j++)
//                        {
//                        BondFile << BondHist_temp[i][j].x << "," << BondHist_temp[i][j].y << "," << BondHist_temp[i][j].z << "\n";
//                        }
//                    }
//                BondHist_temp.clear();
//                BondFile.close();
//            }
//        else
//            {
//            BondFile.open(fileNameBondHist.c_str(), ios::out | ios::app); //"BondHistory.csv",ios::app);
//            for (unsigned int i=0; i < BondHist.size(); i++)
//                {
//                BondFile << BondHist[i].x << "," << BondHist[i].y << "," << BondHist[i].z << "\n";
//                }
//            BondHist.clear();
//            BondFile.close();
//            }
//    #endif
//        }

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

    unsigned int num_colloid; // the number of colloid particles
protected:
    //unsigned int num_colloid; // the number of colloid particles
    std::map<std::pair<unsigned int, unsigned int>, BondInfo> bonds; // Map to store bonds
    ofstream BondFile;
    //unsigned int num_solvent;         //!< the number of solvent particles
    //unsigned int num_slot;            //!< number of j colloids able to pair with a colloid particle i
 //   unsigned int num_rank;            //!< number of MPI ranks
 //   unsigned int my_rank;             //!< current rank
    //unsigned int num_ele_per_process; //!< number of elements per MPI process
    //unsigned int index_start;         //!< elements in this process (start)
    //unsigned int index_end;           //!< elements in this process (end)
    //unsigned int track_size;          //!< size of a process
    //map<unsigned int, unsigned char> Bond_check;
    //vector<unsigned int> Bond_track;
    vector<uint3> BondHist;
    std::string fileNameBondHist;        //!< File name for bond broke history

    // generate dummy data
    void populateDummyData()
    {
        for (unsigned int i = 0; i < 10; ++i)
        {
            for (unsigned int j = i + 1; j < 10; ++j)
            {
                BondInfo bond;
                bond.timestep = 1000 + i * 100 + j;
                bond.id1 = i;
                bond.id2 = j;
                bond.typei = i % 3;
                bond.typej = j % 3;
                bond.formedTimestep = 1000;
                bond.broken = (i + j) % 2 == 0; // Randomly mark some bonds as broken
                bond.brokeTimestep = bond.broken ? bond.timestep + 100 : 0;
                bonds[{i, j}] = bond;
            }
        }
    }

};

} // end namespace md
} // end namespace hoomd

#endif // __BONDMAP_H__

