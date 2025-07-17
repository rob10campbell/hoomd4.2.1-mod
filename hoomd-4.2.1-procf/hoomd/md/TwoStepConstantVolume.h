// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// ########## Modified by Rheoinformatic //~ [RHEOINF] ##########

#ifndef HOOMD_TWOSTEPCONSTANTVOLUME_H
#define HOOMD_TWOSTEPCONSTANTVOLUME_H

#include "ComputeThermo.h"
#include "IntegrationMethodTwoStep.h"
#include "Thermostat.h"
#include "hoomd/Variant.h"
#include <pybind11/pybind11.h>
namespace hoomd::md
    {

/** Perform constant volume simulation.

    Implement the the Velocity-Verlet integration scheme with an optional velocity rescaling
    Thermostat.
*/
class PYBIND11_EXPORT TwoStepConstantVolume : public IntegrationMethodTwoStep
    {
    public:
    /** Construct the constant volume integration method.

        @param sysdef System to work on.
        @param group Subset of particles to integrate.
        @param thermostat Thermostat to use. Set to null for constant energy simulations.
    */
    TwoStepConstantVolume(std::shared_ptr<SystemDefinition> sysdef,
                          std::shared_ptr<ParticleGroup> group,
                          std::shared_ptr<Thermostat> thermostat)
        : IntegrationMethodTwoStep(sysdef, group), m_thermostat(thermostat), m_use_walls(false), m_wall("") //~ add walls [RHEOINF]
        {
        }

    virtual ~TwoStepConstantVolume() { }

    /** Performs the first half-step of the integration.

        @param timestep Current simulation timestep.

        @post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2 per
        the Velocity-Verlet method.
    */
    virtual void integrateStepOne(uint64_t timestep);

    /** Performs the second half-step of the integration.

        @param timestep Current simulation timestep.

        @post Particle velocities are moved forward to timestep+1.
    */
    virtual void integrateStepTwo(uint64_t timestep);

    /** Set the thermostat.

        @param thermostat Thermostat to use for velocity rescaling (may be null).
    */
    void setThermostat(std::shared_ptr<Thermostat> thermostat)
        {
        m_thermostat = thermostat;
        }

    //~ add walls [RHEOINF]
    /** Set the wall option

        @param wall_axes One of 3 wall axes options for bounceback: +x-x, +y-y, +z-z
 
    */
    /* // IF wall_axes IS REQUIRED
     * void setWall(const std::string& wall_axes)
        {
        static const std::set<std::string> valid_axes = {"+x-x", "+y-y", "+z-z"};
	if (wall_axes.empty())
	    {
            m_use_walls = false;
    	    }
	else if (valid_axes.find(wall_axes) == valid_axes.end())
            {
            throw std::invalid_argument("Invalid wall_axes option: " + wall_axes);
            m_use_walls = false;
            }
	else
	    {
            m_wall = wall_axes;
            m_use_walls = true;
    	    }
        }

    /// Get the wall option
    std::string getWall() const
        {
        if (m_use_walls)
            {
            return m_wall;
            }
        else
            {
            return "";
            }
        }
    */
    void setWall(const pybind11::object wall_axes)
    {
        static const std::set<std::string> valid_axes = {"+x-x", "+y-y", "+z-z"};

        if (wall_axes.is_none())
        {
            m_use_walls = false;
        }
        else
        {
            std::string wall_str = pybind11::cast<std::string>(wall_axes);
            if (valid_axes.find(wall_str) == valid_axes.end())
            {
                throw std::invalid_argument("Invalid wall_axes option: " + wall_str);
                m_use_walls = false;
            }
            else
            {
                m_wall = wall_str;
                m_use_walls = true;
            }
        }
    }

    std::string getWall() const
    {
        if (m_use_walls)
        {
            return m_wall;
        }
        else
        {
            return "";
        }
    }
    //~


    /** Set the distance limit applied to particles.

        @param limit Largest distance particles should move in a given timestep (may be null).

        When non-null, simulations do not conserve energy or momentum, but some users find the limit
        allows them to relax high energy initial conditions.
    */
    void setLimit(std::shared_ptr<Variant>& limit)
        {
        m_limit = limit;
        }

    /// Get the current limit.
    [[nodiscard]] auto getLimit() const
        {
        return m_limit;
        }

    /// Get needed pdata flags.
    virtual PDataFlags getRequestedPDataFlags()
        {
        PDataFlags flags;

        // Compute rotational kinetic energy in simulations with anisotropic degrees of freedom.
        if (m_aniso)
            {
            flags[pdata_flag::rotational_kinetic_energy] = 1;
            }
        return flags;
        }

    protected:
    /// Pack the limit values for use in the GPU kernel.
    auto getKernelLimitValues(uint64_t timestep)
        {
        const auto use_limit = static_cast<bool>(m_limit);
        Scalar maximum_displacement = use_limit ? m_limit->operator()(timestep) : 0.0;
        return std::make_pair(use_limit, maximum_displacement);
        }

    /// The thermostat to apply (may be null).
    std::shared_ptr<Thermostat> m_thermostat;

    //~ [RHEOINF]
    /// flag to enable walls in the y-direction 
    bool m_use_walls;
    /// the axes where walls are applied
    std::string m_wall;
    //~

    /// The distance limit to apply (may be null).
    std::shared_ptr<Variant> m_limit;
    };

    } // namespace hoomd::md

#endif // HOOMD_TWOSTEPCONSTANTVOLUME_H
