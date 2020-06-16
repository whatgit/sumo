/****************************************************************************/
// Eclipse SUMO, Simulation of Urban MObility; see https://eclipse.org/sumo
// Copyright (C) 2001-2019 German Aerospace Center (DLR) and others.
// This program and the accompanying materials
// are made available under the terms of the Eclipse Public License v2.0
// which accompanies this distribution, and is available at
// http://www.eclipse.org/legal/epl-v20.html
// SPDX-License-Identifier: EPL-2.0
/****************************************************************************/
/// @file    MSCFModel_Heterogeneous.h
/// @author
/// @date    May 2020
/// @version $Id$
///
// The Heterogeneous car-following model
/****************************************************************************/
#ifndef MSCFMODEL_HETEROGENEOUS_H
#define MSCFMODEL_HETEROGENEOUS_H

// ===========================================================================
// included modules
// ===========================================================================
#include <config.h>

#include "MSCFModel.h"
#include <microsim/MSLane.h>
#include <microsim/MSVehicle.h>
#include <microsim/MSVehicleType.h>
#include <utils/xml/SUMOXMLDefinitions.h>


// ===========================================================================
// class definitions
// ===========================================================================
/** @class MSCFModel_IDM
 * @brief The Intelligent Driver Model (IDM) car-following model
 * @see MSCFModel
 */
class MSCFModel_Heterogeneous : public MSCFModel {
public:
    /** @brief Constructor
     *  @param[in] vtype the type for which this model is built and also the parameter object to configure this model
     *  @param[in] idmm Wether IDM or IDMM shall be built
     */
    MSCFModel_Heterogeneous(const MSVehicleType* vtype);


    /** @brief Constructor
     * @param[in] accel The maximum acceleration
     * @param[in] decel The maximum deceleration
     * @param[in] emergencyDecel The maximum emergency deceleration
     * @param[in] apparentDecel The deceleration as expected by others
     * @param[in] headwayTime the headway gap
     * @param[in] adaptationFactor a model constant
     * @param[in] adaptationTime a model constant
     * @param[in] internalStepping internal time step size
     */
    MSCFModel_Heterogeneous(const MSVehicleType* vtype, double accel, double decel, double emergencyDecel, double apparentDecel,
                  double headwayTime, double adaptationFactor, double adaptationTime,
                  double thwCC, double thwCT, double thwTC, double thwTT,
                  double internalStepping);


    /// @brief Destructor
    ~MSCFModel_Heterogeneous();


    /// @name Implementations of the MSCFModel interface
    /// @{

    /** @brief Applies interaction with stops and lane changing model influences
     * @param[in] veh The ego vehicle
     * @param[in] vPos The possible velocity
     * @return The velocity after applying interactions with stops and lane change model influences
     */
    double finalizeSpeed(MSVehicle* const veh, double vPos) const;

    /** @brief Computes the vehicle's safe speed without a leader
     *
     * Returns the velocity of the vehicle in dependence to the length of the free street and the target
     *  velocity at the end of the free range. If onInsertion is true, the vehicle may still brake
     *  before the next movement.
     * @param[in] veh The vehicle (EGO)
     * @param[in] speed The vehicle's speed
     * @param[in] seen The look ahead distance
     * @param[in] maxSpeed The maximum allowed speed
     * @param[in] onInsertion whether speed at insertion is asked for
     * @return EGO's safe speed
     */
    virtual double freeSpeed(const MSVehicle* const veh, double speed, double seen,
                             double maxSpeed, const bool onInsertion = false) const;

    /** @brief Computes the vehicle's safe speed (no dawdling)
     * @param[in] veh The vehicle (EGO)
     * @param[in] speed The vehicle's speed
     * @param[in] gap2pred The (netto) distance to the LEADER
     * @param[in] predSpeed The speed of LEADER
     * @return EGO's safe speed
     * @see MSCFModel::ffeV
     */
    double followSpeed(const MSVehicle* const veh, double speed, double gap2pred, double predSpeed, double predMaxDecel, const MSVehicle* const pred = 0) const;


    /** @brief Computes the vehicle's safe speed for approaching a non-moving obstacle (no dawdling)
     * @param[in] veh The vehicle (EGO)
     * @param[in] gap2pred The (netto) distance to the the obstacle
     * @return EGO's safe speed for approaching a non-moving obstacle
     * @see MSCFModel::ffeS
     * @todo generic Interface, models can call for the values they need
     */
    double stopSpeed(const MSVehicle* const veh, const double speed, double gap) const;


    /** @brief Returns the maximum gap at which an interaction between both vehicles occurs
     *
     * "interaction" means that the LEADER influences EGO's speed.
     * @param[in] veh The EGO vehicle
     * @param[in] vL LEADER's speed
     * @return The interaction gap
     * @todo evaluate signature
     * @see MSCFModel::interactionGap
     */
    double interactionGap(const MSVehicle* const, double vL) const;


    /** @brief Computes the vehicle's safe speed (no dawdling)
     * This method is used during the insertion stage. Whereas the method
     * followSpeed returns the desired speed which may be lower than the safe
     * speed, this method only considers safety constraints
     *
     * Returns the velocity of the vehicle in dependence to the vehicle's and its leader's values and the distance between them.
     * @param[in] veh The vehicle (EGO)
     * @param[in] speed The vehicle's speed
     * @param[in] gap2pred The (netto) distance to the LEADER
     * @param[in] predSpeed The speed of LEADER
     * @return EGO's safe speed
     */
    double insertionFollowSpeed(const MSVehicle* const veh, double speed, double gap2pred, double predSpeed, double predMaxDecel, const MSVehicle* const pred = 0) const;


    /** @brief Returns the minimum gap to reserve if the leader is braking at maximum (>=0)
     * @param[in] veh The vehicle itself, for obtaining other values
     * @param[in] pred The leader vehicle, for obtaining other values
     * @param[in] speed EGO's speed
     * @param[in] leaderSpeed LEADER's speed
     * @param[in] leaderMaxDecel LEADER's max. deceleration rate
     */
    double getSecureGap(const MSVehicle* const veh, const MSVehicle* const pred, const double speed, const double leaderSpeed, const double leaderMaxDecel) const;

    /** @brief Returns the model's name
     * @return The model's name
     * @see MSCFModel::getModelName
     */
    int getModelID() const {
        return SUMO_TAG_CF_HTG;
    }
    /// @}



    /** @brief Duplicates the car-following model
     * @param[in] vtype The vehicle type this model belongs to (1:1)
     * @return A duplicate of this car-following model
     */
    MSCFModel* duplicate(const MSVehicleType* vtype) const;


    VehicleVariables* createVehicleVariables() const {
        if (myAdaptationFactor != 1.) {
            return new VehicleVariables();
        }
        return 0;
    }


private:
    class VehicleVariables : public MSCFModel::VehicleVariables {
    public:
        VehicleVariables() : levelOfService(1.) {}
        /// @brief state variable for remembering speed deviation history (lambda)
        double levelOfService;
    };


private:
    double _v(const MSVehicle* const veh, const double gap2pred, const double mySpeed,
              const double predSpeed, const double desSpeed, const bool respectMinGap = true) const;


private:

    /// @brief The IDM delta exponent
    const double myDelta;

    /// @brief The IDMM adaptation factor beta
    const double myAdaptationFactor;

    /// @brief The IDMM adaptation time tau
    const double myAdaptationTime;

    /// @brief The number of iterations in speed calculations
    const int myIterations;

    /// @brief A computational shortcut
    const double myTwoSqrtAccelDecel;

    /// @brief A time headway for car following car situation
    const double myHeadwayCC;

    /// @brief A time headway for car following truck situation
    const double myHeadwayCT;

    /// @brief A time headway for truck following car situation
    const double myHeadwayTC;

    /// @brief A time headway for truck following truck situation
    const double myHeadwayTT;



private:
    /// @brief Invalidated assignment operator
    MSCFModel_Heterogeneous& operator=(const MSCFModel_Heterogeneous& s);
};

#endif /* MSCFMODEL_HETEROGENEOUS_H */
