# Checklist for Polydispersity Update!

[this is a test change]

### **In all evaluators:**

1: label the "use diameter" section "[PROCF2023]" (see EvaluatorPairYukawa.h for example)

```c++
    //!~ don't need diameter [PROCF2023]
    DEVICE static bool needsDiameter()
        {
        return false;
        }
    //! Accept the optional diameter values
    /*! \param di Diameter of particle i
        \param dj Diameter of particle j
    */
    DEVICE void setDiameter(Scalar di, Scalar dj) { }
    //~ [PROCF2023]
```

2: remove "contact" (but keep typeid)


### **In EvaluatorPairMorse**

1: add contact force <br>
2: add contact force boolean flag <br>
3: add scaled-D0 boolean flag <br>
4: replace poly parameter with boolean flag


### **In EvaluatorPairDPDThermoDPDMorse**

1: replace the `radsum = contact` with `radsum = 0.5 * (diameter_i + diameter_j);`
3: [optional] add contact force boolean flag <br>
4: add scaled-D0 boolean flag <br>
5: replace poly parameter with boolean flag <br>
6: [optional] remove a1 and a2?


### **In PotentialPair.h, PotentialPairDPDThermo.h, and PotentialPairAlchemical.h**

1: remove contact <br>
2: remove on-off poly param


### **In pair.py**
1: remove on-off poly param and poly parameter, and replace with boolean flag <br>
2: add on/off contact force flag <br>
3: add on/off scaled-D0 flag <br>
4: [optional] remove a1 and a2?


## **File List** (x when changed)
- **scaled_D0**: scaled D0 by particles size ((radius_i + radius_j)/2) for AO-style multimodal depletion potential; activated by an optional boolean (true/false) flag
- **f_contact=0**: removes contact force (and uses built-in Morse repulsion) when f_contact = 0
- **diameter**: adds diameter back (removed by HOOMD-blue devs between hoomdv3 and hoomdv4)
- **typeIDs**: tracks particle typeID to reset solvent radius to zero as needed in DPD force calcs

* [x] `hoomd/`
    * [x] `example_plugins/`
        * [x] `pair_plugin/`
            * [x] EvaluatorPairExample.h : **diameter, typeIDs**
    * [x] `md/`
        * [x] EvaluatorPairBuckingham.h : **diameter, typeIDs**
        * [x] EvaluatorPairDLVO.h : **diameter, typeIDs**
        * [x] EvaluatorPairDPDThermoDPD.h : **diameter, typeIDs**
        * [x] EvaluatorPairDPDThermoDPDMorse.h : **diameter, typeIDs, f_contact=0, scaled_D0**
        * [x] EvaluatorPairDPDThermoLJ.h : **diameter, typeIDs**
        * [x] EvaluatorPairEwald.h : **diameter, typeIDs**
        * [x] EvaluatorPairExpandedGaussian.h : **diameter, typeIDs**
        * [x] EvaluatorPairExpandedLJ.h : **diameter, typeIDs**
        * [x] EvaluatorPairExpandedMie.h : **diameter, typeIDs**
        * [x] EvaluatorPairForceShiftedLJ.h : **diameter, typeIDs**
        * [x] EvaluatorPairFourier.h : **diameter, typeIDs**
        * [x] EvaluatorPairGauss.h : **diameter, typeIDs**
        * [x] EvaluatorPairLJ.h : **diameter, typeIDs**
        * [x] EvaluatorPairLJ0804.h : **diameter, typeIDs**
        * [x] EvaluatorPairLJ1208.h : **diameter, typeIDs**
        * [x] EvaluatorPairLJGauss.h : **diameter, typeIDs**
        * [x] EvaluatorPairMie.h : **diameter, typeIDs**
        * [x] EvaluatorPairMoliere.h : **diameter, typeIDs**
        * [ ] EvaluatorPairMorse.h : **diameter, typeIDs, f_contact=0, scaled_D0**
        * [x] EvaluatorPairOPP.h : **diameter, typeIDs**
        * [x] EvaluatorPairReactionField.h : **diameter, typeIDs**
        * [x] EvaluatorPairTWF.h : **diameter, typeIDs**
        * [x] EvaluatorPairTable.h : **diameter, typeIDs**
        * [x] EvaluatorPairYukawa.h : **diameter, typeIDs**
        * [x] EvaluatorPairZBL.h : **diameter, typeIDs**
        * [x] EvaluatorWalls.h : **diameter, typeIDs**
        * [x] `pair/`
            * [x] pair.py : **optional scaled_D0 [Morse, DPDMorse]**
        * [x] PotentialPair.h : **track diameter, typeIDs**
        * [x] PotentialPairAlchemical.h : **typeIDs**
        * [x] PotentialPairDPDThermo.h : **diameter, typeIDs**
