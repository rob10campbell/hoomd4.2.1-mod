# Checklist for Polydispersity Update!

[this is a test change]

### **In all evaluators:**

- [x] label the "use diameter" section "[PROCF2023]" (see EvaluatorPairYukawa.h for example)
- [x] remove "contact" (but keep typeid)


### **In EvaluatorPairMorse**

- [x] add contact force
- [x] f_contact=0 reverts to Morse built-in repulsion
- [x] add scaled_D0 boolean flag
- [x] mono and poly both use the real diameters in the system (not provided a1,a2) 
- [ ] [optional] remove a1 and a2?


### **In EvaluatorPairDPDThermoDPDMorse**

- [x] replace the `radsum = contact` with `radsum = 0.5 * (diameter_i + diameter_j);`
- [x] f_contact=0 reverts to Morse built-in repulsion
- [x] add scaled_D0 boolean flag
- [x] mono and poly both use the real diameters in the system (not provided a1,a2) 
- [ ] [optional] remove a1 and a2?


### **In PotentialPair.h, PotentialPairDPDThermo.h, and PotentialPairAlchemical.h**

- [x] remove contact distance calc
= [x] remove on-off poly param


### **In pair.py**
- [x] remove on-off poly param and poly parameter, and replace with optional scaled_D0 boolean flag
- [ ] [optional] remove a1 and a2?


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
