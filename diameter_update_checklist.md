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
- **scaled_D0**: scale D0 by particles size ((radius_i + radius_j)/2) for AO-style multimodal depletion potential; activated by an optional boolean (true/false) flag
- **f_contact=0**: when f_contact = 0, removes contact force (and uses built-in Morse repulsion, or DPD Conservative Force if D0=0)
- **radcontact**: passes the radsum from PotentialPairDPDThermo to DPDMorse; must be added to all Evaluator files but not used elsewhere (for some reason diameter did not work for the DPDMorse Evaluator)
- **diameter**: adds diameter back (removed by HOOMD-blue devs between hoomdv3 and hoomdv4)
- **typeIDs**: tracks particle typeID, used to reset solvent radius to zero in DPD force calcs

* [x] `hoomd/`
    * [x] `example_plugins/`
        * [x] `pair_plugin/`
            * [x] EvaluatorPairExample.h : **radcontact, diameter, typeIDs**
    * [x] `md/`
        * [x] EvaluatorPairBuckingham.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairDLVO.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairDPDThermoDPD.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairDPDThermoDPDMorse.h : **radcontact, diameter, typeIDs, f_contact=0, scaled_D0**
        * [x] EvaluatorPairDPDThermoLJ.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairEwald.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairExpandedGaussian.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairExpandedLJ.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairExpandedMie.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairForceShiftedLJ.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairFourier.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairGauss.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairLJ.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairLJ0804.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairLJ1208.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairLJGauss.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairMie.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairMoliere.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairMorse.h : **radcontact, diameter, typeIDs, f_contact=0, scaled_D0**
        * [x] EvaluatorPairOPP.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairReactionField.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairTWF.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairTable.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairYukawa.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorPairZBL.h : **radcontact, diameter, typeIDs**
        * [x] EvaluatorWalls.h : **radcontact, diameter, typeIDs**
        * [x] `pair/`
            * [x] pair.py : **optional scaled_D0 [Morse, DPDMorse]**
        * [x] PotentialPair.h : **radcontact, track diameter, typeIDs**
        * [x] PotentialPairAlchemical.h : **radcontact, typeIDs**
        * [x] PotentialPairDPDThermo.h : **radcontact, diameter, typeIDs**
