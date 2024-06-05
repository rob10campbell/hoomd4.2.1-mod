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

* [ ] `hoomd/`
    * [x] `example_plugins/`
        * [x] `pair_plugin/`
            * [x] EvaluatorPairExample.h : **pass diameter and typeIDs**
    * [ ] `md/`
        * [x] EvaluatorPairBuckingham.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairDLVO.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairDPDThermoDPD.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairDPDThermoDPDMorse.h : **pass diameter and typeIDs, poly param, poly defaults to mono**
        * [x] EvaluatorPairDPDThermoLJ.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairEwald.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairExpandedGaussian.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairExpandedLJ.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairExpandedMie.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairForceShiftedLJ.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairFourier.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairGauss.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairLJ.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairLJ0804.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairLJ1208.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairLJGauss.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairMie.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairMoliere.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairMorse.h : **pass diameter and typeIDs, poly param**
        * [x] EvaluatorPairOPP.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairReactionField.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairTWF.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairTable.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairYukawa.h : **pass diameter and typeIDs**
        * [x] EvaluatorPairZBL.h : **pass diameter and typeIDs**
        * [x] EvaluatorWalls.h : **pass diameter and typeIDs**
        * [ ] `pair/`
            * [ ] pair.py : **change on/off params: poly param [PotentialPairDPDThermo], poly param [Morse, DPDMorse], contact force [Morse], scaled D0 [Morse, DPDMorse]**
        * [x] PotentialPair.h : **track diameter and typeIDs**
        * [x] PotentialPairAlchemical.h : **track typeIDs**
        * [ ] PotentialPairDPDThermo.h : **track diameter, on/off poly param**
