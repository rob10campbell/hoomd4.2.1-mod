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
    * [ ] `example_plugins/`
        * [ ] `pair_plugin/`
            * [ ] EvaluatorPairExample.h : **pass diameter and typeIDs**
    * [ ] `md/`
        * [ ] EvaluatorPairBuckingham.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairDLVO.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairDPDThermoDPD.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairDPDThermoDPDMorse.h : **pass diameter and typeIDs, poly param, poly defaults to mono**
        * [ ] EvaluatorPairDPDThermoLJ.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairEwald.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairExpandedGaussian.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairExpandedLJ.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairExpandedMie.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairForceShiftedLJ.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairFourier.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairGauss.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairLJ.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairLJ0804.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairLJ1208.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairLJGauss.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairMie.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairMoliere.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairMorse.h : **pass diameter and typeIDs, poly param**
        * [ ] EvaluatorPairOPP.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairReactionField.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairTWF.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairTable.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairYukawa.h : **pass diameter and typeIDs**
        * [ ] EvaluatorPairZBL.h : **pass diameter and typeIDs**
        * [ ] EvaluatorWalls.h : **set diameter and typeIDs to zero, pass diameter and typeIDs**
        * [ ] `pair/`
            * [ ] pair.py : **change on/off params: poly param [PotentialPairDPDThermo], poly param [Morse, DPDMorse], contact force [Morse], scaled D0 [Morse, DPDMorse]**
        * [ ] PotentialPair.h : **track diameter and typeIDs**
        * [ ] PotentialPairAlchemical.h : ~**track contact dist**~
        * [ ] PotentialPairDPDThermo.h : **track diameter, on/off poly param**
