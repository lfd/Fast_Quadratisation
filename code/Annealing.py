# Copyright 2025 OTH - Laboratory for Digitalisation (LfD)
# Written by Lukas Schmidbauer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np
import random

def evalPBF(pbf, varAssignement):
    """
    Input:  pbf: a dictionary like function: {(1,2,4): 3.4, (3,): 7, ...} <-> 3.4 * x_1x_2x_4 + 7x_3
            varAssignement: Assignement of variables as a dictionary: {1:0, 2:0, 3:1, 4:0, ...}

    returns: Evaluated pbf as a double
    """

    out = 0
    for monomial in pbf.keys():
        impacts = True
        for var in monomial:
            val = varAssignement[var]
            if val == 0:
                impacts = False
                break

        if impacts:
            out += pbf[monomial]
    return out

def getPBFLandscape(pbf):
    """
    Input:  pbf: a dictionary like function: {(1,2,4): 3.4, (3,): 7, ...} <-> 3.4 * x_1x_2x_4 + 7x_3

    returns: Evaluated pbf at every possible value
    """
    varAssignement = getInitialVarAssignement(pbf)
    varList = []
    for var in varAssignement:
        varAssignement[var] = 0
        varList.append(var)
    
    return getPBFLandscapeHelper(pbf, varAssignement, varList, 0)

def getPBFLandscapeHelper(pbf, varAssignement, varList, varPos):
    if varPos >= len(varList):
        val = evalPBF(pbf, varAssignement)
        varAssignementStr = "".join(str(v) for k,v in sorted(varAssignement.items()))
        return [(varAssignementStr, val)]
    else:
        varAssignement[varList[varPos]] = 0
        left = getPBFLandscapeHelper(pbf, varAssignement, varList, varPos +1)
        varAssignement[varList[varPos]] = 1
        right = getPBFLandscapeHelper(pbf, varAssignement, varList, varPos +1)
        return left + right

def getInitialVarAssignement(pbf, seed=42):
    """
    Generate a random initial variable assignement.

    Input:  pbf: a dictionary like function: {(1,2,4): 3.4, (3,): 7, ...} <-> 3.4 * x_1x_2x_4 + 7x_3

    returns: Assignement of variables as a dictionary: {1:0, 2:0, 3:1, 4:0, ...}
    """
    random.seed(seed)
    varAssignement = dict()
    for monomial in pbf:
        for var in monomial:
            if not varAssignement.__contains__(var):
                varAssignement[var] = random.randint(0,1)
    
    return varAssignement

def allowAnyway(delta, temperature):
    """
    Performs random experiment for simulated annealing, depending on the current temperature and the delta energy
    """
    if temperature <= 0:
        return False
    prob = np.min([np.exp((-delta)/temperature), 1])
    return True if random.random() < prob else False


def simulatedAnnealing(pbf, steps, initialTemperature = 1000, seed=42):
    """
    Performs a linear temperature schedule of simulated annealing and returns a variable assignement and the energy value
    Input:  pbf: a dictionary like function: {(1,2,4): 3.4, (3,): 7, ...} <-> 3.4 * x_1x_2x_4 + 7x_3
            steps: number of iterations
            initialTemperature: higher energy -> higher probability to switch initially even if it does yield a better solution
    """

    varAssignement = getInitialVarAssignement(pbf, seed)
    random.seed(seed)
    varsList = [k for k in varAssignement.keys()]
    former_energy = evalPBF(pbf, varAssignement)
    
    for step in range(1, steps +1):
        T = -(step - steps) * initialTemperature

        var = varsList[random.randint(0, len(varsList)-1)]
        if varAssignement[var] == 0:
            varAssignement[var] = 1
        else:
            varAssignement[var] = 0 
        new_energy = evalPBF(pbf, varAssignement)

        if former_energy >= new_energy:
            former_energy = new_energy 
        elif allowAnyway(new_energy - former_energy, T):
            former_energy = new_energy
        else:
            if varAssignement[var] == 0:
                varAssignement[var] = 1
            else:
                varAssignement[var] = 0
        
    varAssignementStr = "".join(str(v) for k,v in sorted(varAssignement.items()))
    return varAssignementStr, varAssignement, former_energy