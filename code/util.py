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
from quark import PolyBinary
import copy
from multiset import Multiset
from itertools import combinations
import pandas as pd
from LRQAOA import *
from qiskit.primitives import StatevectorSampler

def readDIMACSCNF(filename):
    """
    Input: filename: Name or path to file that encodes a SAT formula in DIMACS CNF format (see https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/satformat.ps)
    returns: python SAT formula in list of lists form (inner list corresponds to clause)
    """

    with open(filename, "r") as f:
        lines = f.readlines()

    v = int(lines[0].replace("\n","").replace("\r","").split(" ")[2])
    c = int(lines[0].replace("\n","").replace("\r","").split(" ")[3])
    
    formula = []
    avg_k = 0
    for i in range(1, len(lines)):
        if any(c.isdigit() and c != "0" for c in lines[i]):
            C = lines[i].rstrip("0 \n").replace("-", "!").split(" ")
            avg_k += len(C)
            formula.append(C)
    
    avg_k /= len(formula)

    return avg_k, v, c, formula

def consecutiveVariableSubstitution(SAT_instance):
    """
    Replaces variables in a SAT instance by consecutive variables starting from 1
    """
    origV = getSATVars(SAT_instance)

    mapping = dict( (str(v),str(k)) for k,v in enumerate(sorted(origV), start=1))

    ret_formula = []
    for Clause in SAT_instance:
        lC = []
        for lit in Clause:
            if "!" in lit:
                lit = "!" + mapping[lit.replace("!","")]
            else:
                lit = mapping[lit]
            lC.append(lit)
        ret_formula.append(lC)

    return ret_formula

def positiveVariableOptimisation(SAT_instance, varPrefix=""):
    """
    If clause C_i has 5 positive variables, a single positive variable is replaced by a new negative variable (in each clause it occurs in).
    Additionally, the new negative variable is constrained with the old positive variable by xor.
    In general: 
    If clause C_i has t > 5 positive variables, t-5 positive variables are replaced by t-5 new negative variables (in each clause they occur in).
    They are then each constrained analogously via xor.
    """

    local_instance = copy.deepcopy(SAT_instance)

    while not checkVariableOptimisationDone(local_instance):
        for C in local_instance:
            pLs = getPositiveLiterals(C)
            if len(pLs) >= 5:            
                vars = getSATVars(local_instance)
                n = np.max([int(v.replace("x", "")) for v in vars], initial=0) + 1
                pL = pLs.pop()

                for C2 in local_instance:
                    if pL in C2:
                        C2.remove(pL)
                        C2.append("!" + varPrefix + str(n))

                local_instance.append([pL, varPrefix + str(n)])
                local_instance.append(["!" + pL, "!" + varPrefix + str(n)])

                break

    return local_instance

def checkVariableOptimisationDone(SAT_instance):
    """
    Are there clauses in SAT_instance that contain 5 or more positive variables?
    """
    for C in SAT_instance:
        if len(getPositiveLiterals(C)) >= 5:
            return False
    return True

def getMaxPositiveLiteralsInClauses(SAT_instance):
    """
    Returns the maximum number of positive variables in clauses of a SAT instance
    """
    pV = 0
    for C in SAT_instance:
        pV = max(pV, len(getPositiveLiterals(C)))
    return pV

def getPositiveLiterals(clause):
    out = set()
    for l in clause:
        if "!" not in l:
            out.add(l)

    return out

def getNegativeLiterals(clause):
    out = set()
    for l in clause:
        if "!" in l:
            out.add(l)

    return out

def getPositiveLiteralsFormula(formula):
    out = set()
    for clause in formula:
        for l in clause:
            if "!" not in l:
                out.add(l)

    return out

def getNegativeLiteralsFormula(formula):
    out = set()
    for clause in formula:
        for l in clause:
            if "!" in l:
                out.add(l)

    return out

def getClauseVars(clause):
    """
    Returns set of used variables from Clause
    """
    vars = set()
    for lit in clause:
        vars.add(lit.replace("!", ""))
    return vars

def getSATVars(SAT_instance):
    """
    Returns set of used variables from SAT instance
    """
    vars = set()
    for C in SAT_instance:
        for lit in C:
            vars.add(lit.replace("!", ""))
    return vars

def getSATAvgK(SAT_instance):
    """
    Returns the average k in a given SAT_instance
    """

    avg_k = 0
    for c in SAT_instance:
        avg_k += len(c)
    avg_k /= len(SAT_instance)

    return avg_k

def getBorosPenaltyScaling(polynomial):
    """
    Returns 1 + 2*sum(abs(alpha_S))
    """
    return 1 + 2*np.sum(np.abs(list(polynomial.values())))

def getDensities(poly):
    """
    Input:
        poly: a dictionary of the form: {(<tuple of variables>): <alpha>, ...}
    Returns all degree-k densities of poly (i.e. actual / possible monomials of degree-k) in the form of a dictionary
        {<k>: <deg-k-density>}
    """

    variables = set()
    degree_dict = dict()
    for monomial in poly.keys():
        for v in monomial:
            variables.add(v)
        current_degree = len(monomial)
        if degree_dict.__contains__(current_degree):
            degree_dict[current_degree].append((poly[monomial], monomial))
        else:
            degree_dict[current_degree] = [(poly[monomial], monomial)]

    densities_dict = dict()
    for i in range(max(len(degree_dict.keys()) + 1, 5)):
        if degree_dict.__contains__(i):
            densities_dict[i] = len(degree_dict[i]) / scipy.special.comb(
                len(variables), i, exact=True
            )
        else:
            densities_dict[i] = 0

    return densities_dict

def recursiveExpand(input:list, outDict:dict, path=[]):
    """
    Recursively expands a term of the form (1-x1)(1-x2)...(1-xn) and writes result in outDict in the form {<monomial>: <alpha>}
    Input: List of variables to expand, for example ["x1", "x3", ...]. They are expanded, assuming above form: (1-x1)...
    Care when using this method with PolyBinary! -> Exponents are ignored and terms are added and subtracted there.
    """
    if len(input) <= 0: 
        if len(path) % 2 == 0:
            outDict.update({tuple(path): 1})
        else:
            outDict.update({tuple(path): -1})
    else:
        nv = input[0]
        p1 = path + [nv]
        p2 = path

        recursiveExpand(input[1:], outDict, p1)
        recursiveExpand(input[1:], outDict, p2)



def PUBOfromSAT(SAT_instance):
    """
    Creates a PUBO from a SAT_instance, via extended DeMorgan's law
    """
    vars = set()
    for C in SAT_instance:
        for lit in C:
            vars.add(lit.replace("!", ""))

    pubo = PolyBinary()
    for C in SAT_instance:
        negativeVars = set()
        positiveVars = set()
        for lit in C:
            if "!" in lit:
                negativeVars.add(int(lit.replace("x", "").replace("!","")))
            else:
                positiveVars.add(int(lit.replace("x", "")))
            
        positiveExpandDict = dict()
        recursiveExpand(list(positiveVars), positiveExpandDict)

        puboDict = dict()
        for key in positiveExpandDict.keys():
            puboDict[tuple(sorted(tuple(negativeVars) + key))] = positiveExpandDict[key]
        

        pubo += 1 - PolyBinary(puboDict)
    
    return cleanPolyBinary(PolyBinary({}) - pubo) 

def cleanPolyBinary(polynomial: PolyBinary):
    """
    Removes all monomials with 0 coefficients
    """
    wDict = {}

    for monomial, alpha in polynomial.items():
        if alpha != 0:
            wDict[monomial] = alpha

    return PolyBinary(wDict)

def polyDictToQubo(polynomial, compactify=True):
    """
    Return upper triangular qubo representation of quadratic polynomial (PolyBinary)
    """
    if compactify:
        polynomial = polynomial.compact()

    qubo = np.zeros(shape=(len(polynomial.variables), len(polynomial.variables)))
    for m in polynomial.keys():
        if len(m) > 2:
            print("Supplied higher order monomial! -> Ignoring")
        elif len(m) == 2:
            qubo[m[0]][m[1]] = polynomial[m]
        elif len(m) == 1:
            qubo[m[0]][m[0]] = polynomial[m]
        else:
            print("Supplied constant monomial! -> Ignoring")
    
    return qubo


def qiskitCountsToDict(counts, ising, pubo, qubo = None):
    """
    Converts Qiskit's measurement results into a dictionary of the form {"bitstr": [<list of strings>], "energy": [<list of energies>]}
    Note that we invert the bitstring, due to qiskits reverse qubit order.
    This function assumes consecutive numerical variable indices in pbf
    Input: 
        pbf: a quark pbf
        counts: a dictionary of the form {"bitstr": <count>, ...}
    """
    out = {}
    out["bitstr"] = []
    out["energy_current_pbf"] = []
    out["energy_original_pbf"] = []
    out["energy_ising"] = []
    for bitstr, cnt in counts.items():
        reverse_bitstr = bitstr[::-1]
        min_var_idx = min(pubo.variables)
        varA = {}
        for i in range(len(bitstr)):
            varA[i+min_var_idx] = int(reverse_bitstr[i])
        energyIsing = ising.evaluate({v: binToIsingValues(val) for v,val in varA.items()})
        energyPubo = pubo.evaluate(varA)
        energyQubo = qubo.evaluate(varA) if qubo != None else energyPubo
        for i in range(cnt):
            out["bitstr"] += [reverse_bitstr]
            out["energy_current_pbf"] += [energyQubo]
            out["energy_original_pbf"] += [energyPubo]
            out["energy_ising"] += [energyIsing]
    
    return out

def binToIsingValues(bin, inverted=True):
    """
    Transforms bin to ising values, according to:
    if inverted:
        0 <-> 1
        1 <-> -1
    else:
        0 <-> -1
        1 <-> 1
    """
    if inverted:
        return 1 if bin == 0 else -1
    else:
        return 1 if bin == 1 else -1

def doQAOASimulation(ising, pubo, qubo, dbeta, dgamma, p, shots, pbf_type, filename):
    """
    Returns a pandas dataframe for results obtained by statevector simulation
    ising: Quark Ising model that originates from pubo.to_ising()
    pubo: Quark Polynomial
    qubo: Quark Polynomial
    dbeta: float for LRQAOA mixer
    dgamma: float for LRQAOA problem hamiltonian
    p: LRQAOA layers
    shots: Number of shots for statevector sampler
    pbf_type: string for pbf type identification in dataframe
    filename: filename of sat formulation to join to other dataframes
    """
    if len(ising.variables) > 23:
        print("Info: Ising model for ", filename, "has ", len(ising.variables), " > 23 variables. No simulation executed!")
        return pd.DataFrame({"SAT_filename": filename, "QA_status": "failed", "QA_pbf_type": pbf_type}, index=[0])
    
    sampler = StatevectorSampler()
    circ = create_LR_QAOA_Circuit(ising, dbeta=dbeta, dgamma=dgamma, p=p)
    circ.measure_all()

    counts = sampler.run([circ], shots=shots).result()[0].data.meas.get_counts()
    circDf = pd.DataFrame(qiskitCountsToDict(counts=counts, ising=ising, pubo=pubo, qubo=qubo))
    circDf["SAT_filename"] = filename
    circDf["QA_pbf_type"] = pbf_type
    circDf["QA_p"] = p
    circDf["QA_shots"] = shots
    circDf["QA_dbeta"] = dbeta
    circDf["QA_dgamma"] = dgamma
    circDf["QA_status"] = "ok"

    return circDf


#######################################################################################
# Visuals
#######################################################################################
def multiGraphToLatexTikz(G, radius, clique_sort=False):
        """
        Input: networkx MultiGraph with compact node names
                radius: Node radius (typ 2, big graphs 10)
                clique_sort: Sorts nodes by clique size
        Returns: latex string of that drawn MultiGraph
        """

        oG = copy.deepcopy(G)
        nG = nx.Graph()
        nG.add_nodes_from(list(G.nodes()))
        nG.add_edges_from(list(G.edges()))
        
        out = "\\resizebox{\\linewidth}{!}{%\n\\begin{tikzpicture}\n\t\\begin{scope}[circle]\n\t\t\draw\n"
        
        size = G.number_of_nodes()
        i = 0

        node_colour = nx.get_node_attributes(oG, "colour", default="lfdblack")

        if clique_sort:
            nodes = []
            while len(G.nodes) > 0:
                n = sorted(list(nx.find_cliques(G)),key=len, reverse=True)[0]
                nodes += n
                G.remove_nodes_from(n)
                
            for n in nodes:
                out += "\t\t(" + str((i*360.0)/size) +  ":" + str(radius) + ") \tnode [" + node_colour[n] + "] (x" + str(n) + ") {$x_{" + str(n) + " }$}\n"
                i+=1
        else:
            for n in G.nodes():
                out += "\t\t(" + str((i*360.0)/size) +  ":" + str(radius) + ") \tnode [" + node_colour[n] + "] (x" + str(n) + ") {$x_{" + str(n) + " }$}\n"
                i+=1
        
        out += "\t\t;\n\t\\end{scope}\n"
        
        edge_colour = nx.get_edge_attributes(oG, "colour", default="lfdblack")
        edge_weight = nx.get_edge_attributes(oG, "weight", default=1)
        
        out += "\t\\begin{scope}[-]\n"
        for e in nG.edges():
            num_ed = edge_weight[e]
            if num_ed == 1:
                out += "\t\t\\draw [" + edge_colour[e] + "] (x" + str(e[0]) + ") to " + "(x" + str(e[1]) + ");\n"
            else:
                out += "\t\t\\draw [" + edge_colour[e] + "] (x" + str(e[0]) + ") to node[above, sloped] {$" + str(num_ed) + "$} (x" + str(e[1]) + ");\n"
            
        
        out += "\t\\end{scope}\n\end{tikzpicture}\n}"
        
        return out

def primalGraphFromPolynomial(polynomial, oldGraph=None):
    """
    Creates a primal (or variable-incidence) graph from a given polynomial (in Quark format).
    Variables are nodes and edges are introduced whenever two variables occur in the same monomial
    If oldGraph is provided, new edges and new nodes will be coloured red
    """
    edges = Multiset()
    vars = set()
    for m in polynomial.keys():
        sub_vars = set()
        for var in m:
            vars.add(var)
        for e in list(combinations(m, 2)):
            edges.add(e)
    
    G = nx.Graph()
    for n in vars:
        if oldGraph is not None and n not in oldGraph.nodes:
            G.add_node(n, colour="lfdred")
        else:
            G.add_node(n, colour="lfdblack")

    for e in edges:
        if oldGraph is not None and e not in oldGraph.edges:
            G.add_edge(e[0], e[1], weight=edges.get(e, 0), colour="lfdred")
        else:
            G.add_edge(e[0], e[1], weight=edges.get(e, 0), colour="lfdblack")

    
    return G

def primalGraphFromkSAT(SAT_instance, oldGraph=None):
    """
    Creates a primal (or variable-incidence) graph from a given SAT instance.
    Variables are nodes and edges are introduced whenever two variables occur in the same clause (doesn't matter if negated or not)
    If oldGraph is provided, new edges and new nodes will be coloured red
    """

    edges = []
    vars = set()
    for C in SAT_instance:
        sub_vars = set()
        for lit in C:
            vars.add(int(lit.replace("!", "").replace("x","")))
            sub_vars.add(int(lit.replace("!", "").replace("x","")))
        edges += combinations(sub_vars, 2)
    
    G = nx.Graph()
    for n in vars:
        if oldGraph is not None and n not in oldGraph.nodes:
            G.add_node(n, colour="lfdred")
        else:
            G.add_node(n, colour="lfdblack")

    for e in edges:
        if oldGraph is not None and e not in oldGraph.edges:
            G.add_edge(e[0], e[1], weight=edges.count(e), colour="lfdred")
        else:
            G.add_edge(e[0], e[1], weight=edges.count(e), colour="lfdblack")

    return G