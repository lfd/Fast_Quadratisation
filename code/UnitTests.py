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

import math
import random
import time
import pandas as pd
import scipy.special

from p_tqdm import p_imap
from itertools import combinations
from quark import PolyBinary
from pathlib import Path
import shutil
import glob
from FastReduc import *
from util import *
from LRQAOA import *
from edge_coloring import *
from Annealing import *

EXP_INDUSTRY_DIRECTORY = "Experiments_Industrial_Test/"
EXP_GRAPH_DIRECTORY = "Graph_eval/"


def testDimacsGraph(clique_sort=False):
    instances = ["Instances/power_7k_13v_13c.cnf",
                 "Instances/power_5k_13v_37c.cnf"]
    for i in instances:
        testDimacsInstancesSingle(i, graph_eval=True)

    for dest in list(glob.glob(EXP_GRAPH_DIRECTORY + "/*.gml.bz2")):
        G = nx.read_gml(dest, destringizer=int)
        print(multiGraphToLatexTikz(G, 10, clique_sort), file=open(dest.replace(EXP_GRAPH_DIRECTORY, "pics/").replace(".gml.bz2", ".tex"), 'w'))

def testDimacsInstances(folder="Instances/"):
    """
    Tests pre generated SAT instances in the Dimacs format
    """
    filenames = glob.glob(folder + "*.cnf")

    shutil.rmtree(EXP_INDUSTRY_DIRECTORY)
    Path(EXP_INDUSTRY_DIRECTORY).mkdir(parents=True, exist_ok=True)
    pbf_landscape = [False for x in filenames]
    pbf_simulatedAnnealing = [False for x in filenames]
    qaoa_simulation = [False for x in filenames]
    completeExpDimacs = p_imap(testDimacsInstancesSingle, filenames, pbf_landscape, pbf_simulatedAnnealing, qaoa_simulation)

    pd.DataFrame(completeExpDimacs).to_csv(EXP_INDUSTRY_DIRECTORY + "Dimacs_results.csv")

def testDimacsInstancesSingle(filename, pbf_landscape=False, pbf_simulatedAnnealing=False, qaoa_simulation=False, graph_eval=False):
    """
    Loads a given k-SAT instance from file given by filename. 
    Creates two paths: (a) optimised and (b) non-optimised.
    From each SAT representation (a) or (b), we create a pubo. 
    It is either directly cast into QAOA circuit (a1, b1) or transformed to qubo and then transformed to a qaoa circuit (a2, b2).
    """

    DIres = {}
    DIres["SAT_filename"] = filename
    DIres["SAT_genVariant"] = filename.split("/")[1].split("_")[0]
    
    #k-SAT
    tmR1 = time.time()
    avg_k, v, c, formula = readDIMACSCNF(filename)
    tmR2 = time.time()
    formula = consecutiveVariableSubstitution(formula)
    tmR3 = time.time()
    DIres["SAT_gen_v"] = v
    v = len(getSATVars(formula))
    DIres["SAT_avgk"] = avg_k
    DIres["SAT_v"] = v
    DIres["SAT_c"] = len(formula)
    DIres["SAT_lit"] = len(getPositiveLiteralsFormula(formula)) + len(getNegativeLiteralsFormula(formula))
    DIres["SAT_pLit"] = len(getPositiveLiteralsFormula(formula))
    DIres["SAT_nLit"] = len(getNegativeLiteralsFormula(formula))
    DIres["SAT_time_read"] = tmR2 - tmR1
    DIres["SAT_time_consVarRep"] = tmR3 - tmR2

    if graph_eval:
        GnpVsat = primalGraphFromkSAT(formula)
        nx.write_gml(GnpVsat, EXP_GRAPH_DIRECTORY + Path(filename).stem + "_npVSat.gml.bz2")

    directSATviable = True if getMaxPositiveLiteralsInClauses(formula) <= 10 else False

    #Optimised k-SAT
    tmpVOpt1 = time.time()
    pVOptFormula = positiveVariableOptimisation(formula)
    tmpVOpt2 = time.time()
    DIres["pVSAT_avgk"] = getSATAvgK(pVOptFormula)
    DIres["pVSAT_v"] = len(getSATVars(pVOptFormula))
    DIres["pVSAT_c"] = len(pVOptFormula)
    DIres["pVSAT_lit"] = len(getPositiveLiteralsFormula(pVOptFormula)) + len(getNegativeLiteralsFormula(pVOptFormula))
    DIres["pVSAT_pLit"] = len(getPositiveLiteralsFormula(pVOptFormula))
    DIres["pVSAT_nLit"] = len(getNegativeLiteralsFormula(pVOptFormula))
    DIres["pVSAT_time_pVOpt"] = tmpVOpt2 - tmpVOpt1

    if graph_eval:
        GpVsat = primalGraphFromkSAT(pVOptFormula, GnpVsat)
        nx.write_gml(GpVsat, EXP_GRAPH_DIRECTORY + Path(filename).stem + "_pVSat.gml.bz2")

    if directSATviable:
        #Pubo from non optimised k-SAT
        tmnpVPubo1 = time.time()
        nonpVOptPubo = PUBOfromSAT(formula)
        tmnpVPubo2 = time.time()
        nonPVOptPenaltyScalingFactor = getBorosPenaltyScaling(nonpVOptPubo)
        nonpVOptPuboDensities = getDensities(nonpVOptPubo)
        DIres["npVPubo_boros_penalty_factor"] = nonPVOptPenaltyScalingFactor
        DIres["npVPubo_deg"] = nonpVOptPubo.degree
        DIres["npVPubo_num_vars"] = len(nonpVOptPubo.variables)
        DIres["npVPubo_num_monomials"] = len(nonpVOptPubo.keys())
        for deg, val in nonpVOptPuboDensities.items():
            DIres[("npVPubo_density_deg_" + str(deg))] = val
        DIres["npVPubo_time_PubofromSAT"] = tmnpVPubo2 - tmnpVPubo1
        if graph_eval:
            GnpVPubo = primalGraphFromPolynomial(nonpVOptPubo, GnpVsat)
            nx.write_gml(GnpVPubo, EXP_GRAPH_DIRECTORY + Path(filename).stem + "_npVPubo.gml.bz2")

    #Pubo from optimised k-SAT
    tmpVPubo1 = time.time()
    pVOptPubo = PUBOfromSAT(pVOptFormula)
    tmpVPubo2 = time.time()
    pVOptPenaltyScalingFactor = getBorosPenaltyScaling(pVOptPubo)
    pVOptPuboDensities = getDensities(pVOptPubo)
    DIres["pVPubo_boros_penalty_factor"] = pVOptPenaltyScalingFactor
    DIres["pVPubo_deg"] = pVOptPubo.degree
    DIres["pVPubo_num_vars"] = len(pVOptPubo.variables)
    DIres["pVPubo_num_monomials"] = len(pVOptPubo.keys())
    for deg, val in pVOptPuboDensities.items():
        DIres[("pVPubo_density_deg_" + str(deg))] = val
    DIres["pVPubo_time_PubofromSAT"] = tmpVPubo2 - tmpVPubo1

    if graph_eval:
        GpVPubo = primalGraphFromPolynomial(pVOptPubo, GpVsat)
        nx.write_gml(GpVPubo, EXP_GRAPH_DIRECTORY + Path(filename).stem + "_pVPubo.gml.bz2")

    if directSATviable:
        #Qubo from Pubo from non optimised k-SAT
        tmnpVQubo1 = time.time()
        nonpVOptQuboPBF, nonpVOptPenalty = fastPolyQuadratisation(nonpVOptPubo, 2, selection_quantile=1.0)
        nonpVOptQubo = PolyBinary(nonpVOptQuboPBF) + (PolyBinary(nonpVOptPenalty) * nonPVOptPenaltyScalingFactor)
        tmnpVQubo2 = time.time()
        nonpVOptQuboDensities = getDensities(nonpVOptQubo)
        DIres["npVQubo_deg"] = nonpVOptQubo.degree
        DIres["npVQubo_num_vars"] = len(nonpVOptQubo.variables)
        DIres["npVQubo_num_monomials"] = len(nonpVOptQubo.keys())
        DIres["npVQubo_density_deg_0"] = nonpVOptQuboDensities[0]
        DIres["npVQubo_density_deg_1"] = nonpVOptQuboDensities[1]
        DIres["npVQubo_density_deg_2"] = nonpVOptQuboDensities[2]
        DIres["npVQubo_time_fastPoly"] = tmnpVQubo2 - tmnpVQubo1
        if graph_eval:
            GnpVQubo = primalGraphFromPolynomial(nonpVOptQubo, GnpVPubo)
            nx.write_gml(GnpVQubo, EXP_GRAPH_DIRECTORY + Path(filename).stem + "_npVQubo.gml.bz2")

    #Qubo from Pubo from optimised k-SAT
    tmpVQubo1 = time.time()
    pVOptQuboPBF, pVOptPenalty = fastPolyQuadratisation(pVOptPubo, 2, selection_quantile=1.0)
    pVOptQubo = PolyBinary(pVOptQuboPBF) + (PolyBinary(pVOptPenalty) * pVOptPenaltyScalingFactor)
    tmpVQubo2 = time.time()
    pVOptQuboDensities = getDensities(pVOptQubo)
    DIres["pVQubo_deg"] = pVOptQubo.degree
    DIres["pVQubo_num_vars"] = len(pVOptQubo.variables)
    DIres["pVQubo_num_monomials"] = len(pVOptQubo.keys())
    DIres["pVQubo_density_deg_0"] = pVOptQuboDensities[0]
    DIres["pVQubo_density_deg_1"] = pVOptQuboDensities[1]
    DIres["pVQubo_density_deg_2"] = pVOptQuboDensities[2]
    DIres["pVQubo_time_fastPoly"] = tmpVQubo2 - tmpVQubo1

    if graph_eval:
        GpVQubo = primalGraphFromPolynomial(pVOptQubo, GpVPubo)
        nx.write_gml(GpVQubo, EXP_GRAPH_DIRECTORY + Path(filename).stem + "_pVQubo.gml.bz2")

    #Do landscape evaluation of pbf's
    if pbf_landscape:
        LSres = {}
        LSres["SAT_filename"] = filename
        LSres["LS_pbf_type"] = []
        LSres["LS_bitstr"] = []
        LSres["LS_energy"] = []
        if directSATviable:
            for vAStr, En in sorted(getPBFLandscape(nonpVOptPubo)):
                LSres["LS_pbf_type"] += ["npVPubo"]
                LSres["LS_bitstr"] += [vAStr]
                LSres["LS_energy"] += [En]
        
        for vAStr, En in sorted(getPBFLandscape(pVOptPubo)):
            LSres["LS_pbf_type"] += ["pVPubo"]
            LSres["LS_bitstr"] += [vAStr]
            LSres["LS_energy"] += [En]
        if directSATviable:
            for vAStr, En in sorted(getPBFLandscape(nonpVOptQubo)):
                LSres["LS_pbf_type"] += ["npVQubo"]
                LSres["LS_bitstr"] += [vAStr]
                LSres["LS_energy"] += [En]

        for vAStr, En in sorted(getPBFLandscape(pVOptQubo)):
            LSres["LS_pbf_type"] += ["pVQubo"]
            LSres["LS_bitstr"] += [vAStr]
            LSres["LS_energy"] += [En]

        csv_string = EXP_INDUSTRY_DIRECTORY + Path(filename).stem + "_LS.csv"
        df = pd.DataFrame(LSres)
        df.to_csv(csv_string)

    #Do simulated annealing of pbf's
    if pbf_simulatedAnnealing:
        SAres = {}
        SAres["SAT_filename"] = filename
        SAres["SA_pbf_type"] = []
        SAres["SA_steps"] = []
        SAres["SA_steps_variant"] = []
        SAres["SA_init_temp"] = []
        SAres["SA_seed"] = []
        SAres["SA_bitstr"] = []
        SAres["SA_energy_current_pbf"] = []
        SAres["SA_energy_original_pbf"] = []
        SAres["SA_time"] = []

        seed_start = 42
        loop_cnt = 0
        for i in range(100):
            for steps in [v, v*v]:
                initTempnonPv = nonPVOptPenaltyScalingFactor
                initTempPv = pVOptPenaltyScalingFactor
                
                while initTempnonPv > 1 or initTempPv > 1:
                    random.seed(loop_cnt+seed_start)
                    loop_cnt += 1

                    seed = random.randint(0, 2e31)
                    
                    if directSATviable and initTempnonPv > 1:
                        #NonpVOpt Pubo
                        tmnpVPuboSA = time.time()
                        vAStr, vA, En = simulatedAnnealing(nonpVOptPubo, steps, initTempnonPv, seed=seed)
                        SAres["SA_pbf_type"] += ["npVPubo"]
                        SAres["SA_steps"] += [steps]
                        SAres["SA_steps_variant"] += ["linear" if steps==v else "quadratic"]
                        SAres["SA_init_temp"] += [initTempnonPv]
                        SAres["SA_seed"] += [seed]
                        SAres["SA_bitstr"] += [vAStr]
                        SAres["SA_energy_current_pbf"] += [En]
                        SAres["SA_energy_original_pbf"] += [En]
                        SAres["SA_time"] += [time.time() - tmnpVPuboSA]


                    if directSATviable and initTempnonPv > 1:
                        #NonpVOpt Qubo
                        tmnpVQuboSA = time.time()
                        vAStr, vA, En = simulatedAnnealing(nonpVOptQubo, steps, initTempnonPv, seed=seed)
                        SAres["SA_pbf_type"] += ["npVQubo"]
                        SAres["SA_steps"] += [steps]
                        SAres["SA_steps_variant"] += ["linear" if steps==v else "quadratic"]
                        SAres["SA_init_temp"] += [initTempnonPv]
                        SAres["SA_seed"] += [seed]
                        SAres["SA_bitstr"] += [vAStr]
                        SAres["SA_energy_current_pbf"] += [En]
                        SAres["SA_energy_original_pbf"] += [nonpVOptPubo.evaluate(vA)]
                        SAres["SA_time"] += [time.time() - tmnpVQuboSA]
                        
                    if initTempPv > 1:
                        #pVOpt Pubo
                        tmpVPuboSA = time.time()
                        vAStr, vA, En = simulatedAnnealing(pVOptPubo, steps, initTempPv, seed=seed)
                        SAres["SA_pbf_type"] += ["pVPubo"]
                        SAres["SA_steps"] += [steps]
                        SAres["SA_steps_variant"] += ["linear" if steps==v else "quadratic"]
                        SAres["SA_init_temp"] += [initTempPv]
                        SAres["SA_seed"] += [seed]
                        SAres["SA_bitstr"] += [vAStr]
                        SAres["SA_energy_current_pbf"] += [En]
                        SAres["SA_energy_original_pbf"] += [En]
                        SAres["SA_time"] += [time.time() - tmpVPuboSA]

                        #pVOpt Pubo
                        tmpVQuboSA = time.time()
                        vAStr, vA, En = simulatedAnnealing(pVOptQubo, steps, initTempPv, seed=seed)
                        SAres["SA_pbf_type"] += ["pVQubo"]
                        SAres["SA_steps"] += [steps]
                        SAres["SA_steps_variant"] += ["linear" if steps==v else "quadratic"]
                        SAres["SA_init_temp"] += [initTempPv]
                        SAres["SA_seed"] += [seed]
                        SAres["SA_bitstr"] += [vAStr]
                        SAres["SA_energy_current_pbf"] += [En]
                        SAres["SA_energy_original_pbf"] += [pVOptPubo.evaluate(vA)]
                        SAres["SA_time"] += [time.time() - tmpVQuboSA]

                    initTempnonPv /= 2
                    initTempPv /= 2
        
        csv_string = EXP_INDUSTRY_DIRECTORY + Path(filename).stem + "_SA.csv"
        df = pd.DataFrame(SAres)
        df.to_csv(csv_string)

    if directSATviable:
        #Ising from Pubo from non optimised k-SAT
        tmnpVPuboIsing1 = time.time()
        nonpVOptPuboIsing = nonpVOptPubo.to_ising(inverted=True)
        tmnpVPuboIsing2 = time.time()
        nonpVOptPuboIsingDensities = getDensities(nonpVOptPuboIsing)
        DIres["npVPuboIsing_deg"] = nonpVOptPuboIsing.degree
        DIres["npVPuboIsing_num_vars"] = len(nonpVOptPuboIsing.variables)
        DIres["npVPuboIsing_num_monomials"] = len(nonpVOptPuboIsing.keys())
        for deg, val in nonpVOptPuboIsingDensities.items():
            DIres[("npVPuboIsing_density_deg_" + str(deg))] = val
        DIres["npVPuboIsing_time_toIsing"] = tmnpVPuboIsing2 - tmnpVPuboIsing1
        if graph_eval:
            GnpVPuboIsing = primalGraphFromPolynomial(nonpVOptPuboIsing, GnpVPubo)
            nx.write_gml(GnpVPuboIsing, EXP_GRAPH_DIRECTORY + Path(filename).stem + "_npVPuboIsing.gml.bz2")
    
    #Ising from Pubo from optimised k-SAT
    tmpVPuboIsing1 = time.time()
    pVOptPuboIsing = pVOptPubo.to_ising(inverted=True)
    tmpVPuboIsing2 = time.time()
    pVOptPuboIsingDensities = getDensities(pVOptPuboIsing)
    DIres["pVPuboIsing_deg"] = pVOptPuboIsing.degree
    DIres["pVPuboIsing_num_vars"] = len(pVOptPuboIsing.variables)
    DIres["pVPuboIsing_num_monomials"] = len(pVOptPuboIsing.keys())
    for deg, val in pVOptPuboIsingDensities.items():
        DIres[("npVPuboIsing_density_deg_" + str(deg))] = val
    DIres["pVPuboIsing_time_toIsing"] = tmpVPuboIsing2 - tmpVPuboIsing1
    if graph_eval:
        GpVPuboIsing = primalGraphFromPolynomial(pVOptPuboIsing, GpVPubo)
        nx.write_gml(GpVPuboIsing, EXP_GRAPH_DIRECTORY + Path(filename).stem + "_pVPuboIsing.gml.bz2")

    if directSATviable:
        #Ising from Qubo from Pubo from non optimised k-SAT
        tmnpVQuboIsing1 = time.time()
        nonpVOptQuboIsing = nonpVOptQubo.to_ising(inverted=True)
        tmnpVQuboIsing2 = time.time()
        nonpVOptQuboIsingDensities = getDensities(nonpVOptQuboIsing)
        DIres["npVQuboIsing_deg"] = nonpVOptQuboIsing.degree
        DIres["npVQuboIsing_num_vars"] = len(nonpVOptQuboIsing.variables)
        DIres["npVQuboIsing_num_monomials"] = len(nonpVOptQuboIsing.keys())
        DIres["npVQuboIsing_density_deg_0"] = nonpVOptQuboIsingDensities[0]
        DIres["npVQuboIsing_density_deg_1"] = nonpVOptQuboIsingDensities[1]
        DIres["npVQuboIsing_density_deg_2"] = nonpVOptQuboIsingDensities[2]
        DIres["npVQuboIsing_time_toIsing"] = tmnpVQuboIsing2 - tmnpVQuboIsing1
        if graph_eval:
            GnpVQuboIsing = primalGraphFromPolynomial(nonpVOptQuboIsing, GnpVQubo)
            nx.write_gml(GnpVQuboIsing, EXP_GRAPH_DIRECTORY + Path(filename).stem + "_npVQuboIsing.gml.bz2")

    #Ising from Qubo from Pubo from optimised k-SAT
    tmpVQuboIsing1 = time.time()
    pVOptQuboIsing = pVOptQubo.to_ising(inverted=True)
    tmpVQuboIsing2 = time.time()
    pVOptQuboIsingDensities = getDensities(pVOptQuboIsing)
    DIres["pVQuboIsing_deg"] = pVOptQuboIsing.degree
    DIres["pVQuboIsing_num_vars"] = len(pVOptQuboIsing.variables)
    DIres["pVQuboIsing_num_monomials"] = len(pVOptQuboIsing.keys())
    DIres["pVQuboIsing_density_deg_0"] = pVOptQuboIsingDensities[0]
    DIres["pVQuboIsing_density_deg_1"] = pVOptQuboIsingDensities[1]
    DIres["pVQuboIsing_density_deg_2"] = pVOptQuboIsingDensities[2]
    DIres["pVQuboIsing_time_toIsing"] = tmpVQuboIsing2 - tmpVQuboIsing1
    if graph_eval:
        GpVQuboIsing = primalGraphFromPolynomial(pVOptQuboIsing, GpVQubo)
        nx.write_gml(GpVQuboIsing, EXP_GRAPH_DIRECTORY + Path(filename).stem + "_pVQuboIsing.gml.bz2")

    if directSATviable:
        #Logical QAOA from Ising from Pubo from non optimised k-SAT
        tmnpVPuboQAOA1 = time.time()
        nonpVOptPuboIsingLQaoa = create_LR_QAOA_Circuit(nonpVOptPuboIsing, 0.5, 0.5, 1)
        tmnpVPuboQAOA2 = time.time()
        DIres["npVPuboIsingLQaoa_num_qubits"] = nonpVOptPuboIsingLQaoa.num_qubits
        DIres["npVPuboIsingLQaoa_depth"] = nonpVOptPuboIsingLQaoa.depth()
        DIres["npVPuboIsingLQaoa_num_non_local_gates"] = nonpVOptPuboIsingLQaoa.num_nonlocal_gates()
        nonpVOptPuboIsingLQaoaOps = nonpVOptPuboIsingLQaoa.count_ops()
        if "rz" in nonpVOptPuboIsingLQaoaOps:
            DIres["npVPuboIsingLQaoa_num_local_gates"] = nonpVOptPuboIsingLQaoaOps["rz"] + nonpVOptPuboIsingLQaoaOps["rx"]
        else:
            DIres["npVPuboIsingLQaoa_num_local_gates"] = nonpVOptPuboIsingLQaoaOps["rx"]
        DIres["npVPuboIsingLQaoa_time_genCirc"] = tmnpVPuboQAOA2 - tmnpVPuboQAOA1

    #Logical QAOA from Ising from Pubo from optimised k-SAT
    tmpVPuboQAOA1 = time.time()
    pVOptPuboIsingLQaoa = create_LR_QAOA_Circuit(pVOptPuboIsing, 0.5, 0.5, 1)
    tmpVPuboQAOA2 = time.time()
    DIres["pVPuboIsingLQaoa_num_qubits"] = pVOptPuboIsingLQaoa.num_qubits
    DIres["pVPuboIsingLQaoa_depth"] = pVOptPuboIsingLQaoa.depth()
    DIres["pVPuboIsingLQaoa_num_non_local_gates"] = pVOptPuboIsingLQaoa.num_nonlocal_gates()
    pVOptPuboIsingLQaoaOps = pVOptPuboIsingLQaoa.count_ops()
    if "rz" in pVOptPuboIsingLQaoaOps:
        DIres["pVPuboIsingLQaoa_num_local_gates"] = pVOptPuboIsingLQaoaOps["rz"] + pVOptPuboIsingLQaoaOps["rx"]
    else:
        DIres["pVPuboIsingLQaoa_num_local_gates"] = pVOptPuboIsingLQaoaOps["rx"]
    DIres["pVPuboIsingLQaoa_time_genCirc"] = tmpVPuboQAOA2 - tmpVPuboQAOA1

    if directSATviable:
        #Logical QAOA Ising from Qubo from Pubo from non optimised k-SAT
        tmnpVQuboQAOA1 = time.time()
        nonpVOptQuboIsingLQaoa = create_LR_QAOA_Circuit(nonpVOptQuboIsing, 0.5, 0.5, 1)
        tmnpVQuboQAOA2 = time.time()
        DIres["npVQuboIsingLQaoa_num_qubits"] = nonpVOptQuboIsingLQaoa.num_qubits
        DIres["npVQuboIsingLQaoa_depth"] = nonpVOptQuboIsingLQaoa.depth()
        DIres["npVQuboIsingLQaoa_num_non_local_gates"] = nonpVOptQuboIsingLQaoa.num_nonlocal_gates()
        nonpVOptQuboIsingLQaoaOps = nonpVOptQuboIsingLQaoa.count_ops()
        if "rz" in nonpVOptQuboIsingLQaoaOps:
            DIres["npVQuboIsingLQaoa_num_local_gates"] = nonpVOptQuboIsingLQaoaOps["rz"] + nonpVOptQuboIsingLQaoaOps["rx"]
        else:
            DIres["npVQuboIsingLQaoa_num_local_gates"] = nonpVOptQuboIsingLQaoaOps["rx"]
        DIres["npVQuboIsingLQaoa_time_genCirc"] = tmnpVQuboQAOA2 - tmnpVQuboQAOA1

    #Logical QAOA Ising from Qubo from Pubo from optimised k-SAT
    tmpVQuboQAOA1 = time.time()
    pVOptQuboIsingLQaoa = create_LR_QAOA_Circuit(pVOptQuboIsing, 0.5, 0.5, 1)
    tmpVQuboQAOA2 = time.time()
    DIres["pVQuboIsingLQaoa_num_qubits"] = pVOptQuboIsingLQaoa.num_qubits
    DIres["pVQuboIsingLQaoa_depth"] = pVOptQuboIsingLQaoa.depth()
    DIres["pVQuboIsingLQaoa_num_non_local_gates"] = pVOptQuboIsingLQaoa.num_nonlocal_gates()
    pVOptQuboIsingLQaoaOps = pVOptQuboIsingLQaoa.count_ops()
    if "rz" in pVOptQuboIsingLQaoaOps:
        DIres["pVQuboIsingLQaoa_num_local_gates"] = pVOptQuboIsingLQaoaOps["rz"] + pVOptQuboIsingLQaoaOps["rx"]
    else:
        DIres["pVQuboIsingLQaoa_num_local_gates"] = pVOptQuboIsingLQaoaOps["rx"]
    DIres["pVQuboIsingLQaoa_time_genCirc"] = tmpVQuboQAOA2 - tmpVQuboQAOA1

    csv_string = EXP_INDUSTRY_DIRECTORY + Path(filename).stem + ".csv"
    df = pd.DataFrame(DIres, index=[0])
    df.to_csv(csv_string)

    #Do noiseless simulation with qaoa
    if qaoa_simulation:        
        shots = 10000
        dbeta = 0.5
        dgamma = 0.5
        
        res = pd.DataFrame()
        for p in [100, 1000]:
            #Logical QAOA from Ising from Pubo from non optimised k-SAT
            if directSATviable:
                nonpVOptPuboIsingLQaoaDf = doQAOASimulation(ising=nonpVOptPuboIsing, pubo=nonpVOptPubo, qubo=None, dbeta=dbeta, dgamma=dgamma, p=p, shots=shots, pbf_type="npVPubo", filename=filename)
                pVOptPuboIsingLQaoaDf = doQAOASimulation(ising=pVOptPuboIsing, pubo=pVOptPubo, qubo=None, dbeta=dbeta, dgamma=dgamma, p=p, shots=shots, pbf_type="pVPubo", filename=filename)
                nonpVOptQuboIsingLQaoaDf = doQAOASimulation(ising=nonpVOptQuboIsing, pubo=nonpVOptPubo, qubo=nonpVOptQubo, dbeta=dbeta, dgamma=dgamma, p=p, shots=shots, pbf_type="npVQubo", filename=filename)
                pVOptQuboIsingLQaoaDf = doQAOASimulation(ising=pVOptQuboIsing, pubo=pVOptPubo, qubo=pVOptQubo, dbeta=dbeta, dgamma=dgamma, p=p, shots=shots, pbf_type="pVQubo", filename=filename)
                res = pd.concat([nonpVOptPuboIsingLQaoaDf, pVOptPuboIsingLQaoaDf, nonpVOptQuboIsingLQaoaDf, pVOptQuboIsingLQaoaDf, res], ignore_index=True)
            else:
                pVOptPuboIsingLQaoaDf = doQAOASimulation(ising=pVOptPuboIsing, pubo=pVOptPubo, qubo=None, dbeta=dbeta, dgamma=dgamma, p=p, shots=shots, pbf_type="pVPubo", filename=filename)
                pVOptQuboIsingLQaoaDf = doQAOASimulation(ising=pVOptQuboIsing, pubo=pVOptPubo, qubo=pVOptQubo, dbeta=dbeta, dgamma=dgamma, p=p, shots=shots, pbf_type="pVQubo", filename=filename)
                res = pd.concat([pVOptPuboIsingLQaoaDf, pVOptQuboIsingLQaoaDf, res], ignore_index=True)
           

        csv_string = EXP_INDUSTRY_DIRECTORY + Path(filename).stem + "_QAOA.csv"
        res.to_csv(csv_string)

    return DIres

def testPerformance():
    """
    Compares the time for a reduction for the graph based reduction vs the standard monomial based reduction
    """
    
    shutil.rmtree("Experiments/")
    Path("Experiments/").mkdir(parents=True, exist_ok=True)

    #Graph based
    selection_type = []
    test_polys = []
    for i in range(6, 90, 2):
        for d in range(1, 11, 1):
            for st in range(5, 11, 1):
                test_polys.append((i, 4, d / 10.0))
                selection_type.append(st / 10.0)

    GBtest_polys2 = test_polys[::-1]
    GBselection_type2 = selection_type[::-1]

    selection_type = []
    test_polys = []
    for i in range(6, 90, 2):
        for d in range(1, 11, 1):
            for st in range(8, 11, 1):
                if i > 39:
                    if st == 8:
                        test_polys.append((i, 4, d / 10.0))
                        selection_type.append(st / 10.0)
                else:
                    test_polys.append((i, 4, d / 10.0))
                    selection_type.append(st / 10.0)

    MBtest_polys2 = test_polys[::-1]
    MBselection_type2 = selection_type[::-1]
    
    completeExpGB = p_imap(testPerformanceSingleGB, GBtest_polys2, GBselection_type2)

    completeExpMB = p_imap(testPerformanceSingleMB, MBtest_polys2, MBselection_type2)

    pdres = pd.concat(
        [pd.DataFrame(completeExpMB), pd.DataFrame(completeExpGB)],
        ignore_index=True,
        axis=0,
    )
    pdres.to_csv("Performance_results.csv")

def testPerformanceSingleMB(poly_config, selection_type):
    """
    Compares the time for a reduction for the standard monomial based reduction
    Input:
    poly_config: (<NumVars>, <MaxDegree>, <Density>) for createPolyFunction
    """
    v, deg, dens = poly_config
    poly_as_dict_flat = createPoly(v, deg, density=dens)
    origPolyDensities = getDensities(poly_as_dict_flat)

    poly = PolyBinary(poly_as_dict_flat)
    num_terms = len(poly.keys())
    num_vars_before = len(poly.variables)

    # Monomial based
    choice = "better"
    if selection_type <= 0.8:
        choice = "simple"
    elif selection_type <= 0.9:
        choice = "medium"

    tm1 = time.time()
    qPoly, penalty = poly.reduce(
        max_degree=2, var_pair_choice=choice, reduction_variable_prefix="r"
    )
    tm2 = time.time()

    MBresPoly = PolyBinary(qPoly) + PenaltyToPolyBinary(penalty)
    MBresDensities = getDensities(MBresPoly)

    monomial_all_vars = set()
    for m in qPoly.keys():
        for var in m:
            monomial_all_vars.add(var)

    # Monomial based
    MBresults = {}
    MBresults["num_variables_before"] = num_vars_before
    MBresults["poly_len"] = num_terms
    MBresults["time"] = tm2 - tm1
    MBresults["monomial_or_graph_based"] = "monomial_based"
    MBresults["Selection_type"] = choice
    MBresults["gen_type"] = str((v, deg, dens))
    MBresults["num_variables_after"] = len(monomial_all_vars)
    MBresults["density_deg1"] = origPolyDensities[1]
    MBresults["density_deg2"] = origPolyDensities[2]
    MBresults["density_deg3"] = origPolyDensities[3]
    MBresults["density_deg4"] = origPolyDensities[4]
    MBresults["density_deg1r"] = MBresDensities[1]
    MBresults["density_deg2r"] = MBresDensities[2]
    MBresults["density_deg3r"] = MBresDensities[3]
    MBresults["density_deg4r"] = MBresDensities[4]

    csv_string = "Experiments/MB" + str(poly_config) + str(selection_type) + ".csv"
    df = pd.DataFrame(MBresults, index=[0])
    df.to_csv(csv_string)

    return MBresults

def testPerformanceSingleGB(poly_config, selection_type):
    """
    Compares the time for a reduction for the graph based reduction
    Input:
    poly_config: (<NumVars>, <MaxDegree>, <Density>) for createPolyFunction
    """
    v, deg, dens = poly_config
    poly_as_dict_flat = createPoly(v, deg, density=dens)
    origPolyDensities = getDensities(poly_as_dict_flat)

    poly = PolyBinary(poly_as_dict_flat)
    num_terms = len(poly.keys())
    num_vars_before = len(poly.variables)

    # Graph based
    tg1 = time.time()
    qPoly, penalty = fastPolyQuadratisation(poly, 2, selection_quantile=selection_type)
    tg2 = time.time()

    GBresPoly = PolyBinary(qPoly) + PolyBinary(penalty)
    GBresDensities = getDensities(GBresPoly)

    graph_all_vars = set()
    for m in qPoly.keys():
        for var in m:
            graph_all_vars.add(var)

    for m in penalty.keys():
        if len(m) == 1:
            graph_all_vars.add(m[0])

    # Graph based
    GBresults = {}
    GBresults["num_variables_before"] = num_vars_before
    GBresults["poly_len"] = num_terms
    GBresults["time"] = tg2 - tg1
    GBresults["monomial_or_graph_based"] = "graph_based"
    GBresults["Selection_type"] = selection_type
    GBresults["gen_type"] = str((v, deg, dens))
    GBresults["num_variables_after"] = len(graph_all_vars)
    GBresults["density_deg1"] = origPolyDensities[1]
    GBresults["density_deg2"] = origPolyDensities[2]
    GBresults["density_deg3"] = origPolyDensities[3]
    GBresults["density_deg4"] = origPolyDensities[4]
    GBresults["density_deg1r"] = GBresDensities[1]
    GBresults["density_deg2r"] = GBresDensities[2]
    GBresults["density_deg3r"] = GBresDensities[3]
    GBresults["density_deg4r"] = GBresDensities[4]

    csv_string = "Experiments/GB" + str(poly_config) + str(selection_type) + ".csv"
    df = pd.DataFrame(GBresults, index=[0])
    df.to_csv(csv_string)

    return GBresults

def createPoly(variables: int, degree: int, density: float = 1.0, seed=42):
    """
    Creates polynomial dicts of the form {<monomial>: <alpha>} (e.g. {(1,2,3,4,5): 7.3})
    variables: number of variables in poly
    degree: polynomial's degree
    density: density for degree-k polynomials
    """
    random.seed(seed)
    out = dict()
    if variables < degree:
        raise ValueError("createPoly: degree must be at most #variables")
    varList = [x for x in range(0, variables)]

    for i in range(1, degree + 1):
        monomials = list(combinations(varList, i))
        for m in monomials:
            if random.random() < density:
                out[tuple(m)] = i + 0.3

    return out

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

def PenaltyToPolyBinary(penalty):
    """
    Input: Penalty terms of the form [(<x_i>, <x_j>, <y_h>), ...]
    Returns: PolyBinary
    """
    ret = PolyBinary()
    for xi, xj, yh in penalty:
        ret = ret + PolyBinary({(yh,): 3, (xi, xj): 1, (xi, yh): -2, (xj, yh): -2})
    return ret
