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
from FastReduc import *


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
        # create degree-i monomials:
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
    # Create a dictionary of the form:
    # {0: [(<alpha>, (<deg0-monomial>)), ...],
    #  1: [(<alpha>, (<deg1-monomial>)), ...],
    # }
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
