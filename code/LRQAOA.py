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

from qiskit import QuantumCircuit
import numpy as np
import networkx as nx
import warnings
import math
from edge_coloring import *
from quark import PolyIsing

from util import *

def create_LR_QAOA_Circuit(poly: PolyIsing, dbeta, dgamma, p):
    """
    The optimal bitstring will be flipped (::-1) (e.g. if x = 10010 is optimal in the IsingModel, s = 01001 is the optimal solution for the QuantumCircuit)
    Creates a depth optimised LR-QAOA circuit with p layers based on the provided IsingModel. It automatically decomposes R_Z^n, n>1 gates into CNOT ladders and R_Z gates. 
    This circuit needs to undergo a hardware transpile step to solve the qubit assignment and qubit routing problem.
    Input:
        poly: a flat (numerical indices) quark PolyIsing 
        dbeta: float
        dgamma: float
        p: int
    """
    normalizationFactor = 0
    for monomial, alpha in poly.items():
        if len(monomial) > 0:
            normalizationFactor = max([normalizationFactor, abs(alpha)])
    zero_idx_offset = min(poly.variables)
    if len(poly.variables) < max(poly.variables) - min(poly.variables) + 1:
        warnings.warn("The input polynomial has missing consecutive variable indices! This will lead to unused qubits in the quantum circuit.")
    qc = QuantumCircuit(max(poly.variables) - min(poly.variables) + 1)
    for i in poly.variables:
        qc.h(i - zero_idx_offset)
    for layer in range(p):
        betai = (1 - (layer / p)) * dbeta
        gammai = ((1 + layer) / p) * dgamma
        
        for monomial, val in poly.items():
            if len(monomial) > 0:
                for i in range(1,len(monomial)):
                    qc.cx(control_qubit=monomial[i-1]-zero_idx_offset, target_qubit=monomial[i]-zero_idx_offset)
                qc.rz(qubit=monomial[-1]-zero_idx_offset, phi=(val*gammai*math.pi)/normalizationFactor)
                for i in range(len(monomial)-1,0, -1):
                    qc.cx(control_qubit=monomial[i-1]-zero_idx_offset, target_qubit=monomial[i]-zero_idx_offset)

        for r in poly.variables:
            qc.rx(theta=-2 * betai, qubit=r-zero_idx_offset)

    return qc
