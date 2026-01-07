#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 09:23:22 2022

@author: mafal
"""
from openfermion import MolecularData
from openfermionpyscf import run_pyscf


def create_h2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        h2 (PyscfMolecularData): the linear H2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['H', [0, 0, 0]], ['H', [0, 0, r]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    h2 = MolecularData(geometry, basis, multiplicity, charge, description='H2')
    h2 = run_pyscf(h2, run_fci=True, run_ccsd=True)

    return h2


def create_h3(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        h3 (PyscfMolecularData): the linear H3 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['H', [0, 0, 0]], ['H', [0, 0, r]], ['H', [0, 0, 2 * r]]]
    basis = 'sto-3g'
    multiplicity = 2  # odd number of electrons
    charge = 0
    h3 = MolecularData(geometry, basis, multiplicity, charge, description='H3')
    h3 = run_pyscf(h3, run_fci=True, run_ccsd=False)  # CCSD doesn't work here?

    return h3


def create_h4(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        h4 (PyscfMolecularData): the linear H4 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [('H', (0, 0, 0)), ('H', (0, 0, r)), ('H', (0, 0, 2 * r)),
                ('H', (0, 0, 3 * r))]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    h4 = MolecularData(geometry, basis, multiplicity, charge, description='H4')
    h4 = run_pyscf(h4, run_fci=True, run_ccsd=True)

    return h4


def create_h5(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        h5 (PyscfMolecularData): the linear H5 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [('H', (0, 0, 0)), ('H', (0, 0, r)), ('H', (0, 0, 2 * r)),
                ('H', (0, 0, 3 * r)), ('H', (0, 0, 4 * r))]
    basis = 'sto-3g'
    multiplicity = 2  # odd number of electrons
    charge = 0
    h5 = MolecularData(geometry, basis, multiplicity, charge, description='H5')
    h5 = run_pyscf(h5, run_fci=True, run_ccsd=False)  # CCSD doesn't work here?

    return h5


def create_h6(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        h6 (PyscfMolecularData): the linear H6 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [('H', (0, 0, 0)), ('H', (0, 0, r)), ('H', (0, 0, 2 * r)),
                ('H', (0, 0, 3 * r)), ('H', (0, 0, 4 * r)), ('H', (0, 0, 5 * r))]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    h6 = MolecularData(geometry, basis, multiplicity, charge, description='H6')
    h6 = run_pyscf(h6, run_fci=True, run_ccsd=True)

    return h6


def create_h7(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        h7 (PyscfMolecularData): the linear H7 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [('H', (0, 0, 0)), ('H', (0, 0, r)), ('H', (0, 0, 2 * r)),
                ('H', (0, 0, 3 * r)), ('H', (0, 0, 4 * r)), ('H', (0, 0, 5 * r)), ('H', (0, 0, 6 * r))]
    basis = 'sto-3g'
    multiplicity = 2  # odd number of electrons
    charge = 0
    h7 = MolecularData(geometry, basis, multiplicity, charge, description='H7')
    h7 = run_pyscf(h7, run_fci=True, run_ccsd=False)

    return h7


def create_lih(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        lih (PyscfMolecularData): the LiH molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Li', [0, 0, 0]], ['H', [0, 0, r]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    li_h = MolecularData(geometry, basis, multiplicity, charge, description='LiH')
    li_h = run_pyscf(li_h, run_fci=True, run_ccsd=True)

    return li_h


def create_beh2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        beh2 (PyscfMolecularData): the BeH2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['H', [0, 0, 0]], ['Be', [0, 0, r]], ['H', [0, 0, 2 * r]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    be_h2 = MolecularData(geometry, basis, multiplicity, charge, description='BeH2')
    be_h2 = run_pyscf(be_h2, run_fci=True, run_ccsd=True)

    return be_h2

def create_ch4():
    """
    Returns:
        ch4 (PyscfMolecularData): the CH4 molecule in the minimal STO-3G basis set
    """
    geometry = [
        ['C', [0.0, 0.0, 0.0]],
        ['H', [0.629, 0.629, 0.629]],
        ['H', [-0.629, -0.629, 0.629]],
        ['H', [0.629, -0.629, -0.629]],
        ['H', [-0.629, 0.629, -0.629]]
    ]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    ch4 = MolecularData(geometry, basis, multiplicity, charge, description='CH4')
    ch4 = run_pyscf(ch4, run_fci=True, run_ccsd=True)

    return ch4

def create_h2o(r):

    """

    Arguments:

        r (float): interatomic distance (angstrom)

    Returns:

        h2o (PyscfMolecularData): the H2O molecule at interatomic distance r, in the minimal STO-3G basis set

    """

    geometry = [['O', [2.5369 * r, -0.1550 * r, 0]], ['H', [3.0739 * r, 0.1550 * r, 0]], ['H', [2.0000 * r, 0.1550 * r, 0]]]

    basis = 'sto-3g'

    multiplicity = 1

    charge = 0

    h2o = MolecularData(geometry, basis, multiplicity, charge, description='H20')

    h2o = run_pyscf(h2o, run_fci=True, run_ccsd=True)

    return h2o

def create_mgh2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        mgh2 (PyscfMolecularData): the MgH2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Mg', [2.5369 * r, -0.1550 * r, 0]], ['H', [3.0739 * r, 0.1550 * r, 0]], ['H', [2.0000 * r, 0.1550 * r, 0]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    mgh2 = MolecularData(geometry, basis, multiplicity, charge, description='MgH2')
    mgh2 = run_pyscf(mgh2, run_fci=True, run_ccsd=True)

    return mgh2

def create_sio2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        sio2 (PyscfMolecularData): the SiO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Si', [0, -0.5548 * r, 0]], ['O', [1.2172 * r, 0.2774 * r, 0]], ['O', [-1.2171 * r, 0.2774 * r, 0]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    sio2 = MolecularData(geometry, basis, multiplicity, charge, description='SiO2')
    sio2 = run_pyscf(sio2, run_fci=True, run_ccsd=True)

    return sio2

def create_mgo(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        mgo (PyscfMolecularData): the linear MgO molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['O', [2*r, 0, 0]], ['Mg', [3*r, 0 , 0]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    mgo = MolecularData(geometry, basis, multiplicity, charge, description='MgO')
    mgo = run_pyscf(mgo, run_fci=True, run_ccsd=True)

    return mgo

def create_co2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['C', [0, 0, 0]], ['O', [-1.1970 * r, 0 , 0]], ['O', [1.1970 * r, 0, 0]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    co2 = MolecularData(geometry, basis, multiplicity, charge, description='CO2')
    co2 = run_pyscf(co2, run_fci=True, run_ccsd=True)

    return co2

def create_tin(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        tin (PyscfMolecularData): the linear TiN molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Ti', [2*r, 0, 0]], ['N', [3*r, 0 , 0]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    tin = MolecularData(geometry, basis, multiplicity, charge, description='TiN')
    tin = run_pyscf(tin, run_fci=True, run_ccsd=True)

    return tin

def create_sic(r):

    """

    Arguments:

        r (float): here factor multiplied to real physical distance

    Returns:

        sic (PyscfMolecularData): the SiC molecule at interatomic distance r, in the minimal STO-3G basis set

    """

    geometry = [['Si', [0.7570, 0, 0]], ['C', [-0.7570, 0, 0]]]

    basis = 'sto-3g'

    multiplicity = 1

    charge = 0

    sic = MolecularData(geometry, basis, multiplicity, charge, description='SiC')

    sic = run_pyscf(sic, run_fci=True, run_ccsd=True)

    return sic

def create_sicl4(r):
    """
    Arguments:
        r (float): here factor multiplied to real physical distance
    Returns:
        sicl4 (PyscfMolecularData): the SiCl4 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Si', [-0.0002, 0.0004, 0.0001]], ['Cl', [-1.8192, 0.8785, 0.1915]], ['Cl', [1.4504, 1.4142, 0.1169]], ['Cl', [0.2563, -1.3559, 1.4872]], ['Cl', [0.1127, -0.9372, -1.7957]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    sicl4 = MolecularData(geometry, basis, multiplicity, charge, description='SiCl4')
    sicl4 = run_pyscf(sicl4, run_fci=True, run_ccsd=True)

    return sicl4

def create_ch2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['C', [2.5369, -0.1550, 0.0000]], ['H', [3.0739*r, 0.1550, 0.0000]], ['H', [2.0000*r,    0.1550,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    ch2 = MolecularData(geometry, basis, multiplicity, charge, description='CH2')
    ch2 = run_pyscf(ch2, run_fci=True, run_ccsd=True)

    return ch2

def create_hf(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['F', [3.0000,    0.0000,    0.0000]], ['H', [2.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    hf = MolecularData(geometry, basis, multiplicity, charge, description='HF')
    hf = run_pyscf(hf, run_fci=True, run_ccsd=True)

    return hf

def create_hn(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['N', [3.0000,    0.0000,    0.0000]], ['H', [2.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    hn = MolecularData(geometry, basis, multiplicity, charge, description='HN')
    hn = run_pyscf(hn, run_fci=True, run_ccsd=True)

    return hn

def create_hb(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['B', [3.0000,    0.0000,    0.0000]], ['H', [2.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    hb = MolecularData(geometry, basis, multiplicity, charge, description='HB')
    hb = run_pyscf(hb, run_fci=True, run_ccsd=True)

    return hb

def create_b2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['B', [2.0000,    0.0000,    0.0000]], ['B', [3.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    b2 = MolecularData(geometry, basis, multiplicity, charge, description='B2')
    b2 = run_pyscf(b2, run_fci=True, run_ccsd=True)

    return b2

def create_b3(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['B', [2.5000,   -0.4330,    0.0000]], ['B', [2.0000*r,    0.4330,    0.0000]], ['B', [3.0000*r,    0.4330,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    b3 = MolecularData(geometry, basis, multiplicity, charge, description='B3')
    b3 = run_pyscf(b3, run_fci=True, run_ccsd=True)

    return b3

def create_be2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Be', [2.0000,    0.0000,    0.0000]], ['Be', [3.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    be2 = MolecularData(geometry, basis, multiplicity, charge, description='Be2')
    be2 = run_pyscf(be2, run_fci=True, run_ccsd=True)

    return be2

def create_be3(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Be', [2.8660,   -0.2500,    0.0000]], ['Be', [3.7320*r,    0.2500,    0.0000]], ['Be', [2.0000*r,    0.2500,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    be3 = MolecularData(geometry, basis, multiplicity, charge, description='Be3')
    be3 = run_pyscf(be3, run_fci=True, run_ccsd=True)

    return be3

def create_be4(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Be', [2.8660,   -0.2500,    0.0000]], ['Be', [3.7320*r,    0.2500,    0.00000]], ['Be', [2.0000*r,    0.2500,    0.0000]], ['Be', [4.5981*r,   -0.2500,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    be4 = MolecularData(geometry, basis, multiplicity, charge, description='Be4')
    be4 = run_pyscf(be4, run_fci=False, run_ccsd=True)

    return create_be4

def create_c2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['C', [3.0000,    0.0000,    0.0000]], ['C', [2.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    c2 = MolecularData(geometry, basis, multiplicity, charge, description='C2')
    c2 = run_pyscf(c2, run_fci=True, run_ccsd=True)

    return c2

def create_c3(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['C', [2.8660,    0.5000,    0.0000]], ['C', [2.8660*r,   -0.5000,    0.0000]], ['C', [2.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    c3 = MolecularData(geometry, basis, multiplicity, charge, description='C3')
    c3 = run_pyscf(c3, run_fci=True, run_ccsd=True)

    return c3

def create_f2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['F', [2.0000,    0.0000,    0.0000]], ['F', [3.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    f2 = MolecularData(geometry, basis, multiplicity, charge, description='F2')
    f2 = run_pyscf(f2, run_fci=True, run_ccsd=True)

    return f2

def create_li2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Li', [0.0000,    0.0000,    0.0000]], ['Li', [1.0000*r,    0.0000,    0.0000 ]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    li2 = MolecularData(geometry, basis, multiplicity, charge, description='Li2')
    li2 = run_pyscf(li2, run_fci=True, run_ccsd=True)

    return li2

def create_li3(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Li', [0.0000,    0.0000,    0.0000]], ['Li', [1.0000*r,    0.0000,    0.0000]], ['Li', [2.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 2
    charge = 0
    li3 = MolecularData(geometry, basis, multiplicity, charge, description='Li3')
    li3 = run_pyscf(li3, run_fci=True, run_ccsd=True)

    return li3

def create_li4(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Li', [0.0000,    0.0000,    0.0000]], ['Li', [1.0000*r,    0.0000,    0.0000]], ['Li', [2.0000*r,    0.0000,    0.0000]], ['Li', [3.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    li4 = MolecularData(geometry, basis, multiplicity, charge, description='Li4')
    li4 = run_pyscf(li4, run_fci=True, run_ccsd=True)

    return li4

def create_li6(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['Li', [0.0000,    0.0000,    0.0000]], ['Li', [1.0000*r,    0.0000,    0.0000]], ['Li', [2.0000*r,    0.0000,    0.0000]], ['Li', [3.0000*r,    0.0000,    0.0000]], ['Li', [4.000*r,    0.0000,    0.0000]], ['Li', [5.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    li6 = MolecularData(geometry, basis, multiplicity, charge, description='Li6')
    li6 = run_pyscf(li6, run_fci=False, run_ccsd=True)

    return li6

def create_n2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['N', [2.0000,    0.0000,    0.0000]], ['N', [3.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    n2 = MolecularData(geometry, basis, multiplicity, charge, description='N2')
    n2 = run_pyscf(n2, run_fci=False, run_ccsd=True)

    return n2

def create_n3(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['N', [2.8660*r,    0.0000,    0.0000]], ['N', [3.7321*r,    0.5000,    0.0000]], ['N', [2.0000*r,   -0.5000,    0.00000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    n3 = MolecularData(geometry, basis, multiplicity, charge, description='N3')
    n3 = run_pyscf(n3, run_fci=True, run_ccsd=True)

    return n3

def create_o2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['O', [2.0000,    0.0000,    0.0000]], ['O', [3.0000*r,    0.0000,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    o2 = MolecularData(geometry, basis, multiplicity, charge, description='O2')
    o2 = run_pyscf(o2, run_fci=True, run_ccsd=True)

    return o2

def create_o3(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['O', [2.8660,   -0.2500,    0.0000]], ['O', [2.0000*r,    0.2500,    0.0000]], ['O', [3.7321*r,    0.2500,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    o3 = MolecularData(geometry, basis, multiplicity, charge, description='O3')
    o3 = run_pyscf(o3, run_fci=True, run_ccsd=True)

    return o3

def create_nh3(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['N', [2.5369,    0.1550,    0.0000]], ['H', [3.0739*r,    0.4650,    0.0000]], ['H', [2.0000*r,    0.4650,    0.0000]], ['H', [2.5369*r,   -0.4650,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    nh3 = MolecularData(geometry, basis, multiplicity, charge, description='NH3')
    nh3 = run_pyscf(nh3, run_fci=True, run_ccsd=True)

    return nh3

def create_bh2(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['B', [2.5369,   -0.1550,    0.0000]], ['H', [3.0739*r,    0.1550,    0.0000]], ['H', [2.0000*r,    0.1550,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 2
    charge = 0
    bh2 = MolecularData(geometry, basis, multiplicity, charge, description='BH2')
    bh2 = run_pyscf(bh2, run_fci=True, run_ccsd=True)

    return bh2

def create_bh3(r):
    """
    Arguments:
        r (float): interatomic distance (angstrom)
    Returns:
        co2 (PyscfMolecularData): the linear CO2 molecule at interatomic distance r, in the minimal STO-3G basis set
    """

    geometry = [['B', [2.5369,    0.1550,    0.0000]], ['H', [3.0739*r,    0.4650,    0.0000]], ['H', [2.0000*r,    0.4650,    0.0000]], ['H', [2.5369*r,   -0.4650,    0.0000]]]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    bh3 = MolecularData(geometry, basis, multiplicity, charge, description='BH3')
    bh3 = run_pyscf(bh3, run_fci=True, run_ccsd=True)

    return bh3

