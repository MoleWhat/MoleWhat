# MoleWhat database initiator
# Copyright (C) 2022  Karime Ochoa Jacinto
#                     Luis Aaron Nieto Cruz
#                     Miriam Guadalupe Valdez Maldonado
#                     Mariela Yael Arias Rojo
#                     Anton Pashkov
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""MoleWhat database generator
A script that initiates the database required for the MoleWhat program. This
database a single table that stores molecule and challenge information.
"""
import sqlite3
from pathlib import Path

# Create database if it doesn't already exists
if Path("Molecules.db").is_file() is False:

    # Connection to database
    con = sqlite3.connect("Molecules.db")
    cur = con.cursor()

    # Molecule table creation command
    TABLE_INIT = (
        "CREATE TABLE IF NOT EXISTS molecule ("
        "id INTERGER, "
        "name_iupac TEXT, "
        "description TEXT, "
        "code_smiles TEXT, "
        "completed INTEGER"
        ");"
    )

    # Creation of molecule table
    cur.execute(TABLE_INIT)

    # Data insertion
    data = [
        (
            1,
            "isopentane;methylbutane;2-methylbutane",
            (
                "The main chain of this alkane has four carbons (butane) with "
                "a methyl group linked to the second carbon, so the correct "
                "name would be 2-methylbutane, or simply methybutane. "
                "Furthermore, any alkane with a methyl group linked to the "
                "second carbon can be named by appending an 'iso' prefix; in "
                "this case, it would isopentane (as it has five carbons)."
            ),
            "CCC(C)C",
            0,
        ),
        (
            2,
            "1,2-dibromobenzene;o-dibromobenzene",
            (
                "Two bromine atoms are joined to a benzene ring on positions "
                "1 and 2, creating the 1,2-dibromobenzene. Because both "
                "bromine atoms are connected to adjacent carbons, this "
                "molecule can also be name o-dibromobenzene."
            ),
            "BrC1CCCCC1Br",
            0,
        ),
        (
            3,
            "ethylmethylpropylamine",
            (
                "Three alkyls (ethyl, methyl and propyl) are linked to a "
                "single nitrogen atom, forming an amine. IUPAC states that "
                "the ordering of the alkyls in the name must be alphabetical, "
                "so this molecule is called ethylmethylpropylamine."
            ),
            "CCCN(C)CC",
            0,
        ),
        (
            4,
            "2-phenyl-1-ethanol;2-phenylethanol;phenethyl alcohol;"
            "2-phenylethan-1-ol",
            (
                "This molecule is an alcohol. The numbering of the carbon "
                "atoms on the main chain (ethane) must begin from the carbon "
                "that is closest to the OH group. Thus, this molecule can be "
                "named in any way of the following: 2-phenyl-1-ethanol, 2-"
                "phenylethanol, phenethyl alcohol, or 2-phenylethan-1-ol."
            ),
            "C1:C:C:C:C:C1CCO",
            0,
        ),
        (
            5,
            "2,4,6-trinitrophenol;2,4,6-trinitro-1-phenol;2-hydroxy-1,3,5-"
            "trinitrobenzene",
            (
                "The carbon atoms on a phenol are numbered starting from the "
                "one that connects to the OH group, in the direction of the "
                "closest carbon atom with an alkyl attached to it. As such, "
                "this molecule is called 2,4,6-trinitrophenol, or, more "
                "explicitly, 2,4,6-trinitro-1-phenol or 2-hydroxy-1,3,5-"
                "trinitrobenzene."
            ),
            "O=[N+]([O-])c1cc(cc([N+]([O-])=O)c1O)[N+]([O-])=O",
            0,
        ),
        (
            6,
            "2-phenoxy-5-hydroxy-4-methylheptane;2-phenoxy-4-methyl-5-"
            "hydroxyheptane",
            (
                "The main chain is a heptane with a phenoxy group on the "
                "second carbon, a hydroxy group on the fifth carbon, and a "
                "methyl group on the fourth carbon: Thus, it can be named as "
                "2-phenoxy-5-hydroxy-4-methylheptane, or, equivalently, "
                "2-phenoxy-4-methyl-5-hydroxyheptane."
            ),
            "CC(OC1:C:C:C:C:C1)CC(C)C(O)CC",
            0,
        ),
        (
            7,
            "3,6-dimethyl-1,4-dioxane-2,5-dione;lactide;dilactide;"
            "(3S,6S)-3,6-dimethyl-1,4-dioxane-2,5-dione",
            (
                "This is a cyclic diester with oxygen atoms on positions 2 "
                "and 5, methyl groups on carbons 3 and 6, and a pair of "
                "oxygen atoms on carbons 1 and 4. Therefore, this molecule "
                "can be named 3,6-dimethyl-1,4-dioxane-2,5-dione; however, it "
                "is more often called lactide or dilactide."
            ),
            "CC1C(=O)OC(C(=O)O1)C",
            0,
        ),
        (
            8,
            "prolyl-leucyl-glycinamide;melanostatin",
            (
                "This molecule is a tripeptid that consists of proline, "
                "leucine, and glycine, so its name is prolyl-leucyl-"
                "glycinamide. It is commonly known as melanostatin."
            ),
            "NC(CNC([C@H](CC(C)C)NC([C@@H]1CCCN1)=O)=O)=O",
            0,
        ),
    ]

    cur.executemany("INSERT INTO molecule VALUES (?,?,?,?,?);", data)

    # Closing remarks
    con.commit()
    cur.close()
    con.close()
