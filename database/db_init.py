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
database consists of two tables: one that stores molecule information, and
another one that stores challenge data.
"""
import sqlite3
from pathlib import Path

# Create database if it doesn't already exists
if Path("Molecules.db").is_file() is False:

    # Connection to database
    con = sqlite3.connect("Molecules.db")
    cur = con.cursor()

    # Molecule table creation command
    MOL_INIT = (
        "CREATE TABLE IF NOT EXISTS molecule ("
        "name_iupac TEXT PRIMARY KEY, "
        "description TEXT, "
        "code_smiles TEXT, "
        "difficulty INTEGER"
        ");"
    )

    # Challenge table creation command
    CHA_INIT = (
        "CREATE TABLE IF NOT EXISTS challenge ("
        "id INTEGER PRIMARY KEY, "
        "completed INTEGER, "
        "explanation TEXT, "
        "FOREIGN KEY (id) REFERENCES molecule (name_iupac)"
        ");"
    )

    # Creation of molecule and challenge tables
    cur.execute(MOL_INIT)
    cur.execute(CHA_INIT)

    # Data insertion
    data = [
        ("Oxidane,Water", "Inorganic chemical compound", "O", "1"),
        (
            "Diatomic oxygen, Oxygen",
            "Colourless, odorless and tasteless gas",
            "O=O",
            "1",
        ),
        ("Ethanol", "Also known as ethyl alcohol", "CCO", "2"),
        ("Ethanoic acid", "Found in vinegar", "CC(O)=O", "2"),
        (
            "Pyridine",
            "Belongs to the family of aromatic heterocyclics",
            "c1cnccc1",
            "3",
        ),
        (
            "D-Galactose,Galactose",
            "The sugar found in milk",
            "C(C1C(C(C(C(O1)O)O)O)O)O",
            "3",
        ),
        (
            "Deae-cellulose,Celusose",
            "It is a biopolymer composed of glucose",
            "C(C1C(C(C(C(O1)OC2C(OC(C(C2O)O)O)CO)O)O)O )O",
            "4",
        ),
        (
            "Deoxynivalenol,Vomitoxin",
            "This mycotoxin occurs predominantly in grains such as wheat, barley, oats, rye and corn",
            "CC4=CC3OC1C(O)CC(C)(C12CO2)C3(CO)C(O)C4= Or",
            "4",
        ),
    ]

    cur.executemany("INSERT INTO molecule VALUES (?,?,?,?);", data)

    # Closing remarks
    con.commit()
    cur.close()
    con.close()
