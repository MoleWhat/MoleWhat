# MoleWhat multiple choice question page
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
"""MoleWhat multiple choice page
Provides a couple of molecule drawings and the user has to choose the one
asked by the program.
"""

import sqlite3
import start
from rdkit.Chem.Draw import MolToImage
from rdkit import Chem
from guizero import App, Box, Text, Picture, Drawing
from random import choice, shuffle

# Global variable to store level number and correct molecule id
LEVEL_ID = 0
CORRECT_ID = None


def do_nothing():
    """A function that does nothing..."""

    return


def lookup(number):
    """Looks up the database for the given challenge number and returns a tuple
    of the row's content."""

    con = sqlite3.connect("../database/Molecules.db")
    cur = con.cursor()
    cur.execute("SELECT * FROM molecule WHERE id = ?", (number,))
    result = cur.fetchall()[0]
    cur.close()
    con.close()

    return result


def set_as_solved(number):
    """Set a level as solved in the database."""

    con = sqlite3.connect("../database/Molecules.db")
    cur = con.cursor()
    cur.execute("UPDATE molecule SET completed = 1 WHERE id = ?", (number,))
    con.commit()
    cur.close()
    con.close()


def check_answer(event_data):
    """Checks the answer of the user's input, and shows a pop-up if the answer
    is correct or incorrect. If the answer is correct, the "completed" column
    for the given challenge number is switched to 1. In case it was incorrect,
    the pop-up will show an explanation."""

    # Get the id of the image chosen
    answer = id(event_data.widget.value)

    # Disable options so they cannot be interacted with
    for obj in event_data.widget.master.children:
        obj.when_clicked = do_nothing
        obj.tk.config(cursor="")

    # Get App object from event_data
    app = event_data.widget
    while type(app) != App:
        app = app.master

    if answer == CORRECT_ID:
        app.info("Congratulations", "You solved this exercise.")
        set_as_solved(LEVEL_ID)
    else:
        app.warn("Wrong answer", "Please, try again.")


def mc_question(event_data):
    """Creates a page with two options of molecules where the user must choose
    the one asked by the program."""

    global LEVEL_ID, CORRECT_ID
    # Fetch level data from the database
    LEVEL_ID, name, info, corr, solved = lookup(int(event_data.widget.value))

    # Randomly choose one of the names to show on-screen
    showed_name = choice(name.split(";"))

    # Get App object from event_data
    app = event_data.widget
    while type(app) != App:
        app = app.master

    # Clear window
    while len(app.children) != 0:
        app.children[0].destroy()

    # Set background color according to level id
    colors = [
        "#8AC926",
        "#8AC926",
        "#FFCA3A",
        "#FFCA3A",
        "#FF595E",
        "#FF595E",
        "#6A4C93",
        "#6A4C93",
    ]
    app.bg = colors[LEVEL_ID - 1]

    # "Home" button container box
    home_box = Box(app, align="bottom", height=130, width=800)

    # Draw the home button
    home = Drawing(home_box, align="right", width=100, height=100)
    home.oval(5, 5, 95, 95, color="#1982C4")
    home.triangle(50, 20, 24, 45, 76, 45, color="#8AC926")
    home.rectangle(35, 45, 65, 70, color="#FFCA3A")
    home.tk.config(cursor="hand1")
    home.when_clicked = start.start

    # Add an empty top box to set the top boundary for the text
    Box(app, align="top", width="fill", height=35)

    # Instructions container box
    top_box = Box(app, align="top", width="fill", height=95)
    Text(
        top_box,
        text=f"From the options below choose the\n{showed_name}",
        size=30,
        color="black",
    )

    # Get all molecules from the database that are not the correct answer
    # Randomly choose one of them
    con = sqlite3.connect("../database/Molecules.db")
    cur = con.cursor()
    cur.execute("SELECT code_smiles FROM molecule WHERE id != ?", (LEVEL_ID,))
    incorr = choice(cur.fetchall())[0]
    cur.close()
    con.close()

    # Molecule options container box
    opt = Box(app, align="top", width=650, height=300)
    options_order = [corr, incorr]
    shuffle(options_order)
    mol1 = MolToImage(Chem.MolFromSmiles(options_order[0]))
    mol2 = MolToImage(Chem.MolFromSmiles(options_order[1]))
    pic1 = Picture(opt, image=mol1, width=300, height=300, align="left")
    pic2 = Picture(opt, image=mol2, width=300, height=300, align="right")
    pic1.tk.config(cursor="hand1")
    pic2.tk.config(cursor="hand1")

    # Store id of the correct answer in CORRECT global variable
    CORRECT_ID = id(mol1) if options_order.index(corr) == 0 else id(mol2)

    # Check answer on mouse click
    pic1.when_clicked = check_answer
    pic2.when_clicked = check_answer
