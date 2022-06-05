# MoleWhat writing type question page
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
"""MoleWhat writing_question
Provides a guizero function that creates a level type that displays a certain
molecule whose name must be typed by the user.
"""

import sqlite3
import start
from rdkit.Chem.Draw import MolToImage
from rdkit import Chem
from guizero import App, Box, Text, Picture, TextBox, Drawing
from tkinter import END

# Global variable to store level number
LEVEL_ID = 0


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
    cur.close()
    con.close()


def check_answer(event_data):
    """Checks the answer of the user's input, and shows a pop-up if the answer
    is correct or incorrect. If the answer is correct, the "completed" column
    for the given challenge number is switched to 1. In case it was incorrect,
    the pop-up will show an explanation."""

    # Get the input text when the user presses ENTER
    if len(event_data.key) != 0:
        if ord(event_data.key) == 13:
            answer = event_data.widget.value.lower()

            # Disable the textbox to avoid further modifications
            event_data.widget.disable()

            # Get App object from event_data
            app = event_data.widget
            while type(app) != App:
                app = app.master

            # Fetch level data from the database
            correct, explanation = lookup(LEVEL_ID)[1:3]
            correct = correct.split(";")

            # If the answer is correct, notify user and set "completed" to 1
            if answer in correct:
                set_as_solved(LEVEL_ID)
                app.info("Congratulations", "You solved this exercise.")
            # If not, show the user an explanation
            else:
                app.warn("Wrong answer", explanation)


def writing_question(event_data):
    """Generate a level where the user must type the name of a given
    molecule."""

    global LEVEL_ID
    # Fetch level data from the database
    LEVEL_ID, name, info, smiles, solved = lookup(int(event_data.widget.value))

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

    # Add an empty box as a boundary of the instructions text
    Box(app, align="top", width="fill", height=35)

    # Create the instructions box
    top_box = Box(app, align="top", width="fill", height=95)
    Text(
        top_box,
        text=f"Level {LEVEL_ID}: Type the name of this molecule",
        size=30,
        color="black",
    )

    # Add the image of the molecule
    mol_top = Box(app, align="top", width=650, height=15)
    mol_top.bg = "white"
    molecule = MolToImage(Chem.MolFromSmiles(smiles))
    mol_box = Box(app, align="top", width=650, height=300)
    mol_box.bg = "white"
    Picture(mol_box, image=molecule, width=300, height=300)
    mol_bottom = Box(app, align="top", width=650, height=15)
    mol_bottom.bg = "white"

    # "Home" button container box
    home_box = Box(app, align="bottom", height=130, width=800)

    # Draw the home button
    home = Drawing(home_box, align="right", width=100, height=100)
    home.oval(5, 5, 95, 95, color="#1982C4")
    home.triangle(50, 20, 24, 45, 76, 45, color="#8AC926")
    home.rectangle(35, 45, 65, 70, color="#FFCA3A")
    home.tk.config(cursor="hand1")
    home.when_clicked = start.start

    # Place the answer box for the user to type in
    Box(app, align="top", width="fill", height=20)
    answer_box = TextBox(app, text="Type here and press ENTER", width=50)
    answer_box.bg = "white"
    answer_box.text_size = 15
    answer_box.tk.selection_range(0, END)
    answer_box.focus()
    answer_box.when_key_pressed = check_answer
