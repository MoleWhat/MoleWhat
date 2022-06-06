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

    print(event_data.tk_event.x)


def mc_question(event_data):
    """Creates a page with two options of molecules where the user must choose
    the one asked by the program."""

    global LEVEL_ID
    # Fetch level data from the database
    LEVEL_ID, name, info, smiles, solved = lookup(int(event_data.widget.value))

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
    incorrect = choice(cur.fetchall())[0]
    cur.close()
    con.close()

    # Molecule options container box
    opt = Box(app, align="top", width=650, height=300, layout="grid")
    options_order = [smiles, incorrect]
    shuffle(options_order)
    mol1 = MolToImage(Chem.MolFromSmiles(options_order[0]))
    mol2 = MolToImage(Chem.MolFromSmiles(options_order[1]))
    pic1 = Picture(opt, image=mol1, width=300, height=300, grid=[0, 0])
    Box(opt, width=50, height=300, grid=[1, 0])
    pic2 = Picture(opt, image=mol2, width=300, height=300, grid=[2, 0])
    pic1.tk.config(cursor="hand1")
    pic2.tk.config(cursor="hand1")


# def mc_question(app, selected_level):
#     # "Home" button container box
#     home_box = Box(app, align="bottom", height=130, width=800)
#
#     # Draw the home button
#     home = Drawing(home_box, align="right", width=100, height=100)
#     home.oval(5, 5, 95, 95, color="#1982C4")
#     home.triangle(50, 20, 24, 45, 76, 45, color="#8AC926")
#     home.rectangle(35, 45, 65, 70, color="#FFCA3A")
#     home.tk.config(cursor="hand1")
#     home.when_clicked = start.start
#
#     # Set the background color according to the level
#
#     if selected_level <= 2:
#         app.bg = "#8AC926"
#     elif selected_level <= 4:
#         app.bg = "#FFCA3A"
#     elif selected_level <= 6:
#         app.bg = "#FF595E"
#     elif selected_level <= 6:
#         app.bg = "#6A4C93"
#     else:
#         app.bg = "#6A4C93"
#
#     # Add an empty top box to set the top boundary for the text
#     Box(app, align="top", width="fill", height=40)
#
#     # Steps to obtain info from the database.
#     conexion = sqlite3.connect("../database/Molecules.db")
#     cursor = conexion.cursor()
#     cursor.execute("SELECT * FROM molecule")
#
#     # Get all the molecules and put them in a list
#     molecules = cursor.fetchall()
#     conexion.close()
#
#     # Instructions container box.
#     text_box = Box(app)
#     molecule_name = []
#     for i in molecules[selected_level - 1][1]:
#         if str(i) == ";":
#             break
#         else:
#             molecule_name.append(i)
#
#     msg = "Level " + str(selected_level) + ": From the options below choose
# the "
#     Text(text_box, text=msg, size=30)
#     name_box = Box(app)
#     Text(name_box, text=molecule_name, size=20)
#
#     # Display the options
#     levels = [0, 1, 2, 3, 4, 5, 6, 7]
#
#     el = Box(app, align="left", height=10, width=100)
#     er = Box(app, align="right", height=10, width=100)
#     a_box = Box(app, align="left", height=300, width=300)
#     b_box = Box(app, align="right", height=300, width=300)
#     boxes = [a_box, b_box]
#     m_a = Chem.MolFromSmiles(molecules[selected_level - 1][3])
#     img = Chem.Draw.MolToImage(m_a)
#     correct = choice(boxes)
#     picture = Picture(correct, image=img)
#     for i in boxes:
#         if i != correct:
#             other = i
#     # Avoid displaying the same molecule twice
#     a = selected_level - 1
#     while a == selected_level - 1:
#         a = choice(levels)
#
#     m_b = Chem.MolFromSmiles(molecules[a][3])
#     img_b = Chem.Draw.MolToImage(m_b)
#     picture_b = Picture(other, image=img_b)
#
#     # Verify if the given answer was correct
#     # =======================================================================
#     # correct.tk.config(cursor="hand1")
#     # correct.when_clicked = app.info
#     # other.tk.config(cursor="hand1")
#     # other.when_clicked = app.warn
#     # =======================================================================
#
#     # picture.tk.config(cursor="hand1")
#     # picture.when_clicked = correct_answ
#     # picture_b.tk.config(cursor="hand1")
#     # picture_b.when_clicked = app.warn()
#
#
# if __name__ == "__main__":
#     """Testing."""
#
#     app = App(title="MoleWhat", width=900, height=650)
#     app.tk.resizable(False, False)
#     mc_question(app, 4)
#     app.display()
