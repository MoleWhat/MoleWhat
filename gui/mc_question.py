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
Provides a series of molecule draws and the user has to choose the one that he 
has been asked for.
"""
from guizero import *
from random import choice
import sqlite3
#import menu
from rdkit import Chem
import rdkit.Chem.Draw
import start



# =============================================================================
# def test3_func():
#     # Test function: replace with start function
#     print("You want to go home?")
# =============================================================================
    
def correct_answ():
    app.info("Great!", "Your answer is correct")
    conexion = sqlite3.connect('../database/Molecules.db')
    cursor = conexion.cursor()
    cursor.execute('UPDATE molecule SET completed = 1 WHERE id = selected_level')
    conexion.close()
    
def wrong_answ():
    app.warn("Wrong answer", molecules[selected_level-1][2])
    
def mc_question(app, selected_level):
    # "Home" button container box
    home_box = Box(app, align="bottom", height=130, width=800)

    # Draw the home button
    home = Drawing(home_box, align="right", width=100, height=100)
    home.oval(5, 5, 95, 95, color="#1982C4")
    home.triangle(50, 20, 24, 45, 76, 45, color="#8AC926")
    home.rectangle(35, 45, 65, 70, color="#FFCA3A")
    home.tk.config(cursor="hand1")
    home.when_clicked = start.start
    
    
    # Set the background color according to the level
    
    if selected_level<=2:
        app.bg = "#8AC926"
    elif selected_level<=4:
        app.bg = "#FFCA3A"
    elif selected_level<=6:
        app.bg = "#FF595E"
    elif selected_level<=6:
        app.bg = "#6A4C93"
    else:
        app.bg = "#6A4C93"
        
    # Add an empty top box to set the top boundary for the text
    Box(app, align="top", width="fill", height=40)
    
    # Steps to obtain info from the database.
    conexion = sqlite3.connect('../database/Molecules.db')
    cursor = conexion.cursor()
    cursor.execute("SELECT * FROM molecule")

    # Get all the molecules and put them in a list
    molecules = cursor.fetchall()
    conexion.close()
    
    
    # Instructions container box.
    text_box = Box(app)
    molecule_name = []
    for i in molecules[selected_level-1][1]:
        if str(i) == ';':
            break
        else:
            molecule_name.append(i)
            
    msg = "Level "+str(selected_level)+ ": From the options below choose the "
    Text(text_box, text=msg, size=30)
    name_box = Box(app)
    Text(name_box, text=molecule_name, size=20)
    
    # Display the options
    levels = [0, 1, 2, 3, 4, 5, 6, 7]
    
    el =Box(app, align="left", height=10, width=100)
    er = Box(app, align="right", height=10, width=100)
    a_box = Box(app, align="left", height=300, width=300)
    b_box = Box(app, align="right", height=300, width=300)
    boxes = [a_box, b_box]
    m_a = Chem.MolFromSmiles( molecules[selected_level-1][3])
    img = Chem.Draw.MolToImage(m_a)
    correct = choice(boxes)
    picture = Picture(correct, image=img)
    for i in boxes:
        if i != correct:
            other= i
    # Avoid displaying the same molecule twice
    a = selected_level-1
    while a==selected_level-1:
        a = choice(levels)

    m_b = Chem.MolFromSmiles(molecules[a][3])
    img_b = Chem.Draw.MolToImage(m_b)
    picture_b = Picture(other, image=img_b)
    
    #Verify if the given answer was correct
    
    picture.tk.config(cursor="hand1")
    picture.when_clicked = correct_answ
    picture_b.tk.config(cursor="hand1")
    picture_b.when_clicked = app.warm()

if __name__ == "__main__":
    """Testing."""

    app = App(title="MoleWhat", width=900, height=650)
    app.tk.resizable(False, False)
    mc_question(app,4)
    app.display()


