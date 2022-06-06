# MoleWhat menu page
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
"""MoleWhat start page
Provides a guizero function that creates the start, allowing the user to begin
to play.
"""

import menu
import sqlite3
from guizero import App, Box, Text, Drawing, event


def get_progress():
    """Connects to database and returns the number of completed exercises."""

    con = sqlite3.connect("../database/Molecules.db")
    cur = con.cursor()
    cur.execute("SELECT * FROM molecule WHERE completed = ?", (1,))
    result = len(cur.fetchall())
    cur.close()
    con.close()

    return round(100 * result / 8)


def start(obj):
    """Clears the window and shows the level menu."""

    # Get App widget from input object
    if type(obj) == event.EventData:
        app = obj.widget.master.master
    else:
        app = obj

    # Clear window
    while len(app.children) != 0:
        app.children[0].destroy()

    # Top box
    Box(app, align="top", width="fill", height=65)

    # Box with progress
    a_box = Box(app, align="top", width="fill", height=170)
    a_box.bg = "white"
    Text(
        a_box,
        f"You have solved {get_progress()}% of the exercises",
        size=20,
        align="top",
        height=3,
    )

    # Box with button
    b_box = Box(app, align="top", width="fill", height=240)
    b_box.bg = "white"
    button = Drawing(b_box, width=240, height=240, align="top")
    button.triangle(60, 10, 60, 150, 200, 80, color="#8AC926")
    button.tk.config(cursor="hand1")
    button.when_clicked = menu.menu

    # Box with name
    c_box = Box(app, align="top", height=110, width="fill")
    c_box.bg = "white"
    Text(
        c_box,
        text="MoleWhat",
        color="black",
        bg="white",
        size=35,
        align="top",
    )
    #
    # for y in [0, 0]:
    #     for x in [0, 2, 3, 4, 5, 6]:
    #         button = Box(a_box, grid=[x, y], width=155, height=155)
    #         button = Drawing(a_box, grid=[3, 2], width=250, height=250)
    #         button.triangle(60, 10, 60, 150, 200, 80, color="#8AC926")
    #
    #         # Change cursor to a "hand" on mouse hover
    #         button.tk.config(cursor="hand1")
    #
    #         b_box = Box(a_box, grid=[3, 3], height=45, width=300)
    #         Text(
    #             b_box,
    #             text="MoleWhat",
    #             color="black",
    #             bg="white",
    #             size=35,
    #             width=100,
    #             height=100,
    #         )
    #
    #         # Call menu function if home button is clicked
    #         button.when_clicked = menu.menu


if __name__ == "__main__":
    """Driver code."""

    # Application initialization
    app = App(title="MoleWhat", height=650, width=900)
    app.tk.resizable(False, False)

    # Set background color to blue (#1982C4)
    app.bg = "#1982C4"

    # Run start page
    start(app)
    app.display()
