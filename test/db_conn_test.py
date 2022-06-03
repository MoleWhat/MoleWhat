import sqlite3

con = sqlite3.connect("../database/Molecules.db")

cur = con.cursor()

cur.execute("SELECT * FROM molecule;")

rows = cur.fetchall()
for row in rows:
    print(row)

con.close()
