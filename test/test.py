from guizero import *

def test():
    app.bg = 'red'

app = App(title="Main window")

test = PushButton(app, text="Open", command=test)

app.display()
