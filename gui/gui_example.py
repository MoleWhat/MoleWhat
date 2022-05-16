from guizero import App, Text, PushButton, Box

app = App('Chemistry')
box = Box(app, layout='grid')
title = Text(app, '¿Qué molécula representa el tolueno?')
i = 0
for x in range(2):
    for y in range(2):
        button = PushButton(box, image=f'molecule{i}.png',
                            grid=[x, y], width=200, height=200)
        i += 1

app.display()
