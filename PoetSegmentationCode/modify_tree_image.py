from PIL import Image, ImageDraw, ImageFont
from colour import Color

# modifies an image of a tree (adds gradient and title)
def modify(fname, title, color_bar_min, color_bar_max, num_color_bar_labels=2):
    # load image of tree
    original = Image.open(fname)

    # fonts
    font_name = "Calibri.ttf"
    number_fnt = ImageFont.truetype(font_name, 15)
    number_text_w, number_text_h = ImageDraw.Draw(original).textsize("test",
        font=number_fnt)
    if title != None:
        title_fnt = ImageFont.truetype(font_name, 20)
        title_text_w, title_text_h = ImageDraw.Draw(original).textsize(title,
            font=title_fnt)

    # visualization parameters
    buf = 15
    bar_height = 50
    lines_per_color = 2
    bar_start = int(original.size[0] / 2) - 100 * lines_per_color / 2
    bar_top = original.size[1] + buf
    new_height = 0
    title_y = bar_top
    if color_bar_max != None:
        new_height += buf*2 + bar_height + number_text_h
        title_y += bar_height + number_text_h + buf
    if title != None:
        new_height += buf*2 + title_text_h
    bar_border = 2

    # resize original tree image
    img = Image.new(original.mode, (original.size[0], original.size[1] + new_height),
        (255, 255, 255, 255))
    img.paste(original, (0, 0, original.size[0], original.size[1]))
    draw = ImageDraw.Draw(img)
    colors = [str(c) for c in list(Color("blue").range_to(Color("red"), 102))[1:-1]]

    if color_bar_max != None:
        # draw color bar border
        draw.rectangle([(bar_start - bar_border, bar_top - bar_border),
            (bar_start + lines_per_color*100 + bar_border, bar_top + bar_height + bar_border)],
            fill='black')
        label_range = color_bar_max - color_bar_min
        # draw gradient line by line
        for c_ndx in range(len(colors)):
            # draw color for multiple pixels
            for i in range(lines_per_color):
                draw.line((bar_start + i, bar_top, bar_start + i, bar_top + bar_height),
                    fill=colors[c_ndx])
            # color bar labels
            if c_ndx % int(100 / num_color_bar_labels) == 0:
                draw.text((bar_start, bar_top + bar_height + bar_border),
                    str(round(color_bar_min+label_range*c_ndx/100, 2)),
                    fill=(0,0,0,255), font=number_fnt)
            bar_start += lines_per_color

        # draw last color bar label (100)
        draw.text((bar_start, bar_top + bar_height + bar_border), str(color_bar_max),
            fill=(0,0,0,255), font=number_fnt)

    if title != None:
        # draw title of image
        draw.text((original.size[0] / 2 - title_text_w / 2,
            title_y), title, fill=(0,0,0,255), font=title_fnt)

    # save file
    img.save(fname)