import imageio
from os import listdir
from reportlab.platypus import SimpleDocTemplate, Image
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch, cm
from reportlab.lib.pagesizes import letter, landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

pdfmetrics.registerFont(TTFont("DejaVuSans", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"))

#create canvas
cw = 24*cm
ch = 20*cm
c = canvas.Canvas("slides.pdf", pagesize=(cw, ch))
c.setFont("DejaVuSans", size=20)

#loop through all iamges and add page for each figure
files = listdir("raw/")
files.sort()
for fp in files:
    print("Adding page for " + fp)
    #image position settings
    iw = 8*cm
    ih = 8*cm
    wm = 2*cm
    hm = 1*cm
    x1 = wm
    x2 = cw - (iw+wm)
    y1 = hm/2
    y2 = ch - (ih+2*hm)

    #draw images
    base_path = fp[:-4]
    c.drawImage("dog/"+base_path+"_dog.png", x1, y1, iw, ih)
    c.drawImage("raw/"+fp, x1, y2, iw, ih)
    c.drawImage("sobel/"+base_path+"_sobel.png", x2, y1, iw, ih)
    c.drawImage("mask/"+base_path+"_mask.png", x2, y2, iw, ih)

    #write titles and image subtitles
    c.setFontSize(size=20)
    c.drawCentredString(x=cw/2,y=ch-1*cm,text=fp)

    c.setFontSize(size=14)
    sx1 = x1+iw/2
    sx2 = x2+iw/2
    sy1 = y1+ih+0.25*cm
    sy2 = y2+ih+0.25*cm
    c.drawCentredString(x=sx1, y=sy1, text="Difference of Gaussians")
    c.drawCentredString(x=sx1, y=sy2, text="Original")
    c.drawCentredString(x=sx2, y=sy1, text="Sobel Edge DoG")
    c.drawCentredString(x=sx2, y=sy2, text="Edge Mask")
    c.showPage()

c.save()