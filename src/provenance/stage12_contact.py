from PIL import Image,ImageOps,ImageDraw
from pathlib import Path
for group in ['main','supp']:
 fs=list(Path('stage12_pdf_preview',group).glob('page-*.png'));w,h=410,550;sheet=Image.new('RGB',(3*w,((len(fs)+2)//3)*(h+25)),'#dddddd');d=ImageDraw.Draw(sheet)
 for i,f in enumerate(fs):
  im=Image.open(f);im.thumbnail((w-10,h));x=(i%3)*w;y=(i//3)*(h+25);sheet.paste(im,(x+5,y+20));d.text((x+10,y+3),f'{group} page {i+1}',fill='black')
 sheet.save(Path('stage12_pdf_preview',group+'_contact.png'))
