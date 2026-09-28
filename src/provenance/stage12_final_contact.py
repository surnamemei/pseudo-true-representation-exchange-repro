from pathlib import Path
from PIL import Image,ImageDraw
for group in ['final_main','final_supp']:
 fs=list(Path('stage12_pdf_preview',group).glob('page-*.png'));w,h=420,560;sheet=Image.new('RGB',(3*w,((len(fs)+2)//3)*(h+22)),'#ddd');d=ImageDraw.Draw(sheet)
 for i,f in enumerate(fs):
  im=Image.open(f);im.thumbnail((w-10,h));x=i%3*w;y=i//3*(h+22);sheet.paste(im,(x+5,y+20));d.text((x+10,y+3),f'{group} {i+1}',fill='black')
 sheet.save(Path('stage12_pdf_preview',group+'_contact.png'))
