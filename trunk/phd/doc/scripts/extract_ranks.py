import xml.etree.ElementTree as ET

root = ET.parse(file('pca_rank.xml','r')).getroot()
pcanode = root.find('pcaindex').find('k-eff-correction')
for v,var in enumerate(pcanode):
  print var.find('dim').text,'&',var.find('index').text,'\\\\'
  if float(var.find('index').text) < 1e-6: break
