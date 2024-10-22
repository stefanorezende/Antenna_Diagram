from PIL import Image
import os

def addIMG_ant(pngFile, dir, imgAnt):
    pathFile = os.path.join(dir, pngFile)
    # Abra a imagem de fundo (a imagem na qual você deseja inserir a camada)
    diagPNG = Image.open(pathFile)

    # Abra a imagem da camada que você deseja inserir
     
    antPNG = Image.open('img_ant\\'+imgAnt)

    # Verifique se as dimensões das imagens são compatíveis, ajustando, se necessário
    if diagPNG.size != antPNG.size:
        antPNG = antPNG.resize(diagPNG.size)

    # Aplique a opacidade (50%)
    diagPNG = diagPNG.convert('RGB')
    antPNG = antPNG.convert('RGB')

    # antPNG.show()
    resultPNG = Image.blend(diagPNG, antPNG, alpha=0.2)

    # Salve a imagem resultante
    resultPNG.save('output\\'+ pngFile, format='PNG')
    # imagem_camada=imagem_camada.apply_transparency()
    # imagem_fundo.paste(imagem_camada, (0, 0), imagem_camada)
    # resultPNG.show()


if __name__ == "__main__":
    dir = 'results\\Sierp_4'
    antenna ='2A'
    tx = 'H'
    diag = 'H'

    fileName = antenna+tx+diag
    imgAnt = 'HH.png' if diag == 'H' else 'VV.png'
    # Lista todos os arquivos na pasta
    folderFiles = os.listdir(dir)

    # Filtra os arquivos que começam com "2AVH" e têm a extensão .png
    pngFiles = [arquivo for arquivo in folderFiles if arquivo.startswith(fileName) and arquivo.endswith('.png')]

    for pngFile in pngFiles:
        addIMG_ant(pngFile, dir, imgAnt)