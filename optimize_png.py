import argparse
import base64
from io import BytesIO
import pathlib
from PIL import Image

def optimizePNG(filename, colors=256, binary=True):
    palette = Image.ADAPTIVE
    with Image.open(filename) as im:
        im = im.convert("P", palette=palette, colors=colors)
        out = BytesIO()
        im.save(out, format="PNG", optimize=True)
    if binary:
        return out.getvalue()
    return base64.b64encode(out.getvalue()).decode("utf-8")

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dirname", type = pathlib.Path)
    args = parser.parse_args()

    for filename in args.dirname.rglob("*.png"):
        print("Optimizing file: %s" % filename)
        res = optimizePNG(filename, colors = 32, binary = True)
        with open(filename, 'wb') as f:
            f.write(res)
