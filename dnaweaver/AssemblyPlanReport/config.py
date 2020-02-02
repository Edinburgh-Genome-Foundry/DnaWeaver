"""Global settings for DnaWeaver

Examples
--------

>>> import dnaweaver.config as config
>>> fontawesome_path = config["fontawesome-ttf-path"]
>>> # To change that path globally for DnaWeaver:
>>> config["fontawesome-ttf-path"] = "some/new/path.ttf"

"""

import os

this_path = os.path.join(os.path.dirname(os.path.realpath(__file__)))
fonts_path = os.path.join(this_path, "assets", "fonts")
fontawesome_path = os.path.join(fonts_path, "fontawesome-webfont.ttf")
opensans_path = os.path.join(fonts_path, "OpenSans-Light.ttf")
SETTINGS = {
    "fontawesome-ttf-path": fontawesome_path,
    "OpenSans-ttf-path": opensans_path,
    "template_path": os.path.join(this_path, "assets", "templates"),
}
