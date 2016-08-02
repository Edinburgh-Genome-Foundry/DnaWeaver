"""Global settings for DnaWeaver

Examples
--------

>>> import dnaweaver.config as config
>>> fontawesome_path = config["fontawesome-ttf-path"]
>>> # To change that path globally for DnaWeaver:
>>> config["fontawesome-ttf-path"] = "some/new/path.ttf"

"""

SETTINGS = {
    "fontawesome-ttf-path": '/usr/share/fonts/truetype/fontawesome-webfont.ttf',
    "OpenSans-ttf-path": '/usr/share/fonts/truetype/OpenSans-Light.ttf'
}
