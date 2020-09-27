import colorsys
import itertools
import matplotlib.colors as cl
import matplotlib.cm as cm


def hls_to_hex(hue, luminance, saturation):
    """Return (R,G,B) equivalent of a hue/staturation/value color."""
    return cl.rgb2hex(colorsys.hls_to_rgb(hue, luminance, saturation))


def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return "#%02x%02x%02x" % (int(red), int(green), int(blue))


class ColorsMixin:
    def autocolor_quote_sources(
        self,
        hues=(0.635, 0.047, 0.117),
        saturations=(0.9, 0.7, 0.5, 0.3),
        min_lum=0.2,
        max_lum=0.8,
    ):
        """Auto-add a `_report_color` field to the sources in in quote.sources.

        Sources at the same depth share the same luminance.
        """

        colors = itertools.cycle(
            [
                rgb_to_hex(*[255 * e ** 0.4 for e in cm.Paired(0.13 * i % 1.0)][:3])
                for i in range(30)
            ]
        )
        for _name, source in sorted(self.sources.items()):
            color = next(colors)
            source._report_color = color
