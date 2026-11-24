"""Show the oscilloscope gain settings stored in one TRC descriptor."""

import argparse
from pathlib import Path
import struct


def read_vertical_gain(path):
    data = path.read_bytes()
    base = data.index(b"WAVEDESC")
    if len(data) - base < 344:
        raise ValueError("Truncated descriptor")
    order = struct.unpack_from("<H", data, base + 34)[0]
    if order not in (0, 1):
        raise ValueError("Invalid COMM_ORDER")
    byte_order = "<" if order == 1 else ">"
    fixed = struct.unpack_from(byte_order + "h", data, base + 332)[0]
    scale = (1, 2, 5)[fixed % 3] * 10 ** (fixed // 3 - 6)
    return {
        "Fixed scale [mV/div]": scale * 1000,
        "Vertical vernier": struct.unpack_from(byte_order + "f", data, base + 336)[0],
        "Probe attenuation": struct.unpack_from(byte_order + "f", data, base + 328)[0],
        "Vertical offset [V]": struct.unpack_from(byte_order + "f", data, base + 160)[0],
        "Acquisition vertical offset [V]": struct.unpack_from(byte_order + "f", data, base + 340)[0],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace", type=Path)
    arguments = parser.parse_args()
    for label, value in read_vertical_gain(arguments.trace).items():
        print(f"{label}: {value:g}")
