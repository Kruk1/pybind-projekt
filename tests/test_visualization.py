from __future__ import annotations

import scikit_build_example as m

import wave
import struct
def read_wave(filename: str) -> list[int]:
    with wave.open(filename) as f:
        data = f.readframes(f.getnframes())
        samples = struct.unpack('{n}h'.format(n=f.getnframes()*f.getnchannels()), data)
    return samples

m.visualization(read_wave("file.wav"))


def test_version():
    assert m.__version__ == "0.0.1"


def test_add():
    assert m.add(1, 2) == 3


def test_sub():
    assert m.subtract(1, 2) == -1
